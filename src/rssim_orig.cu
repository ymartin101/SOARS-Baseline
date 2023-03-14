/// rssim.cpp - Functions which perform the actual simulations
/// Marc Brooker, 30 May 2006
/// Edited by Yaaseen Martin, 02 September 2019

#include <cmath>
#include <limits>
#include <stdexcept>
#include <stdio.h>
#include "rssim.cuh"
#include "rsworld.cuh"
#include "rsradar.cuh"
#include "rsdebug.cuh"
#include "rsparameters.cuh"
#include "rsantenna.cuh"
#include "rstarget.cuh"
#include "rsresponse.cuh"
#include "rsnoise.cuh"

// Check for CUDA errors
#define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while (0)

using namespace rs;

namespace {

    /// Results of solving the bistatic radar equation and friends
    struct REResults {
        rsFloat power;
        rsFloat delay;
        rsFloat doppler;
        rsFloat phase;
        rsFloat noise_temperature;
    };

  /// Kernel to solve the radar equation and friends (doppler, phase, delay) for direct Tx-to-Rx transmission
  void SolveREDirect(Transmitter *trans, Receiver *recv, rsFloat time, rsFloat length, const RadarSignal *wave, REResults &results, rsFloat temp)
  {
      // Calculate the vectors to and from the transmitter
      Vec3 tpos = trans->GetPosition(time);
      Vec3 rpos = recv->GetPosition(time);
      SVec3 transvec = SVec3(tpos - rpos);
      SVec3 recvvec = SVec3(rpos - tpos);

      // Calculate the range
      rsFloat R = transvec.length;

      // Normalize transvec and recvvec for angle calculations
      transvec.length = 1;
      recvvec.length = 1;

      // If the two antennas are not in the same position, this can be calculated
      if (R > std::numeric_limits<rsFloat>::epsilon())
      {
          // Step 1: Calculate the delay
          results.delay = R/rsParameters::c();

          // Calculate the wavelength
          rsFloat Wl = rsParameters::c()/wave->GetCarrier();

          // Get the antenna gains
          rsFloat Gt = trans->GetGain(transvec, trans->GetRotation(time), Wl);
          rsFloat Gr = recv->GetGain(recvvec, recv->GetRotation(time+results.delay), Wl);

          // Step 2: Calculate the received power
          results.power = Gt*Gr*Wl*Wl/(4*M_PI);
          if (!recv->CheckFlag(Receiver::FLAG_NOPROPLOSS))
              results.power *= 1/(4*M_PI*pow(R, 2));

          // Step 3: Calculate the doppler shift (if one of the antennas is moving)
          Vec3 tpos_end = trans->GetPosition(time+length);
          Vec3 rpos_end = recv->GetPosition(time+length);
          Vec3 trpos_end = tpos_end - rpos_end;
          rsFloat R_end = trpos_end.Length();

          // Compute the Doppler shift
          rsFloat vdoppler = (R_end - R)/length;
          results.doppler = (rsParameters::c() + vdoppler)/(rsParameters::c() - vdoppler);

          // Receiver duals do not receive any direct transmissions
          // However, real receivers can receive direct transmissions from a transmitter dual
          if (trans->IsMultipathDual())
              results.power *= trans->MultipathDualFactor();

          // Step 4: Calculate phase shift
          results.phase = fmod(results.delay*2*M_PI*wave->GetCarrier(), 2*M_PI);

          // Step 5: Account for the external plus antenna noise temperature
          rsFloat recv_temp = wave->GetTemp() + temp;
          results.noise_temperature = recv_temp;
      }
  }

  // Function for direct Tx-to-Rx processing - model the pulse which is received directly by a receiver from a transmitter
  void AddDirect(World *world, Receiver *recv, Transmitter *trans, TransmitterPulse *signal, RadarSignal *wave, rsFloat temp, rsFloat pulse_length) {

      // If receiver and transmitter share the same antenna (monostatic), there cannot be a direct pulse between Tx and Rx
      if (trans->IsMonostatic() && (trans->GetAttached() == recv)){
          return;
      }

      // Get the pulse start and end times
      rsFloat start_time = signal->time;
      rsFloat end_time = start_time + pulse_length;

      // Calculate the number of interpolation points we need to add
      rsFloat sample_time = 1.0 / rsParameters::cw_sample_rate(); // Default CW sample rate used as 1000 Hz
      rsFloat point_count = std::ceil(pulse_length/sample_time);

      // Create the response
      Response *response = new Response(wave, trans);

      // Loop through and add interpolation points
      for (int i = 0; i < point_count; i++) {
          rsFloat stime = i*sample_time + start_time;
          REResults results;
          SolveREDirect(trans, recv, stime, sample_time, wave, results, temp);
          InterpPoint point(results.power, results.delay + stime, results.delay, results.doppler, results.phase, results.noise_temperature);
          response->AddInterpPoint(point);
      }

      // Add one more point at the end
      REResults results;
      SolveREDirect(trans, recv, end_time, sample_time, wave, results, temp);
      InterpPoint point(results.power, results.delay+end_time, results.delay, results.doppler, results.phase, results.noise_temperature);
      response->AddInterpPoint(point);

      // Add the response to the receiver
      recv->AddResponse(response);
  }
}

/// Solve the radar equation for a given set of parameters
__global__ void SolveRE(rsFloat length, REResults* results, rsFloat temp, int targsize, int numPoints, int pulses, Vec3* trpos_arr, Vec3* trpos_end_arr,
                        Vec3* repos_arr, Vec3* repos_end_arr, Vec3* tapos_start_arr, Vec3* tapos_end_arr,
                        rsFloat carrier, rsFloat c, rsFloat* rcs_arr, rsFloat* gt_arr, rsFloat* gr_arr, bool noprop,
                        rsFloat trans_dualfactor, rsFloat recv_dualfactor, rsFloat recv_temp)
{
    // Only use threads up until the end of (targsize*numPoints*pulses)
    int tid = int(threadIdx.x + (blockIdx.x * blockDim.x)); // Thread index
    int stride = blockDim.x * gridDim.x;  // Total number of threads spawned
    for (int i = tid; i < (targsize*numPoints*pulses); i += stride){ // Will only iterate again if (targsize*numPoints*pulses) > (MaxBlocks*MaxThreads); very unlikely

        // Get the positions of simulation objects
        int time_idx = floor(float(i/(targsize)));
        Vec3 tapos = tapos_start_arr[i];
        SVec3 transvec = SVec3(tapos - trpos_arr[time_idx]);
        SVec3 recvvec = SVec3(tapos - repos_arr[time_idx]);

        // Calculate the distances
        rsFloat Rt = transvec.length;
        rsFloat Rr = recvvec.length;

        // Step 1: Calculate the delay (in seconds) experienced by the pulse; see "Delay Equation" in doc/equations/equations.tex
        results[i].delay = (Rt + Rr)/c;

        // Get the RCS
        rsFloat rcs = rcs_arr[i];

        // Calculate the wavelength
        rsFloat Wl = c/carrier;

        // Get the system antenna gains (which include loss factors) from arrays
        rsFloat Gt = gt_arr[i];
        rsFloat Gr = gr_arr[i];

        // Step 2: Calculate the received power using the narrowband bistatic radar equation; see "Bistatic Narrowband radar equation" in doc/equations/equations.tex
        results[i].power = Gt*Gr*rcs/(4*M_PI);
        if (!noprop)
            results[i].power *= (Wl*Wl)/(pow(4*M_PI, 2)*Rt*Rt*Rr*Rr);

        // Multiply by the loss factor; neither, one, or both of Tx/Rx can "observe" a reflection
        results[i].power *= trans_dualfactor; // If the transmitter is a dual, account for one reflection
        results[i].power *= recv_dualfactor; // If the receiver is a dual, account for one reflection

        // Step 3: Calculate phase shift; see "Phase Delay Equation" in doc/equations/equations.tex
        results[i].phase = fmod(results.delay*2*M_PI*carrier, 2*M_PI);

        // Step 4: Calculate doppler shift
        
        Vec3 tapos_end = tapos_end_arr[i];  // Need to project onto bistatic plane
        SVec3 transvec_end = SVec3(tapos_end - trpos_end_arr[time_idx]);
        SVec3 recvvec_end = SVec3(tapos_end - repos_end_arr[time_idx]);
        rsFloat Rt_end = transvec_end.length;
        rsFloat Rr_end = recvvec_end.length;

        // Doppler shift equation
        rsFloat V_t = (Rt_end - Rt)/length;
        rsFloat V_r = (Rr_end - Rr)/length;
        results[i].doppler = (std::sqrt((1 + V_r/c)/(1 - V_r/c))*std::sqrt((1 + V_t/c)/(1 - V_t/c))); // Calculates Fr/Fc

        // Step 5: Account for the external + antenna noise temperature
        results[i].noise_temperature = recv_temp + temp;

        rsFloat dopp = results[i].doppler)*1.35e9;
        printf("Doppler: %.15lf, x0: %.15lf, y0: %.15lf, z0: %.15lf\n", dopp, tapos_start_arr[i].x, tapos_start_arr[i].y, tapos_start_arr[i].z);
        printf("Doppler: %.15lf, x1: %.15lf, y1: %.15lf, z1: %.15lf\n", dopp, tapos_end_arr[i].x, tapos_end_arr[i].y, tapos_end_arr[i].z);
    }
}

using namespace std; // For vector usage etc.

/// Simulate a transmitter-receiver pair with a pulsed transmission
void rs::SimulatePair(Transmitter *trans, Receiver *recv, World *world, int MaxThreads, int MaxBlocks)
{
    // Set up general variables
    int pulses = trans->GetPulseCount();  // Number of pulses to transmit
    rsFloat c = rsParameters::c();  // Speed of propagation
    bool noprop = recv->CheckFlag(Receiver::FLAG_NOPROPLOSS); // Account for propagation loss?
    int targsize = (world->targets).size(); // Number of targets
    rsFloat sim_endtime = rsParameters::end_time();  // Endtime of the full simulation

    // Account for multipath of trans and recv
    rsFloat trans_dualfactor; // Transmitter multipath dual factor
    rsFloat recv_dualfactor; // Receiver multipath dual factor
    if (trans->IsMultipathDual())
        trans_dualfactor = trans->MultipathDualFactor(); // If the transmitter is a dual, account for one reflection
    else
        trans_dualfactor = 1;   // Will multiply this with overall power later
    if (recv->IsMultipathDual())
        recv_dualfactor = recv->MultipathDualFactor(); // If the receiver is a dual, account for one reflection
    else
        recv_dualfactor = 1;    // Will multiply this with overall power later

    // Get first pulse signal and set up further variables
    TransmitterPulse* signal = new TransmitterPulse();
    trans->GetPulse(signal, 0); // Pulse signal
    RadarSignal *wave = signal->wave;  // Pulse wave
    rsFloat pulse_length = wave->GetLength(); // Pulse length
    rsFloat carrier = wave->GetCarrier(); // Carrier frequency
    rsFloat Wl = c/carrier; // Wavelength
    rsFloat recv_temp = wave->GetTemp();  // Receiver temperature
    rsFloat temp = recv->GetNoiseTemperature(); // Internal noise
    rsFloat sample_time = 1.0 / rsParameters::cw_sample_rate(); // Default CW sample rate used as 1000 Hz
    int point_count = ceil(pulse_length / sample_time); // Number of interpolation points we need to add
    int numPoints = point_count + 1;  // Number of Points

    // Set up arrays for use on kernel
    Target** targ_arr = (world->targets).data(); // Targets array
    vector<Vec3> tapos_start_arr(targsize*numPoints*pulses); // Target present positions vector
    vector<Vec3> tapos_end_arr(targsize*numPoints*pulses); // Target future positions vector
    vector<rsFloat> rcs_arr(targsize*numPoints*pulses); // Target RCS vector
    vector<rsFloat> gt_arr(targsize*numPoints*pulses); // Transmitter gain vector
    vector<rsFloat> gr_arr(targsize*numPoints*pulses); // Receiver gain vector
    vector<REResults> results_arr(targsize*numPoints*pulses); // Results struct vector
    vector<Vec3> trpos_arr(numPoints*pulses); // Tx present positions vector
    vector<Vec3> repos_arr(numPoints*pulses); // Tx present positions vector
    vector<Vec3> trpos_end_arr(numPoints*pulses); // Tx future positions vector
    vector<Vec3> repos_end_arr(numPoints*pulses); // Rx future positions vector
    vector<rsFloat> start_time_arr(pulses); // start_time vector varying with pulse number
    vector<rsFloat> time_vec(numPoints*pulses); // Time vector

    // Create device copies of arrays for kernel
    Vec3 *d_tapos_start_arr;
    Vec3 *d_tapos_end_arr;
    Vec3 *d_trpos_arr;
    Vec3 *d_repos_arr;
    Vec3 *d_trpos_end_arr;
    Vec3 *d_repos_end_arr;
    rsFloat *d_rcs_arr;
    rsFloat *d_gt_arr;
    rsFloat *d_gr_arr;
    REResults *d_results_arr;

    // Allocate memory for device variable copies
    cudaMalloc((void **)&d_tapos_start_arr, sizeof(rs::Vec3)*targsize*numPoints*pulses);
    cudaMalloc((void **)&d_tapos_end_arr, sizeof(rs::Vec3)*targsize*numPoints*pulses);
    cudaMalloc((void **)&d_trpos_arr, sizeof(rs::Vec3)*numPoints*pulses);
    cudaMalloc((void **)&d_trpos_end_arr, sizeof(rs::Vec3)*numPoints*pulses);
    cudaMalloc((void **)&d_repos_arr, sizeof(rs::Vec3)*numPoints*pulses);
    cudaMalloc((void **)&d_repos_end_arr, sizeof(rs::Vec3)*numPoints*pulses);
    cudaMalloc((void **)&d_rcs_arr, sizeof(rsFloat)*targsize*numPoints*numPoints*pulses);
    cudaMalloc((void **)&d_gt_arr, sizeof(rsFloat)*targsize*numPoints*pulses);
    cudaMalloc((void **)&d_gr_arr, sizeof(rsFloat)*targsize*numPoints*pulses);
    cudaMalloc((void **)&d_results_arr, sizeof(REResults)*targsize*numPoints*pulses);
    cudaCheckErrors("Malloc fail");

    // Fill the arrays for kernel transfer
    for (int k = 0; k < pulses; k++){
        trans->GetPulse(signal, k); // Pulse signal
        start_time_arr[k] = signal->time;

        // Loop through Points while looping through pulses
        for (int n = 0; n < numPoints; n++){
            int nk_index = n + (numPoints*k);

            // Assign time values
            if (n == point_count)
                time_vec[nk_index] = start_time_arr[k] + pulse_length; // Final Point
            else
                time_vec[nk_index] = (n*sample_time) + start_time_arr[k]; // Time of the start of the sample

            // Assign Tx/Rx positions
            trpos_arr[nk_index] = trans->GetPosition(time_vec[nk_index]);
            repos_arr[nk_index] = recv->GetPosition(time_vec[nk_index]);
            trpos_end_arr[nk_index] = trans->GetPosition(time_vec[nk_index] + sample_time);
            repos_end_arr[nk_index] = recv->GetPosition(time_vec[nk_index] + sample_time);

            // Loop through targets while looping through Points and pulses
            for (int j = 0; j < targsize; j++){
                int jnk_index = j + (targsize*n) + (targsize*numPoints*k); // Overall index for (pulse, point, target)
                tapos_start_arr[jnk_index] = targ_arr[j]->GetPosition(time_vec[nk_index]);
                tapos_end_arr[jnk_index] = targ_arr[j]->GetPosition(time_vec[nk_index] + sample_time);
                SVec3 transvec = SVec3(tapos_start_arr[jnk_index] - trpos_arr[nk_index]);
                SVec3 recvvec = SVec3(tapos_start_arr[jnk_index] - repos_arr[nk_index]);
                rsFloat Rt = transvec.length;
                rsFloat Rr = recvvec.length;
                transvec.length = 1;  // Normalise
                recvvec.length = 1;   // Normalise
                rcs_arr[jnk_index] = targ_arr[j]->GetRCS(transvec, recvvec);
                gt_arr[jnk_index] = trans->GetGain(transvec, trans->GetRotation(time_vec[nk_index]), Wl);
                gr_arr[jnk_index] = recv->GetGain(recvvec, recv->GetRotation(((Rt + Rr)/c) + time_vec[nk_index]), Wl);
            }
        }
    }

    // Copy variables from host to device (needs to repeat for every new Point)
    cudaMemcpy(d_tapos_start_arr, tapos_start_arr.data(), sizeof(rs::Vec3)*targsize*numPoints*pulses, cudaMemcpyHostToDevice);
    cudaMemcpy(d_tapos_end_arr, tapos_end_arr.data(), sizeof(rs::Vec3)*targsize*numPoints*pulses, cudaMemcpyHostToDevice);
    cudaMemcpy(d_trpos_arr, trpos_arr.data(), sizeof(rs::Vec3)*numPoints*pulses, cudaMemcpyHostToDevice);
    cudaMemcpy(d_trpos_end_arr, trpos_end_arr.data(), sizeof(rs::Vec3)*numPoints*pulses, cudaMemcpyHostToDevice);
    cudaMemcpy(d_repos_arr, repos_arr.data(), sizeof(rs::Vec3)*numPoints*pulses, cudaMemcpyHostToDevice);
    cudaMemcpy(d_repos_end_arr, repos_end_arr.data(), sizeof(rs::Vec3)*numPoints*pulses, cudaMemcpyHostToDevice);
    cudaMemcpy(d_rcs_arr, rcs_arr.data(), sizeof(rsFloat)*targsize*numPoints*pulses, cudaMemcpyHostToDevice);
    cudaMemcpy(d_gt_arr, gt_arr.data(), sizeof(rsFloat)*targsize*numPoints*pulses, cudaMemcpyHostToDevice);
    cudaMemcpy(d_gr_arr, gr_arr.data(), sizeof(rsFloat)*targsize*numPoints*pulses, cudaMemcpyHostToDevice);
    cudaCheckErrors("Memory fail");

    // Call kernel
    if ((targsize*numPoints*pulses) <= MaxThreads){ // If there are fewer targets than the number of threads in one block (or an equal number of targets)
        SolveRE<<<1, (targsize*numPoints*pulses)>>>(sample_time, d_results_arr, temp, targsize, numPoints, pulses, d_trpos_arr, d_trpos_end_arr,
                                  d_repos_arr, d_repos_end_arr, d_tapos_start_arr, d_tapos_end_arr, carrier, c,
                                  d_rcs_arr, d_gt_arr, d_gr_arr, noprop, trans_dualfactor, recv_dualfactor, recv_temp);
    }
    else if ((targsize*numPoints*pulses) > (MaxThreads*MaxBlocks)){   // If there are more targets than the maximum number of parallel threads (across all blocks)
        SolveRE<<<MaxBlocks, MaxThreads>>>(sample_time, d_results_arr, temp, targsize, numPoints, pulses, d_trpos_arr, d_trpos_end_arr,
                                  d_repos_arr, d_repos_end_arr, d_tapos_start_arr, d_tapos_end_arr, carrier, c,
                                  d_rcs_arr, d_gt_arr, d_gr_arr, noprop, trans_dualfactor, recv_dualfactor, recv_temp);
    }
    else{ // If number of targets requires more than 1 block, but not all of them
        SolveRE<<<(((targsize*numPoints*pulses) + (MaxThreads - 1))/MaxThreads), MaxThreads>>>(sample_time, d_results_arr, temp, targsize, numPoints, pulses,
                                  d_trpos_arr, d_trpos_end_arr, d_repos_arr, d_repos_end_arr, d_tapos_start_arr, d_tapos_end_arr, carrier, c,
                                  d_rcs_arr, d_gt_arr, d_gr_arr, noprop, trans_dualfactor, recv_dualfactor, recv_temp);
    }
    cudaDeviceSynchronize();
    cudaCheckErrors("Kernel fail");
    cudaMemcpy(results_arr.data(), d_results_arr, sizeof(REResults)*targsize*numPoints*pulses, cudaMemcpyDeviceToHost); // Copy from GPU to host memory
    cudaCheckErrors("Memory fail");

    // Add all responses to receiver
    for (int k = 0; k < pulses; k++){
        for (int t = 0; t < targsize; t++){ // Iterate through all targets
            REResults results;  // To store results
            Response *response = new Response(wave, trans); // Create new response for every target
            for (int n = 0; n < point_count; n++){ // Iterate through the points
                int ptk_index = t + (targsize*n) + (targsize*numPoints*k); // Overall index for (pulse, point, target)
                results = results_arr[ptk_index];
                if (((n * sample_time) + start_time_arr[k] + results.delay) < sim_endtime){  // Only add the Point if it will be received within the simulation time
                    InterpPoint point(results.power, (n * sample_time) + start_time_arr[k] + results.delay, results.delay, results.doppler, results.phase, results.noise_temperature);
                    response->AddInterpPoint(point);  // Add the point to the response
                }
            }
            int ptk_index_end = t + (targsize*point_count) + (targsize*numPoints*k);
            results = results_arr[ptk_index_end];
            if ((start_time_arr[k] + pulse_length + results.delay) < sim_endtime){  // Only add the Point if it will be received within the simulation time
                InterpPoint point(results.power, start_time_arr[k] + pulse_length + results.delay, results.delay, results.doppler, results.phase, results.noise_temperature);
                response->AddInterpPoint(point);  // Add the final point to the response
            }

            // Add the response to the receiver for every target (only if there are Points recorded in the response)
            if (response->CountPoints() > 0)
                recv->AddResponse(response);
        }
    }

    // Check if direct pulses are being considered for this receiver
    if (!recv->CheckFlag(Receiver::FLAG_NODIRECT))
        AddDirect(world, recv, trans, signal, wave, temp, pulse_length);  // Add the direct pulses (if any)

    // Set new receiver noise temperature (antenna + external noise)
    recv->SetNoiseTemperature(wave->GetTemp() + temp); // "wave" remains the same after i = 0 in the for loop

    // Free up memory
    cudaFree(d_tapos_start_arr);
    cudaFree(d_tapos_end_arr);
    cudaFree(d_rcs_arr);
    cudaFree(d_gt_arr);
    cudaFree(d_gr_arr);
    cudaFree(d_results_arr);
    cudaFree(d_trpos_arr);
    cudaFree(d_trpos_end_arr);
    cudaFree(d_repos_arr);
    cudaFree(d_repos_end_arr);
    cudaCheckErrors("Delete fail");
}
