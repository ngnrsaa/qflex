#include <vector>
#include <string>
#include <cmath>
#include <complex>
#include <ctime>
#include <chrono>
#include <cassert>

#include <iostream>
#include <fstream>
#include <sstream>

#include <omp.h>
#include <sched.h>

#include <mpi.h>

#include <thread>
#include <stdlib.h>

#include "talshxx.hpp"
#ifdef _51q
#include "contraction_51q.h"
#endif
#ifdef _bris_60
#include "contraction_bris_60.h"
#endif
#ifdef _7x7x40
#include "contraction_7x7x40.h"
#endif
#ifdef _8x8x32
#include "contraction_8x8x32.h"
#endif

using namespace std;
using namespace chrono;


int main(int argc, char *argv[]) {

  // Read input parameters



  // Three things to take as arguments: input file, output file, number of
  // amplitudes.

  // MPI
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);

  // Get number of precesses
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the processes
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Get the name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);


  //////////////////////////////////////////////////////////

  // Timing variables.
  high_resolution_clock::time_point t0, t1;
  duration<double> time_span;
  t0 = high_resolution_clock::now();

  // Start here

  //////////////////////// Process 0
  if (rank==0)
  {
    int num_entries = 50; // How many entries you want to run
    int entries_left = num_entries;
    vector<MPI_Request> requests(world_size, MPI_REQUEST_NULL);
    vector<string> parameter_strings(world_size-1);
    vector<string> input_strings(world_size-1);

    // Output file
    string out_filename = "outputs/output_bris_60.txt";
    ofstream out_file(out_filename);

    // Input variables (from file)
    // Open file
    string in_filename("inputs/input_bris_60.txt");
    auto in_file = ifstream(in_filename);
    assert(in_file.good() && "Cannot open file.");
    // Gotten from the file.
    // Number of arguments per amplitude and number of amplitudes per line
    int num_args, num_amps; 
    // The first element should be the number of arguments per line
    in_file >> num_args;
    in_file >> num_amps;
    int batch_size = num_args-2;
    vector<s_type> amplitudes(batch_size*num_amps*(world_size-1));
    // Grab lines one by one
    string line;
    getline(in_file, line); // "Read" one line to start next time from the second.
    // Send num_args
    for (int p=1; p<world_size; ++p)
    {
      MPI_Send(vector<int>({num_args,num_amps}).data(), 2, MPI_INT, p, 0,
               MPI_COMM_WORLD);
    }
    // Send second line with metadata
    if (getline(in_file, line)) if (line.size() && line[0] != '#')
    {
      for (int p=1; p<world_size; ++p)
      {
        MPI_Send(line.c_str(), line.size(), MPI_CHAR, p, 0,
                 MPI_COMM_WORLD);
        parameter_strings[p-1] = line;
      }
    }

    // Send messages
    while (entries_left>0)
    {
      for (int p=1; p<world_size; ++p)
      {
        if (entries_left<=0) break;
        int flag = true;
        if (requests[p]!=MPI_REQUEST_NULL) // First time is flag=true;
        {
          MPI_Request_get_status(requests[p], &flag, MPI_STATUS_IGNORE);
          if (flag)
          {
            out_file << "\nProcess " << p << ": " << parameter_strings[p-1]
                     << "\n";
            out_file << input_strings[p-1] << "\n";
            for (int a=0; a<batch_size*num_amps; ++a)
            {
              out_file << amplitudes[batch_size*num_amps*(p-1)+a].real()
                       << " "<< amplitudes[batch_size*num_amps*(p-1)+a].imag()
                       << " ";
            }
            out_file << "\n";
            out_file << flush;
          }
        }
        if (flag)
        {
          if (getline(in_file, line)) if (line.size() && line[0] != '#')
          {
            MPI_Send(line.c_str(), line.size(), MPI_CHAR, p, 0,
                     MPI_COMM_WORLD);
            input_strings[p-1] = line;
            --entries_left;
            // Time
            t1 = high_resolution_clock::now();
            duration<double> span = duration_cast<duration<double>>(t1 - t0);
            // Time
            //cout << span.count() << "s. Sent batch to process " << p
            //     << ". Entries left = " << entries_left << "\n";
            MPI_Irecv(amplitudes.data()+batch_size*num_amps*(p-1),
                      batch_size*num_amps, MPI_COMPLEX, p, 0,
                      MPI_COMM_WORLD, &requests[p]);
          }
        }
      }
    }

    // Clean up
    for (int p=0; p<world_size; ++p)
    {
      int flag;
      MPI_Request_get_status(requests[p], &flag, MPI_STATUS_IGNORE);
      if (!flag)
      {
        MPI_Wait(&requests[p], MPI_STATUS_IGNORE);
        out_file << "\nProcess " << p << ": " << parameter_strings[p-1]
                 << "\n";
        out_file << input_strings[p-1] << "\n";
        for (int a=0; a<batch_size*num_amps; ++a)
        {
          out_file << amplitudes[batch_size*num_amps*(p-1)+a].real()
                   << " "<< amplitudes[batch_size*num_amps*(p-1)+a].imag()
                   << " ";
        }
        out_file << "\n";
        out_file << flush;
      }
      // Send a dummy pointer to show that the we are done!
      MPI_Send(NULL, 0, MPI_CHAR, p, 0,
               MPI_COMM_WORLD);
    }

    // Close files
    out_file.close();
    in_file.close();

    // Report total time
    t1 = high_resolution_clock::now();
    duration<double> span = duration_cast<duration<double>>(t1 - t0);
    cout << "Total time = " << span.count() << endl;
  }

  ///////////////////////////// Process > 0
  else
  {
    MPI_Status status;
    int length;
    string local_line;
    vector<s_type> amplitudes;
    double amplitude = 0.0;
    int num_args;
    int num_amps;

    // get num_args
    vector<int> num_args_amps(2);
    MPI_Recv(num_args_amps.data(), 2, MPI_INT, 0, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    num_args = num_args_amps[0];
    num_amps = num_args_amps[1];
    // Get first string
    MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_CHAR, &length);
    char * char_ptr = new char[length];
    MPI_Recv(char_ptr, length, MPI_INT, 0, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    local_line = string(char_ptr);
    delete[] char_ptr;
    char_ptr = nullptr;

    unsigned long mem_size(16000000000);
    talsh::initialize(&mem_size);
    {
      Contraction contraction(local_line, num_args, num_amps);

      bool load_circuit(true);
      t0 = high_resolution_clock::now();
      while (true)
      {
        // Start timer
        ///////////////// Get string
        MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_CHAR, &length);
        char_ptr = new char[length];
        if (length==0) break;
        MPI_Recv(char_ptr, length, MPI_INT, 0, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        local_line = string(char_ptr);
        delete[] char_ptr;
        char_ptr = nullptr;

        // Computation
        if (load_circuit) // load only once
        {
          contraction.load_circuit(local_line);
          load_circuit = false;
        }
        contraction.contract(local_line);
        t1 = high_resolution_clock::now();
        duration<double> span = duration_cast<duration<double>>(t1 - t0);
        t0 = high_resolution_clock::now();
        cout << "A round took " << span.count() << "s" << endl;

        ///////////////// Message back
        amplitudes = contraction.get_amplitudes();
        MPI_Send(amplitudes.data(), amplitudes.size(), MPI_COMPLEX,
                 0, 0, MPI_COMM_WORLD);
      }
    }
    talsh::shutdown();

  }

  // Finalize the MPI environment.
  MPI_Finalize();

  return 0;
} 
