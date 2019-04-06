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

#include "contraction.h"
#include "talshxx.hpp"

using namespace std;
using namespace chrono;


int main(int argc, char *argv[]) {

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
    int num_entries = 40; // How many entries you want to run
    int entries_left = num_entries;
    vector<MPI_Request> requests(world_size, MPI_REQUEST_NULL);

    // Input variables (from file)
    // Open file
    string filename("generate_input/input.txt");
    auto io = ifstream(filename);
    assert(io.good() && "Cannot open file.");
    // Gotten from the file.
    // Number of arguments per amplitude and number of amplitudes per line
    int num_args, num_amps; 
    // The first element should be the number of arguments per line
    io >> num_args;
    io >> num_amps;
    int batch_size = num_args-2;
    vector<s_type> amplitudes(batch_size*num_amps*(world_size-1));
    // Grab lines one by one
    string line;
    getline(io, line); // "Read" one line to start next time from the second.
    // Send num_args
    for (int p=1; p<world_size; ++p)
    {
      MPI_Send(vector<int>({num_args,num_amps}).data(), 2, MPI_INT, p, 0,
               MPI_COMM_WORLD);
    }
    // Send second line with metadata
    if (getline(io, line)) if (line.size() && line[0] != '#')
    {
      for (int p=1; p<world_size; ++p)
      {
        MPI_Send(line.c_str(), line.size(), MPI_CHAR, p, 0,
                 MPI_COMM_WORLD);
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
        }
        if (flag)
        {
          if (getline(io, line)) if (line.size() && line[0] != '#')
          {
            MPI_Send(line.c_str(), line.size(), MPI_CHAR, p, 0,
                     MPI_COMM_WORLD);
            --entries_left;
            // Time
            t1 = high_resolution_clock::now();
            duration<double> span = duration_cast<duration<double>>(t1 - t0);
            // Time
            cout << span.count() << "s. Sent batch to process " << p
                 << ". Entries left = " << entries_left << "\n";
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
      }
      // Send a dummy pointer to show that the we are done!
      MPI_Send(NULL, 0, MPI_CHAR, p, 0,
               MPI_COMM_WORLD);
    }

    // Close file
    io.close();

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

    unsigned long mem_size(12000000000);
    talsh::initialize(&mem_size);
    {
      Contraction contraction(local_line, num_args, num_amps);

      while (true)
      {
        // Start timer
        t0 = high_resolution_clock::now();
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
        contraction.contract(local_line);
        t1 = high_resolution_clock::now();
        duration<double> span = duration_cast<duration<double>>(t1 - t0);

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
