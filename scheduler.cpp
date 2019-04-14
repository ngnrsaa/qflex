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
#ifdef _11x11x24
#include "contraction_11x11x24.h"
#endif

using namespace std;
using namespace chrono;


int main(int argc, char *argv[]) {


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

  // Read input parameters
  if (argc<5) throw logic_error("ERROR: Not enough arguments.");
  string in_filename = string(argv[1]);
  string out_filename = string(argv[2]);
  int mem_size_GB = atoi(argv[3]);
  unsigned long mem_size = mem_size_GB * size_t(1073741824);
  int num_entries = atoi(argv[4]); // How many entries from input file to run


  // Timing variables.
  high_resolution_clock::time_point t0, t1;
  high_resolution_clock::time_point t0_init_total, t0_total, t1_total;
  duration<double> span_total;
  t0 = high_resolution_clock::now();

  // Start here

  //////////////////////// Process 0
  if (rank==0)
  {

    int entries_left = num_entries;
    vector<MPI_Request> requests(world_size, MPI_REQUEST_NULL);
    vector<string> parameter_strings(world_size);
    vector<string> input_strings(world_size);

    // Written?
    vector<bool> written(world_size, false);
    // Input lines
    vector<string> input_lines(num_entries);

    // Output file
    ofstream out_file(out_filename);
    out_file.precision(7);

    // Input variables (from file)
    // Open file
    auto in_file = ifstream(in_filename);
    assert(in_file.good() && "Cannot open file.");
    // Gotten from the file.
    // Number of arguments per amplitude and number of amplitudes per line
    int num_args, num_amps; 
    // The first element should be the number of arguments per line
    in_file >> num_args;
    in_file >> num_amps;
    int batch_size = num_args-2;

    // Grab lines one by one
    string line;
    getline(in_file, line); // "Read" one line to start next from the second.
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
        parameter_strings[p] = line;
      }
    }

    // Read all other lines
    t0 = high_resolution_clock::now();
    int line_idx = 0;
    while (getline(in_file, line) && line_idx<num_entries) if (line.size() && line[0] != '#')
    {
      input_lines[line_idx] = line;
      ++line_idx;
    }
    in_file.close();
    line_idx = 0;
    t1 = high_resolution_clock::now();
    duration<double> span = duration_cast<duration<double>>(t1 - t0);
    cout << "Reading took " << span.count() << endl;

    // Set up for messages and storage
    vector<s_type> amplitudes_out(batch_size*num_amps*num_entries);
    vector<int> tasks_out(num_entries);
    vector<double> peak_times_out(num_entries);
    vector<string> input_strings_out(num_entries);
    int idx_entry = 0;
    // amplitudes buffer. +1 is to receive the time_largest_contraction.
    vector<s_type> amplitudes_buffer((batch_size*num_amps+1)*world_size);

    // BEGIN TIMER
    ////// MEASURE BEGIN TIMER
    {
      time_t t = time(NULL);
      tm* timePtr = localtime(&t);
      cout << "Initializing started at: "
           << timePtr->tm_mon << "/"
           << timePtr->tm_mday << "/"
           << (timePtr->tm_year)+1900 << " "
           << timePtr->tm_hour << ":"
           << timePtr->tm_min << ":"
           << timePtr->tm_sec << "\n" << flush;
    }
    ////// END MEASURE BEGIN TIMER
    MPI_Barrier(MPI_COMM_WORLD);
    ////// MEASURE COMPUTATION BEGIN TIMER
    {
      time_t t = time(NULL);
      tm* timePtr = localtime(&t);
      cout << "Computation started at: "
           << timePtr->tm_mon << "/"
           << timePtr->tm_mday << "/"
           << (timePtr->tm_year)+1900 << " "
           << timePtr->tm_hour << ":"
           << timePtr->tm_min << ":"
           << timePtr->tm_sec << "\n" << flush;
      t0_init_total = high_resolution_clock::now();
    }
    ////// END MEASURE COMPUTATION BEGIN TIMER


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
            for (int a=0; a<batch_size*num_amps; ++a)
            {
              amplitudes_out[batch_size*num_amps*idx_entry+a] =
                amplitudes_buffer[(batch_size*num_amps+1)*p+a];
            }
            // Populate out vectors
            tasks_out[idx_entry] = p;
            input_strings_out[idx_entry] = input_strings[p];
            peak_times_out[idx_entry] =
              (double)amplitudes_buffer[(batch_size*num_amps+1)*(p+1)-1].real();
            ++idx_entry;
            
            written[p] = true;
          }
        }
        if (flag)
        {
          line = input_lines[line_idx];
          ++line_idx;
          {
            MPI_Send(line.c_str(), line.size(), MPI_CHAR, p, 0,
                     MPI_COMM_WORLD);
            
            written[p] = false;
            input_strings[p] = line;
            --entries_left;
            // Time
            t1 = high_resolution_clock::now();
            span = duration_cast<duration<double>>(t1 - t0);
            // Time
            //cout << span.count() << "s. Sent batch to process " << p
            //     << ". Entries left = " << entries_left << "\n";
            MPI_Irecv(amplitudes_buffer.data()+(batch_size*num_amps+1)*(p),
                      batch_size*num_amps+1, MPI_COMPLEX, p, 0,
                      MPI_COMM_WORLD, &requests[p]);
          }
        }
      }
    }

    // Clean up
    for (int p=1; p<world_size; ++p)
    {
      int flag;
      MPI_Request_get_status(requests[p], &flag, MPI_STATUS_IGNORE);
      if (!flag || !written[p])
      {
        MPI_Wait(&requests[p], MPI_STATUS_IGNORE);
        for (int a=0; a<batch_size*num_amps; ++a)
        {
          amplitudes_out[batch_size*num_amps*idx_entry+a] =
            amplitudes_buffer[(batch_size*num_amps+1)*p+a];
        }

        // Populate out vectors
        tasks_out[idx_entry] = p;
        input_strings_out[idx_entry] = input_strings[p];
        peak_times_out[idx_entry] =
          (double)amplitudes_buffer[(batch_size*num_amps+1)*(p+1)-1].real();
        ++idx_entry;

        written[p] = true;
      }
      // Send a dummy pointer to show that the we are done!
      MPI_Send(NULL, 0, MPI_CHAR, p, 0,
               MPI_COMM_WORLD);
    }

    // END TIMER
    MPI_Barrier(MPI_COMM_WORLD);
    {
      time_t t = time(NULL);
      tm* timePtr = localtime(&t);
      cout << "Computation ended at: "
           << timePtr->tm_mon << "/"
           << timePtr->tm_mday << "/"
           << (timePtr->tm_year)+1900 << " "
           << timePtr->tm_hour << ":"
           << timePtr->tm_min << ":"
           << timePtr->tm_sec << "\n" << flush;
    }
    t1_total = high_resolution_clock::now();
    span_total = duration_cast<duration<double>>(t1_total - t0_init_total);
    double total_time = span_total.count();
    cout << "Total time (s): " << total_time << endl;

    // Write everything to file
    t0 = high_resolution_clock::now();
    out_file << "\n---------- BEGIN OUTPUT ------------\n";
    out_file << "\nTotal time (s): " << total_time  << "\n\n" << endl;
    for (int l=0; l<tasks_out.size(); ++l)
    {
      out_file << "\nProcess " << tasks_out[l]
               << ": " << parameter_strings[tasks_out[l]] << "\n"
               << "Peak time (s): " << peak_times_out[l] << "\n"
               << input_strings_out[l] << "\n";
      for (int a=0; a<batch_size*num_amps; ++a)
      {
        out_file << amplitudes_out[batch_size*num_amps*l+a].real()
                 << " "<< amplitudes_out[batch_size*num_amps*l+a].imag()
                 << " ";
      }
      out_file << "\n";
    }
    out_file << "\n------------ END OUTPUT ------------\n" << flush;
    t1 = high_resolution_clock::now();
    span = duration_cast<duration<double>>(t1 - t0);
    cout << "Writing took " << span.count() << endl;

    // Close files
    out_file.close();

    // Report total time
    t1 = high_resolution_clock::now();
    span = duration_cast<duration<double>>(t1 - t0);
    //cout << "Total time = " << span.count() << endl;
  }

  ///////////////////////////// Process > 0
  else
  {
    MPI_Status status;
    int length(0);
    string local_line;
    vector<s_type> amplitudes;
    double time_largest_contraction;
    vector<s_type> amplitudes_p_time;
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
    local_line = "";
    for (int l=0; l<length; ++l)
      local_line += string(1,char_ptr[l]);
      
    delete[] char_ptr;
    char_ptr = nullptr;

    talsh::initialize(&mem_size);
    cout << mem_size << endl << flush;
    MPI_Barrier(MPI_COMM_WORLD);
    {
      Contraction contraction(local_line, num_args, num_amps);
      //cout << mem_size << endl << flush;

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
        local_line = "";
        for (int l=0; l<length; ++l)
          local_line += string(1,char_ptr[l]);
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
        //cout << "A round took " << span.count() << "s" << endl;

        ///////////////// Message back
        amplitudes = contraction.get_amplitudes();
        time_largest_contraction = contraction.get_time_largest_contraction();
        amplitudes_p_time = vector<s_type>(amplitudes);
        amplitudes_p_time.push_back(s_type(time_largest_contraction));
        MPI_Send(amplitudes_p_time.data(), amplitudes_p_time.size(),
                 MPI_COMPLEX, 0, 0, MPI_COMM_WORLD);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    talsh::shutdown();

  }

  // Finalize the MPI environment.
  MPI_Finalize();

  return 0;
} 
