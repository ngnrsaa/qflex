#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <complex>
#include <ctime>
#include <chrono>

#include <iostream>
#include <fstream>
#include <sstream>

#include "mkl_tensor.h"
#include "read_circuit.h"

#include "talshxx.hpp"
#include "talsh_wrapper.h"

#include <omp.h>
#include "pin.c"

using namespace std;
using namespace chrono;

// Input: I J K fidelity filename initial_conf (optional) final_conf (optional)
int main(int argc, char **argv) {

  // Pinning threads to cores. Gets some speedup.
  if (getenv("OMP_NUM_THREADS") == NULL) omp_set_num_threads (1);
  #pragma omp parallel
  { pin_relative (omp_get_thread_num()); }

  // Set precision for the printed floats.
  cout.precision(12);

  // Timing variables.
  high_resolution_clock::time_point t_output_0, t_output_1;
  t_output_0 = high_resolution_clock::now();
  high_resolution_clock::time_point t0, t1;
  duration<double> time_span;

  // Reading input.
  t0 = high_resolution_clock::now();
  if (argc<6) throw logic_error("ERROR: Not enough arguments.");
  const int I = atoi(argv[1]);
  const int J = atoi(argv[2]);
  const int K = atoi(argv[3]);
  double fidelity = atof(argv[4]);
  const int super_dim = (int)pow(DIM,K); // Watch out with the *2 factor!
  const string filename = string(argv[5]);
  string initial_conf(52, '0'), final_conf_B(46, '0');
  vector<string> final_conf_A(1, string(6, '0'));
  if (argc>6)
    initial_conf = string(argv[6]);
  if (argc>7)
    final_conf_B = string(argv[7]);
  if (argc>8)
  {
    final_conf_A = vector<string>(argc-8);
    for (int s=0; s<final_conf_A.size(); ++s)
      final_conf_A[s] = string(argv[s+8]);
  }
  const int num_Cs = final_conf_A.size();
  int num_qubits = I*J;
  t1 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t1 - t0);
  //cout << "Time spent reading input: "
  //     << time_span.count()
  //     << "s\n\n";

  // List of qubits to remove (off) and qubits in A.
  vector<vector<int>> qubits_off({
    {0,0},{0,1},{0,2},{0,3},{0,4}            ,{0,7},{0,8},{0,9},{0,10},{0,11},
    {1,0},{1,1},{1,2},{1,3}                        ,{1,8},{1,9},{1,10},{1,11},
    {2,0},{2,1},{2,2}                                    ,{2,9},{2,10},{2,11},
    {3,0},{3,1}                                                ,{3,10},{3,11},
    {4,0}                                                      ,{4,10},{4,11},
    {5,0}                                                ,{5,9},{5,10},{5,11},
    {6,0}                                          ,{6,8},{6,9},{6,10},{6,11},
    {7,0},{7,1}                              ,{7,7},{7,8},{7,9},{7,10},{7,11},
    {8,0},{8,1},{8,2}                  ,{8,6},{8,7},{8,8},{8,9},{8,10},{8,11},
    {9,0},{9,1},{9,2},{9,3},{9,4},{9,5},{9,6},{9,7},{9,8},{9,9},{9,10},{9,11},
    {10,0},{10,1},{10,2},{10,3},{10,4},{10,5},{10,6},{10,7},{10,8},{10,9},{10,10},{10,11}
    });
  vector<vector<int>> qubits_A({{2,8},{3,8},{4,8},{5,8},
                                      {3,9},{4,9}       });


  // Scratch space to be reused for operations.
  t0 = high_resolution_clock::now();
  s_type * scratch = new s_type[(int)pow(super_dim,7)];
  t1 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t1 - t0);
  //cout << "Time spent reading allocating scratch space: "
  //     << time_span.count()
  //     << "s\n\n";


  // Declaring and then filling 2D grid of tensors.
  vector<vector<MKLTensor>> tensor_grid(I);
  for (int i=0; i<I; ++i)
  {
    tensor_grid[i] = vector<MKLTensor>(J);
  }
  // Scope so that the 3D grid of tensors is destructed.
  {
    // Creating 3D grid of tensors from file.
    t0 = high_resolution_clock::now();
    vector<vector<vector<MKLTensor>>> tensor_grid_3D;
    google_circuit_file_to_grid_of_tensors(filename, I, J, K, initial_conf,
               final_conf_B, qubits_A, qubits_off, tensor_grid_3D, scratch);
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    //cout << "Time spent creating 3D grid of tensors from file: "
    //     << time_span.count()
    //     << "s\n\n";

    // Contract 3D grid onto 2D grid of tensors, as usual.
    t0 = high_resolution_clock::now();
    grid_of_tensors_3D_to_2D(tensor_grid_3D, tensor_grid,
                             qubits_A, qubits_off, scratch);
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    //cout << "Time spent creating 2D grid of tensors from 3D one: "
    //     << time_span.count()
    //     << "s\n\n";
  }
  // Renormalize, in order to work in the range of floats.
  complex<double> norm_factor(1.);
  //cout << "Tensor norms:\n";
  /*
  for (int i=0; i<I; ++i) for (int j=0; j<J; ++j)
  { 
    if (find(qubits_off.begin(),qubits_off.end(),vector<int>({i,j})) !=qubits_off.end()) { continue; }
    cout << "Tensor [" << i << "][" << j << "] = ";
    cout << tensor_grid[i][j].tensor_norm() << endl;
  }
  s_type scalar;
  for (int i=0; i<I; ++i) for (int j=0; j<J; ++j)
  {
    if (find(qubits_off.begin(),qubits_off.end(),vector<int>({i,j}))
        !=qubits_off.end())
    { continue; }
    //scalar = 10./sqrt(tensor_grid[i][j].tensor_norm());
    scalar = 1.0;
    norm_factor *= scalar;
    tensor_grid[i][j].scalar_multiply(scalar);
  }
  cout << endl;
  cout << "Normalized tensor norms:\n";
  for (int i=0; i<I; ++i) for (int j=0; j<J; ++j)
  {
    if (find(qubits_off.begin(),qubits_off.end(),vector<int>({i,j}))
        !=qubits_off.end())
    { continue; }
    cout << "Tensor [" << i << "][" << j << "] = ";
    cout << tensor_grid[i][j].tensor_norm() << endl;
  }
  cout << endl;
  cout << "Global normalization factor = " << norm_factor << endl;
  cout << endl;
  */


  // Leave reorderings for later. Although the first time they multiply
  // they'll be already ordered for the next one :). Left and right might have
  // and impact though. Think about that.

  // Allocating tensors for the projections. Two cuts.
  t0 = high_resolution_clock::now();
  MKLTensor t32({""}, {(int)pow(super_dim,1)});
  MKLTensor t42({""}, {(int)pow(super_dim,3)});
  // And for the c open indices projections.
  MKLTensor c28({""}, {(int)pow(super_dim,2)});
  MKLTensor c38({""}, {(int)pow(super_dim,4)});
  MKLTensor c48({""}, {(int)pow(super_dim,4)});
  MKLTensor c58({""}, {(int)pow(super_dim,2)});
  MKLTensor c39({""}, {(int)pow(super_dim,2)});
  MKLTensor c49({""}, {(int)pow(super_dim,2)});
  t1 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t1 - t0);
  //cout << "Time spent allocating new tensors to hold projections: "
  //     << time_span.count() << "s\n\n";
  // No projections so far.

  /* Hardcode the randomly chosen combinations.
  // Choosing randomly the cuts that will be taken into account for a given
  // fidelity.
  int num_cuts = 1;
  num_cuts *= tensor_grid[0][5].get_index_to_dimension()
                                  .at("(0,5),(0,6)");
  num_cuts *= tensor_grid[1][5].get_index_to_dimension()
                                  .at("(1,5),(1,6)");
  num_cuts *= tensor_grid[2][5].get_index_to_dimension()
                                  .at("(2,5),(2,6)");
  //cout << "Number of cuts is: " << num_cuts << endl;
  int num_cuts_taken = ceil(num_cuts * fidelity);
  double fidelity_taken = (double)num_cuts_taken / (double)num_cuts;
  //cout << "Take randomly " << num_cuts_taken << " out of " << num_cuts
  //     << " cuts for a fidelity of " << fidelity_taken << "." << endl;
  vector<vector<int>> cut_combinations(num_cuts);
  int j=0;
  for (int i0=0; i0<tensor_grid[0][5]
       .get_index_to_dimension().at("(0,5),(0,6)"); ++i0)
  {
    for (int i1=0; i1<tensor_grid[1][5]
         .get_index_to_dimension().at("(1,5),(1,6)"); ++i1)
    {
      for (int i2=0; i2<tensor_grid[2][5]
           .get_index_to_dimension().at("(2,5),(2,6)"); ++i2)
      {
        cut_combinations[j] = vector<int>({i0,i1,i2});
        ++j;
      }
    }
  }
  random_shuffle(cut_combinations.begin(), cut_combinations.end());
  // These are the cuts chosen at random for the fidelity given.
  vector<vector<int>> cut_combinations_taken(cut_combinations.begin(),
                                   cut_combinations.begin()+num_cuts_taken);
  //cout << "The chosen cut combinations are (i0, i1):\n";
  for (auto v : cut_combinations_taken)
  {
    //for (auto w : v)
    //  cout << w << " ";
    //cout << endl;
  }
  */
  vector<vector<int>> cut_combinations_taken({{0}});

  // Allocating tensors to be reused.
  // First, helper tensors (H_...) with particular sizes.
  // Sometimes I want to allocate helper tensors (or others) in a {} scope,
  // so that the destructor is called and memory freed. It all depends on the
  // particular case.
  t0 = high_resolution_clock::now();
  MKLTensor H_2_legs_a({""}, {(int)pow(super_dim,2)});
  MKLTensor H_2_legs_b({""}, {(int)pow(super_dim,2)});
  MKLTensor H_3_legs_a({""}, {(int)pow(super_dim,3)});
  MKLTensor H_3_legs_b({""}, {(int)pow(super_dim,3)});
  MKLTensor H_4_legs_a({""}, {(int)pow(super_dim,4)});
  MKLTensor H_4_legs_b({""}, {(int)pow(super_dim,4)});
  MKLTensor H_4_legs_c({""}, {(int)pow(super_dim,4)});
  MKLTensor H_5_legs_a({""}, {(int)pow(super_dim,5)});
  MKLTensor H_5_legs_b({""}, {(int)pow(super_dim,5)});
  MKLTensor H_6_legs_a({""}, {(int)pow(super_dim,6)});
  MKLTensor H_6_legs_b({""}, {(int)pow(super_dim,6)});
  MKLTensor H_7_legs_a({""}, {(int)pow(super_dim,7)});
  MKLTensor H_7_legs_b({""}, {(int)pow(super_dim,7)});
  MKLTensor H_7_legs_c({""}, {(int)pow(super_dim,7)});
  // Reused tensors that are not helpers.
  MKLTensor C({""}, {(int)pow(super_dim,4)});
  MKLTensor S({""}, {1}); // Final scalar!
  t1 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t1 - t0);
  //cout << "\nTime spent allocating tensors that will be reused: "
  //     << time_span.count() << "s\n\n";
  // These, plus the tensor_grid and the t... tensors where projections will
  // be stored, are ALL the tensors being held at once. This gives the memory
  // footprint (well, this plus the scratch space allocated in the beginning).


  double time_in_loops = 0.0; // Keeps on adding time spent on different steps.
  // Reorder tensors in C.
  t0 = high_resolution_clock::now();
  tensor_grid[3][9].reorder({"(3,9),(o)","(3,9),(4,9)","(3,8),(3,9)"},
                             scratch);
  tensor_grid[4][9].reorder({"(4,9),(o)","(3,9),(4,9)","(4,8),(4,9)"},
                             scratch);
  tensor_grid[2][8].reorder({"(2,8),(o)","(2,8),(3,8)","(2,7),(2,8)"},
                             scratch);
  tensor_grid[3][8].reorder({"(3,8),(o)","(2,8),(3,8)","(3,8),(3,9)",
                             "(3,8),(4,8)","(3,7),(3,8)"}, scratch);
  tensor_grid[4][8].reorder({"(4,8),(o)","(3,8),(4,8)","(4,8),(4,9)",
                             "(4,8),(5,8)","(4,7),(4,8)"}, scratch);
  tensor_grid[5][8].reorder({"(5,8),(o)","(4,8),(5,8)","(5,7),(5,8)"},
                             scratch);
  t1 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t1 - t0);
  //cout << "Time spent reordering tensors in C: "
  //     << time_span.count() << "s\n\n";
  time_in_loops += time_span.count();

  // First and only loop deals with (i0).
  int i0;
  // Push back cut contributions.
  vector<vector<complex<double>>> amplitudes(num_Cs);
  // Reorder first two layers of projections.
  t0 = high_resolution_clock::now();
  tensor_grid[3][2].reorder({"(3,2),(4,2)","(3,2),(3,3)"}, scratch);
  tensor_grid[4][2].reorder({"(3,2),(4,2)","(4,1),(4,2)",
                             "(4,2),(4,3)","(4,2),(5,2)"}, scratch);
  t1 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t1 - t0);
  //cout << "Time spent reordering two layers to project before iterations: "
  //     << time_span.count() << "s\n\n";
  time_in_loops += time_span.count();
  
  for (int cut=0; cut<cut_combinations_taken.size(); ++cut)
  {
    i0 = cut_combinations_taken[cut][0];

    //cout << "Contracting case " << i0 << " " << i1 << " " << i2 << "\n\n";
    t0 = high_resolution_clock::now();
    // Project.
    tensor_grid[3][2].project("(3,2),(4,2)", i0, t32);
    tensor_grid[4][2].project("(3,2),(4,2)", i0, t42);

    // Build A.
    multiply(tensor_grid[4][1], tensor_grid[5][1], H_3_legs_a, scratch);
    multiply(H_3_legs_a, tensor_grid[6][1], H_3_legs_b, scratch);
    multiply(H_3_legs_b, t42,               H_4_legs_a, scratch);
    multiply(H_4_legs_a, tensor_grid[5][2], H_4_legs_b, scratch);
    multiply(H_4_legs_b, tensor_grid[6][2], H_4_legs_a, scratch);
    multiply(H_4_legs_a, tensor_grid[7][2], H_4_legs_b, scratch);
    multiply(H_4_legs_b, tensor_grid[4][3], H_6_legs_a, scratch);
    multiply(H_6_legs_a, tensor_grid[5][3], H_6_legs_b, scratch);
    multiply(H_6_legs_b, tensor_grid[6][3], H_6_legs_a, scratch);
    multiply(H_6_legs_a, tensor_grid[7][3], H_6_legs_b, scratch);
    multiply(H_6_legs_b, tensor_grid[8][3], H_6_legs_a, scratch);
    multiply(H_6_legs_a, tensor_grid[8][4], H_7_legs_a, scratch);
    multiply(H_7_legs_a, tensor_grid[7][4], H_7_legs_b, scratch);
    multiply(H_7_legs_b, tensor_grid[6][4], H_7_legs_a, scratch);
    multiply(H_7_legs_a, tensor_grid[5][4], H_7_legs_b, scratch);
    multiply(H_7_legs_b, tensor_grid[4][4], H_7_legs_a, scratch);
    multiply(H_7_legs_a, tensor_grid[8][5], H_7_legs_b, scratch);
    multiply(H_7_legs_b, tensor_grid[7][5], H_7_legs_a, scratch);
    multiply(H_7_legs_a, tensor_grid[6][5], H_7_legs_b, scratch);
    multiply(H_7_legs_b, tensor_grid[5][5], H_7_legs_a, scratch);
    multiply(H_7_legs_a, tensor_grid[4][5], H_7_legs_b, scratch);
    multiply(H_7_legs_b, tensor_grid[7][6], H_7_legs_a, scratch);
    multiply(H_7_legs_a, tensor_grid[6][6], H_7_legs_b, scratch);
    multiply(H_7_legs_b, tensor_grid[5][6], H_7_legs_a, scratch);
    multiply(H_7_legs_a, tensor_grid[4][6], H_7_legs_b, scratch);
    multiply(H_7_legs_b, tensor_grid[6][7], H_7_legs_a, scratch);
    multiply(H_7_legs_a, tensor_grid[5][7], H_7_legs_b, scratch);
    multiply(H_7_legs_b, tensor_grid[4][7], H_7_legs_c, scratch);
    //cout << "Built A.\n\n";

    // Build B.
    multiply(t32, tensor_grid[3][3], H_3_legs_a, scratch);
    multiply(H_3_legs_a, tensor_grid[2][3], H_3_legs_b, scratch);
    multiply(H_3_legs_b, tensor_grid[3][4], H_5_legs_a, scratch);
    multiply(H_5_legs_a, tensor_grid[2][4], H_5_legs_b, scratch);
    multiply(H_5_legs_b, tensor_grid[1][4], H_5_legs_a, scratch);
    multiply(H_5_legs_a, tensor_grid[3][5], H_7_legs_a, scratch);
    multiply(H_7_legs_a, tensor_grid[2][5], H_7_legs_b, scratch);
    multiply(H_7_legs_b, tensor_grid[1][5], H_7_legs_a, scratch);
    multiply(H_7_legs_a, tensor_grid[0][5], H_7_legs_b, scratch);
    multiply(H_7_legs_b, tensor_grid[0][6], H_7_legs_a, scratch);
    multiply(H_7_legs_a, tensor_grid[1][6], H_7_legs_b, scratch);
    multiply(H_7_legs_b, tensor_grid[2][6], H_7_legs_a, scratch);
    multiply(H_7_legs_a, tensor_grid[3][6], H_7_legs_b, scratch);
    multiply(H_7_legs_b, tensor_grid[1][7], H_7_legs_a, scratch);
    multiply(H_7_legs_a, tensor_grid[2][7], H_7_legs_b, scratch);
    multiply(H_7_legs_b, tensor_grid[3][7], H_7_legs_a, scratch);
    //cout << "Built B.\n\n";

    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    //cout << "Time spent building A and B: "
    //     << time_span.count() << "s\n\n";
    time_in_loops += time_span.count();
    // Final scalar.
    t0 = high_resolution_clock::now();
    // AB
    multiply(H_7_legs_c, H_7_legs_a, H_4_legs_c, scratch);
    //cout << "Built AB.\n\n";
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    //cout << "Time spent contracting AB: "
    //   << time_span.count() << "s\n\n";
    time_in_loops += time_span.count();
    // Finally, contract AB and C into S (scalar).
    //cout << "Number of zeros in A = ";
    //cout << H_7_legs_c.num_zeros() << endl;
    //cout << "Number of zeros in B = ";
    //cout << H_7_legs_b.num_zeros() << endl;
    //cout << "Number of zeros in AB = ";
    //cout << H_6_legs_c.num_zeros() << endl;
    for (int c=0; c<num_Cs; ++c)
    {
      t0 = high_resolution_clock::now();
      string string_A = final_conf_A[c];
      tensor_grid[2][8].project("(2,8),(o)", (string_A[0]=='0')?0:1, c28);
      tensor_grid[3][8].project("(3,8),(o)", (string_A[1]=='0')?0:1, c38);
      tensor_grid[4][8].project("(4,8),(o)", (string_A[2]=='0')?0:1, c48);
      tensor_grid[5][8].project("(5,8),(o)", (string_A[3]=='0')?0:1, c58);
      tensor_grid[3][9].project("(3,9),(o)", (string_A[4]=='0')?0:1, c39);
      tensor_grid[4][9].project("(4,9),(o)", (string_A[5]=='0')?0:1, c49);
      
      multiply(c39, c49, H_2_legs_a, scratch);
      multiply(c28, c38, H_4_legs_a, scratch);
      multiply(H_4_legs_a, H_2_legs_a, H_4_legs_b, scratch);
      multiply(c48, c58, H_4_legs_a, scratch);
      multiply(H_4_legs_b, H_4_legs_a, C, scratch);

      t1 = high_resolution_clock::now();
      time_span = duration_cast<duration<double>>(t1 - t0);
      //cout << "Time spent contracting C: "
      //   << time_span.count() << "s\n\n";
      time_in_loops += time_span.count();

      t0 = high_resolution_clock::now();
      multiply(H_4_legs_c, C, S, scratch);
      // Print resulting scalar.
      //S.print();
      //S.print_data();
      amplitudes[c].push_back(*S.data());
      t1 = high_resolution_clock::now();
      time_span = duration_cast<duration<double>>(t1 - t0);
      //cout << "Time spent contracting AB with C: "
      //   << time_span.count() << "s\n\n";
      time_in_loops += time_span.count();
    }
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    //cout << "Time spent contracting ABCs: "
    //   << time_span.count() << "s\n\n";
    time_in_loops += time_span.count();
  } // End of first loop (i0, i1, i2, i3, chosen cuts).
  //cout << "All contractions took: " << time_in_loops << "s." << endl;

  // Renormalize results.
  for (int c=0; c<amplitudes.size(); ++c)
  {
    for (int i=0; i<amplitudes[c].size(); ++i)
      amplitudes[c][i] /= norm_factor;
  }

  // Add up amplitudes.
  vector<complex<double>> final_result(num_Cs, 0.0);
  for (int c=0; c<amplitudes.size(); ++c)
  {
    for (int i=0; i<amplitudes[c].size(); ++i)
      final_result[c] += amplitudes[c][i];
  }

  // Printing output
  for (int c=0; c<num_Cs; ++c)
  {
    cout << initial_conf << " ";
    cout << final_conf_B << " ";
    cout << final_conf_A[c] << " ";
    int num_paths = cut_combinations_taken.size();
    for (int i=0; i<num_paths; ++i)
    {
      cout << cut_combinations_taken[i][0] << " "
           << real(amplitudes[c][i]) << " "
           << imag(amplitudes[c][i]) << " ";
    }
    cout << real(final_result[c]) << " " << imag(final_result[c]);
    cout << endl;
  }


  // Freeing scratch data: delete and NULL.
  delete[] scratch;
  scratch = NULL;

  // Final time
  t_output_1 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t_output_1 - t_output_0);
  cout << time_span.count() << " s\n\n";

  return 0;
} 
