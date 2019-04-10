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

#include <omp.h>
#include "pin.c"

using namespace std;
using namespace chrono;

// Input: I J K fidelity filename initial_conf (optional) final_conf (optional)
int main(int argc, char **argv) {

  if (getenv("OMP_NUM_THREADS") == NULL) omp_set_num_threads (1);
  #pragma omp parallel
  { pin_relative (omp_get_thread_num()); }

  // Timing variables.
  high_resolution_clock::time_point t0, t1;
  duration<double> time_span;

  // Reading input.
  t0 = high_resolution_clock::now();
  if (argc<6) throw logic_error("ERROR: Not enough arguments.");
  const int I = atoi(argv[1]);
  const int J = atoi(argv[2]);
  const int K = atoi(argv[3]);
  double fidelity = atof(argv[4]);
  if (fidelity>1./32.) throw logic_error("ERROR: fidelity is too large");
  const int super_dim = (int)pow(DIM,K);
  const string filename = string(argv[5]);
  string initial_conf(I*J, '0'), final_conf(I*J, '0');
  if (argc > 6)
  {
    initial_conf = string(argv[6]);
  }
  if (argc > 7)
  {
    final_conf = string(argv[7]);
  }
  t1 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t1 - t0);
  cout << "Time spent reading input: "
       << time_span.count()
       << "s\n\n";

  // Scratch space to be reused for operations.
  t0 = high_resolution_clock::now();
  s_type * scratch = new s_type[(int)pow(super_dim,6)];
  t1 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t1 - t0);
  cout << "Time spent reading allocating scratch space: "
       << time_span.count()
       << "s\n\n";

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
                                       final_conf, tensor_grid_3D, scratch);
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    cout << "Time spent creating 3D grid of tensors from file: "
         << time_span.count()
         << "s\n\n";

    // Contract 3D grid onto 2D grid of tensors, as usual.
    t0 = high_resolution_clock::now();
    grid_of_tensors_3D_to_2D(tensor_grid_3D, tensor_grid, scratch);
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    cout << "Time spent creating 2D grid of tensors from 3D one: "
         << time_span.count()
         << "s\n\n";
  }
  complex<double> norm_factor(1.);
  cout << "Tensor norms:\n";
  for (int i=0; i<I; ++i) for (int j=0; j<J; ++j)
  {
    cout << "Tensor [" << i << "][" << j << "] = ";
    cout << tensor_grid[i][j].tensor_norm() << endl;
  }
  s_type scalar;
  for (int i=0; i<I; ++i) for (int j=0; j<J; ++j)
  {
    scalar = 10./sqrt(tensor_grid[i][j].tensor_norm());
    norm_factor *= scalar;
    tensor_grid[i][j].scalar_multiply(scalar);
  }
  cout << endl;
  cout << "Normalized tensor norms:\n";
  for (int i=0; i<I; ++i) for (int j=0; j<J; ++j)
  {
    cout << "Tensor [" << i << "][" << j << "] = ";
    cout << tensor_grid[i][j].tensor_norm() << endl;
  }
  cout << endl;
  cout << "Global normalization factor = " << norm_factor << endl;
  cout << endl;
  

  // Leave reorderings for later. Although the first time they multiply
  // they'll be already ordered for the next one :). Left and right might have
  // and impact though. Think about that.

  // Allocating tensors for the projections. Four cuts.
  t0 = high_resolution_clock::now();
  MKLTensor t26({""}, {(int)pow(super_dim,2)});
  MKLTensor t36({""}, {(int)pow(super_dim,2)});
  MKLTensor t62({""}, {(int)pow(super_dim,2)});
  MKLTensor t63({""}, {(int)pow(super_dim,2)});
  t1 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t1 - t0);
  cout << "Time spent allocating new tensors to hold projections: "
       << time_span.count() << "s\n\n";
  // No projections so far.

  // Choosing randomly the cuts that will be taken into account for a given
  // fidelity.
  int num_cuts = 1;
  num_cuts *= tensor_grid[2][6].get_index_to_dimension().at("(2,6),(3,6)");
  num_cuts *= tensor_grid[6][2].get_index_to_dimension().at("(6,2),(6,3)");
  // Only consider two right cuts for randomly choosing a
  // fraction of them. The bottom and left one will always be exhausted.
  // It doesn't have to be like that, but there is a better reuse this way.
  cout << "Number of cuts is: " << num_cuts << endl;
  int num_cuts_taken = ceil(num_cuts * fidelity);
  double fidelity_taken = (double)num_cuts_taken / (double)num_cuts;
  cout << "Take randomly " << num_cuts_taken << " out of " << num_cuts
       << " cuts for a fidelity of " << fidelity_taken << "." << endl;
  int cut_r = rand() % tensor_grid[2][6]
                       .get_index_to_dimension().at("(2,6),(3,6)");
  vector<int> cuts_b(tensor_grid[6][2].get_index_to_dimension()
                     .at("(6,2),(6,3)"));
  int j=0;
  for (int ib=0; ib<tensor_grid[6][2]
       .get_index_to_dimension().at("(6,2),(6,3)"); ++ib)
  {
    cuts_b[j] = ib;
    ++j;
  }
  random_shuffle(cuts_b.begin(), cuts_b.end());
  // These are the cuts chosen at random for the fidelity given.
  vector<int> cuts_b_taken(cuts_b.begin(), cuts_b.begin()+num_cuts_taken);
  cout << "The chosen cut combinations are (ir, ib):\n";
  for (auto v : cuts_b_taken)
  {
    cout << cut_r << " " << v;
    cout << endl;
  }
  cout << endl;


  // Allocating tensors to be reused.
  // First, helper tensors (H_...) with particular sizes.
  // Sometimes I want to allocate helper tensors (or others) in a {} scope,
  // so that the destructor is called and memory freed. It all depends on the
  // particular case.
  t0 = high_resolution_clock::now();
  MKLTensor H_2_legs_a({""}, {(int)pow(super_dim,2)});
  MKLTensor H_3_legs_a({""}, {(int)pow(super_dim,3)});
  MKLTensor H_3_legs_b({""}, {(int)pow(super_dim,3)});
  MKLTensor H_4_legs_a({""}, {(int)pow(super_dim,4)});
  MKLTensor H_4_legs_b({""}, {(int)pow(super_dim,4)});
  MKLTensor H_5_legs_a({""}, {(int)pow(super_dim,5)});
  MKLTensor H_5_legs_b({""}, {(int)pow(super_dim,5)});
  MKLTensor H_6_legs_a({""}, {(int)pow(super_dim,6)});
  MKLTensor H_6_legs_b({""}, {(int)pow(super_dim,6)});
  MKLTensor H_6_legs_c({""}, {(int)pow(super_dim,6)});
  // Reused tensors that are not helpers.
  MKLTensor AB({""}, {(int)pow(super_dim,6)});
  MKLTensor pC({""}, {(int)pow(super_dim,6)});
  MKLTensor pD({""}, {(int)pow(super_dim,6)});
  MKLTensor S({""}, {1}); // Final scalar!
  t1 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t1 - t0);
  cout << "Time spent allocating tensors that will be reused: "
       << time_span.count() << "s\n\n";
  // These, plus the tensor_grid and the t... tensors where projections will
  // be stored, are ALL the tensors being held at once. This gives the memory
  // footprint (well, this plus the scratch space allocated in the beginning).


  double time_in_loops = 0.0; // Keeps on adding time spent on the contraction.
  // Before the first loop, AB and pC, pD can be built.
  // A
  t0 = high_resolution_clock::now();
  multiply(tensor_grid[0][0], tensor_grid[0][1], H_3_legs_a, scratch);
  multiply(H_3_legs_a, tensor_grid[1][0], H_4_legs_a, scratch);
  multiply(H_4_legs_a, tensor_grid[1][1], H_4_legs_b, scratch);
  multiply(H_4_legs_b, tensor_grid[0][2], H_5_legs_a, scratch);
  multiply(H_5_legs_a, tensor_grid[1][2], H_5_legs_b, scratch);
  multiply(H_5_legs_b, tensor_grid[2][0], H_6_legs_a, scratch);
  multiply(H_6_legs_a, tensor_grid[2][1], H_6_legs_b, scratch);
  multiply(H_6_legs_b, tensor_grid[2][2], H_6_legs_c, scratch);
  cout << "Built A.\n\n";
  t1 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t1 - t0);
  cout << "Overhead time for the contraction of A, before any iteration: "
       << time_span.count() << "s\n\n";
  time_in_loops += time_span.count();

  // B
  t0 = high_resolution_clock::now();
  // Reorder cut.
  tensor_grid[2][6].reorder({"(2,6),(3,6)","(2,5),(2,6)","(1,6),(2,6)"},
                            scratch);
  tensor_grid[3][6].reorder({"(2,6),(3,6)","(3,5),(3,6)","(3,6),(4,6)"},
                            scratch);
  // Project cut.
  tensor_grid[2][6].project("(2,6),(3,6)", cut_r, t26);
  tensor_grid[3][6].project("(2,6),(3,6)", cut_r, t36);
  // Contract.
  multiply(tensor_grid[0][6], tensor_grid[1][6], H_3_legs_a, scratch);
  multiply(H_3_legs_a, t26, H_3_legs_b, scratch);
  multiply(H_3_legs_b, tensor_grid[0][5], H_4_legs_a, scratch);
  multiply(H_4_legs_a, tensor_grid[1][5], H_4_legs_b, scratch);
  multiply(H_4_legs_b, tensor_grid[2][5], H_4_legs_a, scratch);
  multiply(H_4_legs_a, tensor_grid[0][4], H_5_legs_a, scratch);
  multiply(H_5_legs_a, tensor_grid[1][4], H_5_legs_b, scratch);
  multiply(H_5_legs_b, tensor_grid[2][4], H_5_legs_a, scratch);
  multiply(H_5_legs_a, tensor_grid[0][3], H_6_legs_a, scratch);
  multiply(H_6_legs_a, tensor_grid[1][3], H_6_legs_b, scratch);
  multiply(H_6_legs_b, tensor_grid[2][3], H_6_legs_a, scratch);
  cout << "Built B.\n\n";
  t1 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t1 - t0);
  cout << "Overhead time for the contraction of B, before any iteration: "
       << time_span.count() << "s\n\n";
  time_in_loops += time_span.count();
  // AB.
  t0 = high_resolution_clock::now();
  multiply(H_6_legs_a, H_6_legs_c, AB, scratch);
  cout << "Contracted AB.\n\n";
  t1 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t1 - t0);
  cout << "Overhead time for the contraction of AB, before any iteration: "
       << time_span.count() << "s\n\n";
  time_in_loops += time_span.count();

  // pC.
  t0 = high_resolution_clock::now();
  multiply(tensor_grid[6][0], tensor_grid[6][1], H_3_legs_a, scratch);
  multiply(H_3_legs_a, tensor_grid[5][0], H_4_legs_a, scratch);
  multiply(H_4_legs_a, tensor_grid[5][1], H_4_legs_b, scratch);
  multiply(H_4_legs_b, tensor_grid[4][0], H_5_legs_a, scratch);
  multiply(H_5_legs_a, tensor_grid[4][1], H_5_legs_b, scratch);
  multiply(H_5_legs_b, tensor_grid[3][0], H_6_legs_a, scratch);
  multiply(H_6_legs_a, tensor_grid[3][1], pC, scratch);
  t1 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t1 - t0);
  cout << "Overhead time for the contraction of pC, before any iteration: "
       << time_span.count() << "s\n\n";
  time_in_loops += time_span.count();

  // pD.
  t0 = high_resolution_clock::now();
  multiply(tensor_grid[6][6], tensor_grid[5][6], H_3_legs_a, scratch);
  multiply(H_3_legs_a, tensor_grid[4][6], H_4_legs_a, scratch);
  multiply(H_4_legs_a, t36, H_4_legs_b, scratch);
  multiply(H_4_legs_b, tensor_grid[6][5], H_5_legs_a, scratch);
  multiply(H_5_legs_a, tensor_grid[5][5], H_5_legs_b, scratch);
  multiply(H_5_legs_b, tensor_grid[4][5], H_5_legs_a, scratch);
  multiply(H_5_legs_a, tensor_grid[3][5], H_5_legs_b, scratch);
  multiply(H_5_legs_b, tensor_grid[6][4], H_6_legs_a, scratch);
  multiply(H_6_legs_a, tensor_grid[5][4], H_6_legs_b, scratch);
  multiply(H_6_legs_b, tensor_grid[4][4], H_6_legs_a, scratch);
  multiply(H_6_legs_a, tensor_grid[3][4], pD, scratch);
  t1 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t1 - t0);
  cout << "Overhead time for the contraction of pD, before any iteration: "
       << time_span.count() << "s\n\n";
  time_in_loops += time_span.count();


  // The only loop deals with the bottom cut. ib iterates over its
  // values, and t62 and t63 are the projections of the tensors involved.
  vector<s_type> amplitudes; // Keeps on pushing back cut contributions
  for (int ib=0; ib<tensor_grid[6][2]
       .get_index_to_dimension().at("(6,2),(6,3)"); ++ib)
  {
    // If this cut was not randomly chosen, then continue.
    if (find(cuts_b_taken.begin(),cuts_b_taken.end(),ib)==cuts_b_taken.end())
    {
      continue;
    }
    cout << "Contracting case " << cut_r << " " << ib << "\n\n";
    t0 = high_resolution_clock::now();
    // Prepare (reorder) and project.
    tensor_grid[6][2].reorder({"(6,2),(6,3)","(6,1),(6,2)","(5,2),(6,2)"},
                              scratch);
    tensor_grid[6][3].reorder({"(6,2),(6,3)","(6,3),(6,4)","(5,3),(6,3)"},
                              scratch);
    // Project.
    tensor_grid[6][2].project("(6,2),(6,3)", ib, t62);
    tensor_grid[6][3].project("(6,2),(6,3)", ib, t63);

    // Build 
    // Build C.
    multiply(t62, tensor_grid[5][2], H_4_legs_a, scratch);
    multiply(pC, H_4_legs_a, H_6_legs_a, scratch);
    multiply(H_6_legs_a, tensor_grid[4][2], H_6_legs_b, scratch);
    multiply(H_6_legs_b, tensor_grid[3][2], H_6_legs_c, scratch);
    // Build D.
    multiply(t63, tensor_grid[5][3], H_4_legs_a, scratch);
    multiply(pD, H_4_legs_a, H_6_legs_a, scratch);
    multiply(H_6_legs_a, tensor_grid[4][3], H_6_legs_b, scratch);
    multiply(H_6_legs_b, tensor_grid[3][3], H_6_legs_a, scratch);
    // Contract CD.
    multiply(H_6_legs_c, H_6_legs_a, H_6_legs_b, scratch);
    // Get final scalar.
    cout << "Number of zeros in AB = ";
    cout << AB.num_zeros() << endl;
    cout << "Number of zeros in CD = ";
    cout << H_6_legs_b.num_zeros() << endl;
    multiply(AB, H_6_legs_b, S, scratch);
    // Print resulting scalar.
    S.print();
    S.print_data();
    amplitudes.push_back(*S.data());
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    cout << "Time spent on each ib completion: "
       << time_span.count() << "s\n\n";
    time_in_loops += time_span.count();
  } // End of only loop (ib).
  cout << "All contractions took: " << time_in_loops << "s." << endl;
  s_type final_result = 0.0;
  for (int i=0; i<amplitudes.size(); ++i)
    final_result += amplitudes[i];
  for (auto v : amplitudes)
    cout << v << " ";
  cout << endl;
  cout << "Final result (amplitude) is: " << final_result << endl;

  // Freeing scratch data: delete and NULL.
  delete[] scratch;
  scratch = NULL;

  return 0;
}
