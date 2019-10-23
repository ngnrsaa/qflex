# Input File Formatting

For each circuit simulated on qFlex, three input files are required:
1. A circuit file detailing the operations to perform at each timestep.
1. An ordering file indicating the order in which tensors should be combined.
1. A grid file representing the positions of active qubits in a 2D lattice.

The purpose of this document is to outline the proper formatting of these files.

## Circuit files

Circuit files express a series of quantum gates for qFlex to simulate. In these
files, each gate is represented by a __cycle__ (the timestep in which the gate
is performed), an __opcode__ (the type of gate to perform), and a list of
__indices__ denoting which qubits the gate affects.

Sample circuit files can be found under [qflex/config/circuits](/config/circuits).

### File format

__Requirement:__ First line of the file MUST be an integer containing the number of active
qubits in the circuit.  

Each line after that contains a gate in the following format:
```
<cycle> <opcode> <indices>

cycle:      The integer value of the timestep in which the gate is performed, starting at 0.
opcode:     The type of gate to perform from the following list:
            Gates without arguments:
                - 'h': Hadamard gate
                - 'cz': Controlled-Z gate
                - 't': T gate (equivalent to Z^1/4)
                - 'x_1_2': X^1/2 gate
                - 'y_1_2': Y^1/2 gate
            Gates with arguments:
                - 'rz(theta)': Z-rotation by theta radians
                - 'fsim(theta,phi)': fermionic simulation gate
                    - theta: |01> to |10> swap angle
                    - phi: conditional phase angle
                    - Examples: fsim(pi/2,0) = iSWAP; fsim(0,pi) = CZ
indices:     The integer indices of the qubit or qubits that the gate affect. 
             Currently, the indices are zero indexed and run from left to right, top to bottom.
```

Sample circuit files can be found under [qflex/circuits](/circuits).

## Ordering files

Circuit-ordering files allow fine-tuned optimization of how qFlex simulates a
given circuit. As defined in
[the original paper](https://arxiv.org/abs/1905.00444), every simulation begins
with contraction of all qubit worldlines to a 2D grid; the steps taken after
are defined in this file. Each simulation step is either a __patch-expansion__, a __patch-merge__, or a __cut__.
In this notation, a "patch" refers to a section of the tensor grid that has been contracted into a single tensor. 
Multiple patches are used to improve parallelism and reduce the maximum required tensor rank.

Sample ordering files can be found under [qflex/config/ordering](/config/ordering).

### File format

__patch-expansion:__ contracting an index onto a given patch.

Performing this operation expands the section of the grid represented by the
target patch - thus the name.
```
expand <patch_name> <index>

patch_name: A text string representing a tensor-contraction patch.
index:      The integer index of the qubit being contracted onto the patch.
```
Example:
```
expand A 24 -> This contracts the tensor with the index 24 onto the A patch.
```

__patch-merge:__ contracting a source patch into a target patch. The resulting
                 patch takes the name of the target patch.

Merging a patch into a previously-unused patch can be used to copy or rename
patches; note however that the memory cost of a simulation scales with the
number of distinct patches in the ordering file.
```
merge <patch_name> <patch_name>

patch_name: A text string representing a tensor-contraction patch.
```
Example:
```
merge A B -> This contracts the patches A and B; the resulting patch is named B.
```

__cut:__ projection of all indices shared by two tensors onto a single value
         OR projection of final state of a single tensor onto an index

__Restriction:__ Any patches that are modified prior to a cut must not be
modified after that cut. (Note that a patch can act as the source of a
patch-merge without being modified.) To work around this restriction, merge the
old patch onto a previously-unused patch after the cut. See the
`ExampleOrdering` test in [contraction_utils_test](/tests/src/contraction_utils_test.cpp)
for an example of this.

*This restriction is due to a memory-saving feature of qFlex: tensors are not
copied during a cut operation, but are reused for each projected value of the
cut index. Merging with an empty patch explicitly indicates to qFlex that the
tensor should be copied.*
```
cut <values> <indices>

values:     A comma-separated list in parentheses with no spaces of integer values
            to be assigned to the index during the cut. Can be empty.
indices:    The index or indices to apply the cut.
```
Example:
```
cut (0,1) 24 42 -> This performs a cut betwen tensors 24 and 42 and projects the values 1 and 0 on the cut.

cut (0) 64 -> This performs a projection of the final state '0' onto the tensor at index 64.
```

Sample ordering files can be found under [qflex/config/ordering](/config/ordering).



## Grid files

Grid files are a simple list of 1s and 0s indicating the positions of active
qubits in a 2D lattice. These files are whitespace-agnostic, so it is
recommended to arrange them in a format matching the chip they represent.

Sample grid files can be found under [qflex/config/grid](/config/grid).

### File format

Each file should contain a grid of 1s and 0s, where the 1s signify a qubit
that is active, and the 0s signify inactive qubits.

Example:
```
0 0 1 1 0 0
0 1 1 1 1 0
0 1 1 1 1 0
0 0 1 1 0 0
```

