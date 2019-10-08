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

### File format

First line of the file should contain the number of qubits in the circuit.
Each line after that contains a gate in the following format:
```
<cycle> <opcode> <indices>

cycle:      The integer value of the timestep in which the gate is performed.
opcode:       The type of gate to perform from the following list:
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
qubits:     The integer indices of the qubit or qubits that the gate affect. 
            Currently, the indices are zero indexed and run from left to right, top to bottom.
```

Sample circuit files can be found under [qflex/circuits](/circuits).

### Formal grammar

```
file:      gates
gates:     gate {newline} gates | gate
gate:      cycle operator qubits | {comment}
cycle:     {index}
operator:  {opcode_no_args} | {opcode_with_args}
qubits:    qubit qubits | qubit
qubit:     {index}
```

### Terminal symbol definitions

```
{newline}: a single line-break. Avoid trailing spaces at the end of lines.
{comment}: a one-line text string whose first character must be "#".
{index}: any integer value, i.e. \[0-9\]+
{opcode_no_args}: the name of an operator which does not take arguments.
These are defined as follows:

- 'h': Hadamard gate
- 'cz': Controlled-Z gate
- 't': T gate (equivalent to Z^1/4)
- 'x_1_2': X^1/2 gate
- 'y_1_2': Y^1/2 gate

{opcode_with_args} the name of an operator that takes arguments, followed by a comma-separated list (no spaces) of arguments in parentheses.
These are defined as follows:

- 'rz(theta)': Z-rotation by theta radians
- 'fsim(theta,phi)': fermionic simulation gate
  - theta: |01> to |10> swap angle
  - phi: conditional phase angle
  - Examples: fsim(pi/2,0) = iSWAP; fsim(0,pi) = CZ
```

## Ordering files

Circuit-ordering files allow fine-tuned optimization of how qFlex simulates a
given circuit. As defined in
[the original paper](https://arxiv.org/abs/1905.00444), every simulation begins
with contraction of all qubit worldlines to a 2D grid; the steps taken after
are defined in this file. Each simulation step is either a __patch-expansion__, a __patch-merge__, or a __cut__. 

### File format

__patch-expansion:__ contracting an index onto a given patch.
```
expand <patch_name> <index>

patch_name: A text string representing a tensor-contraction patch.
index:      The integer index of the qubit being contracted onto the patch.
```
__patch-merge:__ contracting two patches.
```
merge <patch_name> <patch_name>

patch_name: A text string representing a tensor-contraction patch.
```
__cut:__ cutting a tensor or between two tensors
```
cut <values> <indices>

values:     A comma-separated list in parentheses with no spaces of integer values
            to be assigned to the index during the cut. Can be empty.
indices:    The index or indices to apply the cut.
```

Sample ordering files can be found under [qflex/ordering](/ordering).

### Formal grammar

```
file:          steps
steps:         step {newline} steps | step
step:          patch_expand | patch_merge | cut | {comment}
patch_expand:  {expand} {patch_name} {index}
patch_merge:   {merge} {patch_name} {patch_name}
cut:           {cut} {values} {index} | {cut} {values} {index} {index}
```

### Terminal symbol definitions

```
{newline}: a single line-break. Avoid trailing spaces at the end of lines.
{comment}: a one-line text string whose first character must be "#".
{expand}: the string "expand". 
{merge}: the string "merge".
{cut}: the string "cut".
{patch_name}: any text string representing a tensor-contraction patch; e.g. "pB".
{index}: any integer value, i.e. \[0-9\]+
{values}: a comma-separated list (no spaces) of integer values in parentheses. Can be empty.
```

## Grid files

Grid files are a simple list of 1s and 0s indicating the positions of active
qubits in a 2D lattice. These files are whitespace-agnostic, so it is
recommended to arrange them in a format matching the chip they represent.

Sample grid files can be found under [qflex/grid](/grid).

### File format
```
0 0 1 1 0 0
0 1 1 1 1 0
0 1 1 1 1 0
0 0 1 1 0 0
```

### Formal grammar

```
file:    qubits
qubits:  qubit {whitespace} qubits | qubit
qubit:   {zero} | {one}
```

### Terminal symbol definitions

```
{whitespace}: any number of consecutive whitespace characters (tab, space, or newline).
{zero}: the integer, 0.
{one}: the integer, 1.
```

## Notes

There are a few things to take note of before writing these input files:
* Qubits in grid are zero-indexed going left to right, top to bottom. 
* Coordinate system in examples is zero indexed, with the upper left corner as the origin. 
    * (2, 3): Three down from origin, four right of origin.

## Definitions

Defining the terms used in the formal grammar.
* Expand: This is used in the ordering when we are contracting a tensor. 
    * expand A 42: This would contract the tensor with the index 42 onto the A patch.
* Merge: This is used in the ordering where we contract two patches of contracted tensors.
    * merge A B: This contracts the patches A and B
* Cut: This is used in the ordering where we cut a bond between two tensors.
    * cut (0,1) 24 42: This would cut betwen tensors 24 and 42 and project the values 1 and 0 on the cut. 
