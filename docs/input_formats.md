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

Sample circuit files can be found under [qflex/share/circuits](/share/circuits) with the
following format: `[circuit_type]_[circuit_size]_[circuit_depth]_[idx].txt`.

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
that are defined in this file. Each simulation step has an __operation__
(either a patch-expansion, a patch-merge, or a cut) and some combination of
__indices__ or __patch names__ to which the operation applies.

Sample ordering files can be found under [qflex/share/ordering](/share/ordering).

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

Sample grid files can be found under [qflex/share/grid](/share/grid).

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
