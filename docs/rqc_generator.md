# Random Quantum Circuit (RQC) Generator

A Python-based RQC generator is provided in [generator.py](/qflexcirq/circuits/generator.py).

The purpose of this document is to provide information on how to use the provided
RQC generator.

## Usage

### Input

This tool can be called in the following manner:
```
Example usage:
    $ python generator.py --pattern-file=config/patterns/ibm_rochester.txt \
                          --single_qubit_gates=x_1_2,y_1_2,hz_1_2 \
                          --two_qubit_gate=cx \
                          --sequence=ABAC \
                          --depth=20 \
                          --seed=0 \
                          --output=rochester_abac_m20.qsim
```

### Tags
* --pattern_file
    * Name of pattern file to be used. See examples in [/patterns](/patterns)
* --single_qubit_gates
    * Comma-separated list of single-qubit gates known to qsim
* --two_qubit_gate
    * Name of two-qubit gates known to qsim
* --sequence
    * Sequence of device-specific coupler activation patterns
* --depth
    * Number of layers of two-qubit gates
* --seed
    * Random seed to initialize PRNG
* --output
    * Name of output file to write the generated RQC in qsim format
    
### Output

Currently, the generated RQC is not directly compatible with qFlex input. The qubit 
indices generated need to be mapped to the indexing used by qFlex as described in the 
[input_formats](/docs/input_formats.md) file.

## Patterns

To generate a random quantum circuit, a file containing the activation pattern
must be provided (see [/patterns](/patterns) for examples). Activation patterns
must be given in a dictionary-like format (see for instance
[test.txt](/patterns/test.txt)):
```
{
  'A': {(0, 2), (1, 3)},
  'B': {(0, 1), (2, 3)}
}
```
with `A` and `B` being tags for the patterns. Each pattern must be a list of
pairs of qubits where two-qubit gates will be applied. The RQC is created by
applying a layer of single-qubit gates, followed by two-qubit gates accordingly
to the sequence given in `--sequence`.
