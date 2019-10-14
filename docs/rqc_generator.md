# Random Quantum Circuit Generator

As part of input to qFlex, a random quantum circuit (RQC) file must be generated. 
To help with this RQC generation, a python RQC generator has been provided and can be found in here: [generator.py](/python/generator.py).

The purpose of this document is to provide information on how to use the provided
RQC generator.

## Usage

### Input

This RQC can be called in the following manner:
```
Example usage:
    $ python generator.py --pattern-file=patterns/ibm_rochester.txt \
                          --single_qubit_gates=x_1_2,y_1_2,hz_1_2 \
                          --two_qubit_gate=cx \
                          --sequence=ABAC \
                          --depth=20 \
                          --seed=0 \
                          --output=rochester_abac_m20.qsim
```

### Tags
* --pattern_file
    * Name of pattern file to be used, must be stored in [/patterns](/patterns)
* --single_qubit_gates
    * Comma-separated list of single-qubit gates known to qsim
* --two_qubit_gate
    * Name of two-qubit gates known to qsim
* --sequence
    * Sequence of device-specific coupler activation patterns
* --depth
    * Number of layers of two-qubit gates
* --seed
    * Random seed to initialize P
* --output
    * Name of output file to write the generated RQC in qsim format

## Patterns

To generate a random quantum circuit, a file containing activation patterns for a device must be provided in [/patterns](/patterns).
Activation patterns must be given in a specific format as shown below, taken from [test.txt](/patterns/test.txt)
```
{
  'A': {(0, 1)}, 
  'B': {(1, 2)} 
}
```
