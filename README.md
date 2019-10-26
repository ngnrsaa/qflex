# qFlex: A flexible high-performance simulator for the verification and benchmarking of quantum circuits implemented on real hardware

## Description

Flexible Quantum Circuit Simulator (qFlex) implements an efficient tensor
network, CPU-based simulator of large quantum circuits. qFlex computes exact
probability amplitudes, a task that proves essential for the verification of
quantum hardware, as well as mimics quantum machines by computing amplitudes
with low fidelity. qFlex targets quantum circuits in the range of sizes expected
for supremacy experiments based on random quantum circuits, in order to verify
and benchmark such experiments.

## Documentation

qFlex scientific documentation and results can be found in
[[1]](https://arxiv.org/abs/1811.09599) and
[[2]](https://arxiv.org/abs/1905.00444). For technical documentation on how to
use qFlex, see [qflex/docs](/docs).

[[1]](https://arxiv.org/abs/1811.09599) B. Villalonga, *et al.*, **"A flexible
high-performance simulator for verifying and benchmarking quantum circuits
implemented on real hardware"**, NPJ Quantum Information 5, 86 (2019)

[[2]](https://arxiv.org/abs/1905.00444) B. Villalonga, *et al.*, **"Establishing
the Quantum Supremacy Frontier with a 281 Pflop/s Simulation"**,
arXiv:1905.00444 (2019)

## Build methods

To ensure cross-platform viability, qFlex supports multiple different build
methods. If one of the build methods below does not work on your system, try
using one of the other methods listed.

Known incompatibilities:

-   MacOS only supports building via Docker.

### Local installation

To build qFlex on your machine, simply run:

```
$ autoreconf -i && autoconf && ./configure
$ make && make run-tests
```

To disable qFlex python interface, use `./configure --disable-pybind11`.
By default, qFlex is installed in `$HOME/local/`. To change the
installation folder, use `./configure --prefix=/new/installation/folder/`. 

After running these commands, qFlex can be installed by running `make install`.
By default, this installs qFlex in `$HOME/local/`. To change the installation
folder, use `./configure --prefix=/new/installation/folder/`.

qFlex provides an experimental support for `OpenMP`. To activate `OpenMP`, run
`./configure` with the extra-option `--enable-openmp`.

For more information, see [the installation guide](/docs/install.md).

### Build Using Docker

[Docker](https://docker.com) allows you to run qFlex in an isolated environment,
without having to worry about managing dependencies.

To build qFlex with Docker, see [the Docker guide](/docs/docker.md).

### Build Using Rootless Containers

Rootless containers are an alternative to Docker targeted at users who may not
have `root` privileges on their machine.

To run qFlex in a rootless container, see
[qflex/docs/rootless-container.md](/docs/rootless-container.md).

## Using Google Cirq

[Cirq](https://github.com/quantumlib/cirq) is a framework for modeling and
invoking Noisy Intermediate Scale Quantum (NISQ) circuits.

To run qFlex on Google Cirq circuits, or just to call the simulator from Python,
see [qflex/docs/cirq_interface.md](/docs/cirq_interface.md).

## License

Copyright © 2019, United States Government as represented by the Administrator
of the National Aeronautics and Space Administration. All rights reserved.

The Flexible Quantum Circuit Simulator (qFlex) framework is licensed under the
Apache License, Version 2.0 (the "License"); you may not use this application
except in compliance with the License. You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0.

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
