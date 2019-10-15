# qFlex: A flexible high-performance simulator for the verification and benchmarking of quantum circuits implemented on real hardware

## Description

Flexible Quantum Circuit Simulator (qFlex) implements an efficient tensor
network, CPU-based simulator of large quantum circuits. qFlex computes exact
probability amplitudes, a task that proves essential for the verification of
quantum hardware, as well as mimics quantum machines by computing amplitudes
with low fidelity.  qFlex targets quantum circuits in the range of sizes
expected for supremacy experiments based on random quantum circuits, in order to
verify and benchmark such experiments.

## Documentation

qFlex scientific documentation and results can be found in 
[[1]](https://arxiv.org/abs/1811.09599) and [[2]](https://arxiv.org/abs/1905.00444).
For technical documentation on how to use qFlex, see [qflex/docs](/docs).

[[1]](https://arxiv.org/abs/1811.09599) B. Villalonga, *et al.*, **"A flexible 
high-performance simulator for verifying and benchmarking quantum circuits 
implemented on real hardware"**, NPJ Quantum Information 5, 86 (2019)

[[2]](https://arxiv.org/abs/1905.00444) B. Villalonga, *et al.*, **"Establishing
the Quantum Supremacy Frontier with a 281 Pflop/s Simulation"**, arXiv:1905.00444 (2019)

## Build Using Docker

To run qFlex in a docker container, see [qflex/docs/docker.md](/docs/docker.md).

## Build Using Rootless Containers

To run qFlex in a rootless container, see
[qflex/docs/rootless-container.md](/docs/rootless-container.md).

## Licence

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
