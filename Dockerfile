# Base OS
FROM debian:stable-slim

# Install dependencies
RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get -y install make g++ libgsl-dev libgslcblas0 git

# Copy qflex
COPY ./ /qflex/

WORKDIR /qflex/

# Install dependences
RUN git submodule update --init --recursive

# Compile qflex
RUN make

# Compile tests
WORKDIR /qflex/tests/
RUN make

# Create tests script
RUN echo "#!/bin/bash" >> run_all.sh
RUN echo /qflex/tests/contraction_utils_test.x >> run_all.sh
RUN echo /qflex/tests/evaluate_circuit_test.x >> run_all.sh
RUN echo /qflex/tests/read_circuit_test.x >> run_all.sh
RUN echo /qflex/tests/tensor_test.x >> run_all.sh
RUN chmod 755 run_all.sh

CMD ["/qflex/tests/run_all.sh"]
