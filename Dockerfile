# Base OS
FROM ngnrsaa/cirq-alpine:latest

# Install baseline
RUN apk update
RUN apk add git g++ make gsl-dev git autoconf automake python3-dev py3-pybind11 py3-packaging py3-pytest py3-docopt

# Arguments from docker-compose
ARG OMP_NUM_THREADS

# Copy qflex
COPY ./ /qflex

# Move to the right folder
WORKDIR /qflex/

# Install dependences
RUN autoreconf -i && autoconf && ./configure --disable-all_checks --disable-error_on_warnings

# Make sure that there are no compiled libraries from the host
RUN make clean

# Compile qflex
RUN make -j${OMP_NUM_THREADS:-8}

ENTRYPOINT ["/qflex/src/qflex.x"]
