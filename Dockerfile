# Base OS
FROM ngnrsaa/cirq-alpine:latest

# Install baseline
RUN apk update
RUN apk add git g++ make gsl-dev git autoconf automake python3-dev py3-pybind11 py3-packaging py3-pytest py3-docopt

# Clone qflex
RUN git clone --depth 1 --shallow-submodules https://github.com/ngnrsaa/qflex.git /qflex/

WORKDIR /qflex/

# Install dependences
RUN autoreconf -i && autoconf && ./configure --disable-all_checks

# Compile qflex
RUN make -j${OMP_NUM_THREADS:-4}

ENTRYPOINT ["/qflex/src/qflex.x"]
