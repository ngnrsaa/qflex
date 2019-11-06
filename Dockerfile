# Base OS
FROM ngnrsaa/cirq-alpine:latest

# Install baseline
RUN apk update
RUN apk add git openssh g++ make gsl-dev git autoconf automake python3-dev py3-pybind11 py3-packaging py3-pytest py3-docopt

# Copy qflex
COPY ./ /qflex/

WORKDIR /qflex/

# Install dependences
RUN autoreconf -i && autoconf && ./configure --disable-all_checks

# Compile qflex
RUN make -j${OMP_NUM_THREADS:-4}

ENTRYPOINT ["/qflex/src/qflex.x"]
