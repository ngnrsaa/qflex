# Base OS
FROM alpine:latest

# Install dependencies
RUN apk update
RUN apk add g++ make gsl-dev git

# Copy qflex
COPY ./ /qflex/

WORKDIR /qflex/

# Install dependences
RUN git submodule update --init --recursive

# Compile qflex
RUN make -j${OMP_NUM_THREADS}

ENTRYPOINT ["/qflex/qflex.x"]
