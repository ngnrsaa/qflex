# Base OS
FROM alpine:latest

# Install dependencies
RUN apk update
RUN apk add g++ make gsl-dev git autoconf

# Copy qflex
COPY ./ /qflex/

WORKDIR /qflex/

# Install dependences
RUN autoconf && ./configure

# Compile qflex
RUN make -j${OMP_NUM_THREADS}

ENTRYPOINT ["/qflex/src/qflex.x"]
