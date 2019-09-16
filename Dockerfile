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

ENTRYPOINT ["/qflex/qflex.x"]
CMD []
