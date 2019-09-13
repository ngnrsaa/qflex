# Base OS
FROM debian:stable-slim

# Working Directory
WORKDIR /qflex/

# Install dependencies
RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get -y install make g++ libgsl-dev libgslcblas0

# Copy to new location
COPY ./ /qflex/
RUN cd /qflex/

# Run Makefile
RUN make

# Start with bash
CMD ["bash"]
