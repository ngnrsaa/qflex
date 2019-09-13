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
RUN for file in *.x; do echo /qflex/tests/$file >> run_all.sh; done
RUN chmod 755 run_all.sh

CMD ["/qflex/tests/run_all.sh"]
