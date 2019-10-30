# Base OS
FROM dagart/cirq-alpine:latest

# Install baseline
RUN apk update
RUN apk add g++ make gsl-dev git autoconf automake python3-dev py3-pybind11 py3-packaging py3-pytest py3-docopt

ENTRYPOINT ["/bin/sh"]
