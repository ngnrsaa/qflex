# Base OS
FROM alpine:latest

# Install baseline
RUN apk update
RUN apk add g++ make gsl-dev git autoconf automake python3-dev py3-pybind11 py3-packaging py3-pytest py3-docopt
RUN /sbin/apk add g++ make libpng-dev freetype-dev python3-dev py-scipy py-numpy-dev py-cffi glib gtk+3.0 gobject-introspection
RUN /usr/bin/pip3 install pgi
RUN /usr/bin/pip3 install cairocffi
RUN /usr/bin/pip3 install cirq

ENTRYPOINT ["/bin/sh"]
