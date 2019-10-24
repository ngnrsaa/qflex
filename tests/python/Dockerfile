# Use qFlex as baseline
FROM qflex

# Install python components
# TODO: this is pretty slow (~2 minutes). Find a way to install these faster.
# Install dependencies
RUN /sbin/apk add g++ make libpng-dev freetype-dev python3-dev py-scipy py-numpy-dev py-cffi glib gtk+3.0 gobject-introspection
RUN /usr/bin/pip3 install pgi
RUN /usr/bin/pip3 install cairocffi

# Install Cirq
RUN /usr/bin/pip3 install cirq

# Compile tests
WORKDIR /qflex/

ENTRYPOINT make -j${OMP_NUM_THREADS} run-py-tests
