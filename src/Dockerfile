# Base OS
FROM qflex-base

# Copy qflex
COPY ./ /qflex/

WORKDIR /qflex/

# Install dependences
RUN autoreconf -i && autoconf && ./configure --disable-all_checks

# Compile qflex
RUN make -j${OMP_NUM_THREADS}

ENTRYPOINT ["/qflex/src/qflex.x"]
