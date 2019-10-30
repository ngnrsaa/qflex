# Use qFlex as baseline
FROM qflex

# Copy qflex tests
COPY ./python/ /qflex/python/
COPY ./tests/python/ /qflex/tests/python/

# Compile tests
WORKDIR /qflex/

ENTRYPOINT make -j${OMP_NUM_THREADS} run-py-tests
