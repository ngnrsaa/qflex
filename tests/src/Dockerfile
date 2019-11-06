# Use qFlex as baseline
FROM qflex

# Copy qflex tests
COPY ./tests/src/ /qflex/tests/src/

# Compile tests
WORKDIR /qflex/

ENTRYPOINT make -j${OMP_NUM_THREADS:-4} run-cxx-tests
