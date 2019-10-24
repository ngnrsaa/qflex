# Use qFlex as baseline
FROM qflex

# Compile tests
WORKDIR /qflex/

ENTRYPOINT make -j${OMP_NUM_THREADS} run-cxx-tests
