# Use qFlex as baseline
FROM qflex

# Arguments from docker-compose
ARG OMP_NUM_THREADS
ARG QFLEX_BRANCH

# Compile tests
WORKDIR /qflex/

ENTRYPOINT make -j${OMP_NUM_THREADS:-8} run-py-tests
