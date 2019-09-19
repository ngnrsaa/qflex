# Use qFlex as baseline
FROM qflex

# Compile tests
WORKDIR /qflex/tests/
RUN make

# Create tests script
RUN echo "#!/bin/bash" >> run_all.sh
RUN for file in *.x; do echo /qflex/tests/$file >> run_all.sh; done
RUN chmod 755 run_all.sh

ENTRYPOINT ["/qflex/tests/run_all.sh"]
CMD []
