# Use qFlex as baseline
FROM qflex

# Compile tests
WORKDIR /qflex/tests/
RUN make -j${OMP_NUM_THREADS}

# Create tests script
RUN echo "#!/bin/sh" > /qflex/tests/run_all.sh; \
    COUNT=$(ls -1 /qflex/tests/*.x | wc -l); \
    IDX=1; \
    for file in /qflex/tests/*.x; do \
      if [[ $IDX -lt $COUNT ]]; then \
        if [[ $IDX -gt 1 ]]; then \
          echo -ne "\t"; \
        fi; \
        echo "$file && \ "; \
      else \
        echo -e "\t$file"; \
      fi; \
      IDX=$(echo $IDX+1 | bc -l); \
    done >> /qflex/tests/run_all.sh; \
    chmod 755 /qflex/tests/run_all.sh

ENTRYPOINT ["/qflex/tests/run_all.sh"]
CMD []
