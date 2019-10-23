# Use qFlex as baseline
FROM qflex

# Install python components
# TODO: this is pretty slow (~2 minutes). Find a way to install these faster.
RUN apk add libpng-dev freetype-dev
RUN apk add python3-dev python3-tkinter
RUN apk add py3-scipy py-numpy-dev
RUN apk add openblas-dev
RUN pip3 install --upgrade pip
RUN pip3 install -U cirq
RUN pip3 install -U pytest

# Run the test binary
WORKDIR /qflex/
ENTRYPOINT ["pytest"]
CMD ["python/order_circuit_simulation_test.py"]
