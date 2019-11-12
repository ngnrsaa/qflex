# Docker

[Docker](https://docker.com) can be used to run both simulations and tests,
without the necessity to install all the required dependencies to run qFlex.

## Build Docker Images

To build qFlex and run all the tests:

```
# docker-compose up --build
```

`docker-compose` will create three images: `qflex`, `qflex-cxx-tests`, and
`qflex-py-tests`. To build qFlex images without running tests, use the
following command:

```
# docker-compose build
```

## Run Tests

Once `qflex-cxx-tests` is created, use the following command to run all C++ tests:

```
# docker run -ti --rm qflex-cxx-tests
```

Similarly, once `qflex-py-tests` is created the following command will run all
Python tests:

```
# docker run -ti --rm qflex-py-tests
```

## Run Simulations

Once `qflex:latest` is created, use the following command to run a simulation
(see [Input File Formatting](input_formats.md) for more information regarding
the simulation parameters):

```
# docker run -ti --rm -v $PWD/config/circuits:/qflex/config/circuits:ro \
                      -v $PWD/config/ordering:/qflex/config/ordering:ro \
                      -v $PWD/config/grid:/qflex/config/grid:ro \
                      qflex:latest /qflex/config/circuits/bristlecone_48_1-24-1_0.txt \
                                   /qflex/config/ordering/bristlecone_48.txt \
                                   /qflex/config/grid/bristlecone_48.txt
```

The flag `-v [orig]:[dest]:[attr]` is required to allow qFlex image access to
the host folders. However, the original `config/circuits`, `config/ordering` and
`config/grid` folders are copied during the building phase (but not updated if
changes are made after qFlex images are built). To run simulations using the
default input files, `-v` flags may be dropped:

```
# docker run -ti --rm qflex:latest /qflex/config/circuits/bristlecone_48_1-24-1_0.txt \
                                   /qflex/config/ordering/bristlecone_48.txt \
                                   /qflex/config/grid/bristlecone_48.txt
```
