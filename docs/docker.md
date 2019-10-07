# Docker

[Docker](https://docker.com) can be used to run both simulations and tests,
without the necessity to install all the required dependencies to run qFlex.

## Build Docker Images

To build qFlex and all the tests:

```
# docker-compose up
```

`docker-compose` will create two images: `qflex:latest` and `qflex:tests`. To
just build qFlex images, use the following command:

```
# docker-compose build
```

## Run Tests

Once `qflex:tests` is created, use the following command to run all the tests:

```
# docker run -ti --rm qflex:tests
```

## Run Simulations

Once `qflex:latest` is created, use the following command to run a simulation
(see [Input File Formatting](input_format.md) for more information regarding
the simulation parameters):

```
# docker run -ti --rm -v $PWD/circuits:/qflex/circuits:ro \
                      -v $PWD/ordering:/qflex/ordering:ro \
                      -v $PWD/grid:/qflex/grid:ro \
                      qflex:latest 11 12 2 0.005 /qflex/circuits/bristlecone_48_1-40-1_0.txt \
                                                 /qflex/ordering/bristlecone_48.txt \
                                                 /qflex/grid/bristlecone_48.txt
```

The flag `-v [orig]:[dest]:[attr]` is required to allow qFlex image access to
the host folders. However, the original `qflex/circuits`, `qflex/ordering` and
`qflex/grid` folders are copied during the building phase (but not updated if
changes are made after qFlex images are built). To run simulations using the
default input files, `-v` flags may be dropped:

```
# docker run -ti --rm qflex:latest 11 12 2 0.005 /qflex/circuits/ben_11_16_0.txt \
                                                 /qflex/ordering/bristlecone_48.txt \
                                                 /qflex/grid/bristlecone_48.txt
```
