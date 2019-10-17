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
# docker run -ti --rm -v $PWD/config/circuits:/qflex/config/circuits:ro \
                      -v $PWD/config/ordering:/qflex/config/ordering:ro \
                      -v $PWD/config/grid:/qflex/config/grid:ro \
                      qflex:latest 2 /qflex/config/circuits/bristlecone_48_1-24-1_0.txt \
                                     /qflex/config/ordering/bristlecone_48.txt \
                                     /qflex/config/grid/bristlecone_48.txt
```

The flag `-v [orig]:[dest]:[attr]` is required to allow qFlex image access to
the host folders. However, the original `qflex/circuits`, `qflex/ordering` and
`qflex/grid` folders are copied during the building phase (but not updated if
changes are made after qFlex images are built). To run simulations using the
default input files, `-v` flags may be dropped:

```
# docker run -ti --rm qflex:latest 2 /qflex/config/circuits/bristlecone_48_1-24-1_0.txt \
                                     /qflex/config/ordering/bristlecone_48.txt \
                                     /qflex/config/grid/bristlecone_48.txt
```
