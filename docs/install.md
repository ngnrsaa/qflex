# Install qFlex

qFlex can be built and installed locally using [GNU Autotools](https://www.gnu.org/software/automake/faq/autotools-faq.html).
GNU Autotools will automatically check prerequisites and can install software to
any target folder. By default, qFlex installs to `$HOME/local`.

## Required software

To install qFlex on `Debian/Ubuntu` systems, the following packages must be installed in the system:
```
$ sudo apt-get install autoconf automake make g++ libopenblas-dev
```
followed by the installation of required python modules:
```
$ python3 -m pip install --user -r scripts/requirements.txt
```

## Check prerequisites

Prerequisites can be checked by running:
```
$ autoreconf -i && autoconf && ./configure
```
`autoreconf` and `autoconf` generate a suitable `configure` file for the system and `configure`
will check for prerequisites. If prerequisites are satisfied, a `Makefile` is
created to compile and install qFlex. The default installation folder is
`$HOME/local`. If a different folder is desired, `--prefix` can be passed as an
argument to `configure`:
```
$ autoreconf -i && autoconf && ./configure --prefix=/new/installation/folder/
```
By default, `g++` is used to compile qFlex. To use a different compiler (such as
Intel `icpc`), the user may provide the bash variable `CXX`:
```
$ autoreconf -i && autoconf && CXX=icpc ./configure
```
When `CXX=icpc` is used, qFlex will try to use the MKL library and return an
error if not available. For further info:
```
$ autoreconf -i && autoconf && ./configure --help
```

The following options for `configure` are also available:
```
--enable-openmp               Enable experimental support for OpenMP
--disable-python_checks       Disable python checks while configuring
--disable-pybind11            Disable installation of Python porting of qFlex
--disable-python_tests        Disable python tests
--disable-error_on_warnings   Warnings are not considered as errors
```
*Note: If you are using ICPC to build qFlex, you **must** specify
--disable-pybind11 in this step. Construction of the python wrappers with ICPC
is currently unsupported.*

## Prerequisites on Mac OS X

qFlex can be natively compiled on Mac OS X by installing few dependencies through
`brew` and `pip`. More precisely, the following libraries must be installed
using `brew`:
```
$ brew install coreutils openblas autoconf automake libomp pybind11
```
and through `pip`:
```
$ python3 -m pip install --user -r scripts/requirements.txt
```
followed by:
```
$ autoreconf -i && autoconf && ./configure
```

## Compile qFlex

To compile qFlex, it suffices to run:
```
$ make && make run-tests
```
The first command will compile a qFlex executable (`src/qflex.x`) that can be
directly used without installation. The second command will run all the
available tests. While not required, it is suggested to make sure that qFlex has
been properly compiled. To speed-up the compilation, it is possible to specify
the number of parallel processes as:
```
$ make -j8 && make run-tests -j8
```
where `-j8` means that `8` parallel processes are used.

## Running circuits

To run a sample simulation, use the following command:

```
$ ./src/qflex.x config/circuits/bristlecone_48_1-24-1_0.txt \
config/ordering/bristlecone_48.txt config/grid/bristlecone_48.txt
```

## Install qFlex

Finally, qFlex can be installed as:
```
$ make install
```
By default, qFlex is installed in `$HOME/local`. See
["Check prerequisites"](#check-prerequisites) to change the default folder.

## Clean-up

All the installation files can be cleaned-up by running:
```
$ make clean-all
```
