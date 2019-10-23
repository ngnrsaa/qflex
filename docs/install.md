# Install qFlex

qFlex can be built and installed locally using [GNU Autotools](https://www.gnu.org/software/automake/faq/autotools-faq.html).
GNU Autotools will automatically check prerequisites and can install software to
any target folder. By default, qFlex installs to `$HOME/local`.

## Required software

To install qFlex, `autoconf` must be installed in the system:
```
$ sudo apt-get install autoconf
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
--disable-pybind11    Disable installation of Python porting of qFlex
--disable-cirq_tests  Disable tests which depends on cirq
```

## Compile qFlex

To compile qFlex, it suffices to run:
```
make && make run-tests
```
The first command will compile a qFlex executable (`src/qflex.x`) that can be
directly used without installation. The second command will run all the
available tests. While not required, it is suggested to make sure that qFlex has
been properly compiled. To speed-up the compilation, it is possible to specify
the number of parallel processes as:
```
make -j8 && make run-tests -j8
```
where `-j8` means that `8` parallel processes are used.

## Install qFlex

Finally, qFlex can be installed as:
```
make install
```
By default, qFlex is installed in `$HOME/local`. See
["Check PreRequisites"](#check-prerequisites) to change the default folder.

## Clean-up

All the installation files can be cleaned-up by running:
```
make clean-all
```
