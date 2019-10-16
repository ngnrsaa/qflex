# Install qFlex

qFlex can be installed locally using [GNU Autotools](https://www.gnu.org/software/automake/faq/autotools-faq.html). 
GNU Autotools automatically check pre-requisitions and install software to any
given target folder (qFlex uses `$HOME/local` as default folder).

## Required Software

To install qFlex, `autoconf` must be installed in the system.

## Check Pre-Requisites 

Pre-requisites can be checked by running:
```
$ autoconf && ./configure
```
`autoconf` generates a suitable `configure` file for the system and `configure`
will check for pre-requisites. If pre-requisites are satisfied, a `Makefile` is
created to compile and install qFlex. The default installation folder is
`$HOME/local`. If a different folder is desired, `--prefix` can be passed as an
argument to `configure`:
```
$ autoconf && ./configure --prefix=/new/installation/folder/
```
By default, `g++` is used to compile qFlex. To use a different compiler (as
Intel `icpc`), the user may provide the bash variable `CXX`:
```
$ autoconf && CXX=icpc ./configure
```
When `CXX=icpc` is used, qFlex will try to use the MKL library and return an
error if not available. For further info:
```
$ autoconf && ./configure --help
```

## Compile qFlex

To compile qFlex, it suffices to run:
```
make && make run-tests
```
The first command will compile a qFlex executable (`src/qflex.x`) that can be
directly used without installation. The second command will run all the
available tests. While not requires, it is suggested to make sure that qFlex has
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
["Check Pre-Requisites "](#check-pre-requisites) to change the default folder.

## Clean-up

All the installation files can be cleaned-up by running:
```
make clean-all
```
