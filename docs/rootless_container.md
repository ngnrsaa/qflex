# Rootless Containers

As an alternative to [Docker](https://docker.com), qFlex can be run in rootless
containers that do not require `root` privileges.

## Build Rootless Image

To build a rootless qFlex image and run all tests in it, run:

```
$ bash scripts/create_rootless_container.sh (-|image_folder) [-h -j <p> -x -r -c]
```
The script will create a container in `image_folder`, or in a temporary
directory if `-` is specified. The script also supports the following options:
```
-h        Print this help.
-j <p>    Number of parallel processes (where possible).
-x        Do not install qFlex (just create container).
-r        Run rootless container immediately after creation.
-c        Disable installation of Cirq in the rootless container.
```

## Namespaces

To create a rootless container, the user must be allowed to unshare namespaces
from parent. If not enabled:
```
$ sudo echo 1 > /proc/sys/kernel/unprivileged_userns_clone
```

## Run Rootless Image

To run the rootless image:
```
$ bash scripts/run_rootless_container.sh image_folder
```

Inside the rootless image, qFlex can be installed by executing:
```
/install_qflex.sh
```
After the installation of qFlex, tests can be run as:
```
make -C /qflex/ run-tests
```

To run a sample simulation inside the image, use the following command:
```
/qflex/src/qflex.x /qflex/config/circuits/bristlecone_48_1-24-1_0.txt \
/qflex/config/ordering/bristlecone_48.txt /qflex/config/grid/bristlecone_48.txt
```
