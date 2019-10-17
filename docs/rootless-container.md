# Rootless Containers

Alternatively to [Docker](https://docker.com), qFlex can be run in rootless
containers that do not require `root` privileges.

## Build Rootless Image

To build a rootless qFlex image:

```
$ bash scripts/create-rootless-container.sh (-|image_folder) [-x -c]
```
The script will create a container in `image_folder`, or in a temporary
directory if `-` is specified. The flag `-x` create an rootless image without
installing qFlex while the flag `-c` instructs the script to install Cirq in the
rootless container.

## Namespaces

To create a rootless container, the user must be allowed to unshare namespaces
from parent. If not enabled:
```
$ sudo echo 1 > /proc/sys/kernel/unprivileged_userns_clone
```

## Run Rootless Image

To run the rootless image:
```
$ bash scripts/run-rootless-image.sh image_folder
```

Inside the rootless image, qFlex and Cirq can be installed by executing:
```
/install_qflex.sh
```
and
```
/install_cirq.sh
```
respectively. After the installation of qFlex, tests can be run as:
```
make -C /qflex/ run-tests
```
