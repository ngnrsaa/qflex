#!/usr/bin/env bash

if [[ $# -ne 1 ]]; then
  echo -e "\n\n\tUsage: $0 [-|folder]\n\n" >&2
  exit -1
fi

user_root=$1

get_location() {
  if whereis --version >/dev/null 2>/dev/null; then
    location=$(whereis $1)
    location=($location)
    echo ${location[1]}
  else
    echo $1
  fi
}

if $(get_location curl) -V >/dev/null 2>/dev/null; then
  echo "[OK] curl is installed."
else
  echo "[ERROR] curl is required."
  exit -1
fi

for cmd in tar git sed grep mktemp chroot unshare; do
  if $(get_location $cmd) --version >/dev/null 2>/dev/null; then
    echo "[OK] $cmd is installed."
  else
    echo "[ERROR] $cmd is required."
    exit -1
  fi
done

alpine_url="http://dl-cdn.alpinelinux.org/alpine/v3.10/releases/$(uname -m)/"
latest_miniroot=$(curl $alpine_url/latest-releases.yaml 2>/dev/null | grep 'file:' | grep miniroot | sed 's/ *file: *//g')

# Temporary folder
root=$(mktemp -d -t qflex-XXXXXXXXXXX)

# Create folder
echo "[CHROOR] Create folder $root." >&2
mkdir $root

# Get commands with absolute path
unshare="$(get_location unshare) -muipUCrf"
chroot="$(get_location chroot) $root/"

# Download alpine
echo "[CHROOT] Download $alpine_url/$latest_miniroot." >&2
curl $alpine_url/$latest_miniroot --output $root/rootfs.tar.gz

# Extract rootfs
echo "[CHROOT] Extract rootfs." >&2
tar xvf $root/rootfs.tar.gz -C $root >/dev/null

# Copy /etc/resolv.conf for internet access
echo "[CHROOT] Copy /etc/resolv.conf." >&2
cp -fv /etc/resolv.conf $root/etc/

# Clone qflex
echo "[CHROOT] Clone qFlex." >&2
git clone ../ $root/qflex

echo "[CHROOT] Update APK." >&2
$unshare $chroot /sbin/apk update

echo "[CHROOT] Install g++ make gsl-dev git" >&2
$unshare $chroot /sbin/apk add g++ make gsl-dev git

echo "[CHROOT] Create installation script." >&2
cat > $root/install_qflex.sh << EOF
#!/bin/sh

# Change folder
cd /qflex

# Update submodules
git submodule update --init --recursive

# Make qFlex
make -j$OMP_NUM_THREADS

# Change folder
cd tests

# Make tests
make -j$OMP_NUM_THREADS

# Create script to run tests
echo "#!/bin/sh" > /qflex/tests/run_all.sh
for file in /qflex/tests/*.x; do echo \$file; done >> /qflex/tests/run_all.sh
chmod 755 /qflex/tests/run_all.sh
EOF
chmod 755 $root/install_qflex.sh

echo "[CHROOT] Install qFlex." >&2
$unshare $chroot /install_qflex.sh

echo "[CHROOT] Run tests." >&2
$unshare $chroot /qflex/tests/run_all.sh

echo "[CHROOT] Container in: $root" >&2

if [[ x$user_root != "x-" ]]; then
  echo "[CHROOT] Moving container --> $user_root." >&2
  mv $root $user_root
fi
