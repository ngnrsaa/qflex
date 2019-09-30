#!/bin/bash

OMP_NUM_THREADS=4

get_location() {
  location=$(whereis $1)
  location=($location)
  echo ${location[1]}
}

# Check if whereis exists
if whereis whereis >/dev/null 2>/dev/null; then
  echo "[OK] Whereis is installed."
else
  echo "[ERROR] Whereis is required."
  exit -1
fi

# Check if unshare exists and that user namespaces are enabled
if $(get_location unshare) -r echo >/dev/null 2>/dev/null; then
  echo "[OK] User Namespaces are enabled." >&2
else
  echo "[ERROR] User Namespaces are not enabled." >&2
  exit -1
fi

# Check if curl exists
if $(get_location curl) -V >/dev/null 2>/dev/null; then
  echo "[OK] Curl is installed."
else
  echo "[ERROR] Curl is required."
  exit -1
fi

# Check if tar exists
if $(get_location tar) --version >/dev/null 2>/dev/null; then
  echo "[OK] Tar is installed."
else
  echo "[ERROR] Tar is required."
  exit -1
fi

# Check if chroot exists
if $(get_location chroot) --version >/dev/null 2>/dev/null; then
  echo "[OK] Chroot is installed."
else
  echo "[ERROR] Chroot is required."
  exit -1
fi

alpine_url="http://dl-cdn.alpinelinux.org/alpine/v3.10/releases/$(uname -m)/"
latest_miniroot=$(curl $alpine_url/latest-releases.yaml 2>/dev/null | grep 'file:' | grep miniroot | sed 's/ *file: *//g')

# Download alpine
curl $alpine_url/$latest_miniroot --output rootfs.tar.gz

# Create container
#root=/tmp/qflex_rootless_${RANDOM}
root=/tmp/qflex_rootless_123

echo "[CHROOT] $root" >&2

# Get commands with absolute path
unshare="$(get_location unshare) -rf"
chroot="$(get_location chroot) $root/"

# Extract rootfs
mkdir $root
tar xvf rootfs.tar.gz -C $root >/dev/null

# Copy /etc/resolv.conf for internet access
cp -fv /etc/resolv.conf $root/etc/

# Copy qflex
echo "[CHROOT] Copy qFlex." >&2
mkdir $root/qflex
cp -r ../ $root/qflex

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

echo "[CHROOT] $root" >&2
