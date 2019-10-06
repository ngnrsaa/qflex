#!/usr/bin/env bash

# Get Script path
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

print_help() {
  echo -e "\n\tUsage: $(basename $0) {-|folder} [-h -x]\n"            >&2
  echo -e "\t\tCreate rootless container for qflex in [folder]."      >&2
  echo -e "\t\tIf [folder] == \"-\", use a temporary folder instead." >&2
  echo                                                                >&2
  echo -e "\tOptions:"                                                >&2
  echo -e "\t\t-h   : Print this help."                               >&2
  echo -e "\t\t-x   : Do not run tests (just create container)."      >&2
  echo                                                                >&2
}

get_location() {
  if whereis --version >/dev/null 2>/dev/null; then
    location=$(whereis $1)
    location=($location)
    echo ${location[1]}
  else
    echo $1
  fi
}

user_root=""
no_tests=""
num_args=$#

# Check that args are given
if [[ $num_args == 0 ]]; then
  print_help
  exit -1
fi

# Parse options
for((idx=1; idx<=$num_args; ++idx)); do
  if [[ $idx == 1 ]]; then
    if [[ $idx == 1 && ${1:0:1} == "-" && $1 != "-" ]]; then
      print_help
      exit -1
    fi
    user_root=$1
  else
    if [[ ${1:0:1} == "-" ]]; then
      case $1 in
        -h)
          print_help
          exit -1
          ;;
        -x)
          no_tests=1
          ;;
        *)
          print_help
          exit -1
          ;;
      esac
    else
      print_help
      exit -1
    fi
  fi
  shift
done

# Check if folder exists
if [[ $user_root != "-" ]]; then

  if [[ ! -w $(realpath $(dirname $user_root)) ]]; then
    echo "Directory $(realpath $(dirname $user_root)) does not exist or not writable by you." >&2
    print_help
    exit -1
  fi
  
  # Check if root folder does not exists
  if [[ -d $(realpath $user_root) ]]; then
    echo "Directory $(realpath $user_root) already exists." >&2
    print_help
    exit -1
  fi

fi

if $(get_location curl) -V >/dev/null 2>/dev/null; then
  echo "[OK] curl is installed." >&2
  curl="$(get_location curl) --output"
else
  if $(get_location wget) --version >/dev/null 2>/dev/null; then
    echo "[OK] wget is installed." >&2
    curl="$(get_location wget) -O"
  else
    echo "[ERROR] Either curl or wget is required." >&2
    exit -1
  fi
fi

for cmd in tar git sed grep mktemp chroot unshare; do
  if $(get_location $cmd) --version >/dev/null 2>/dev/null; then
    echo "[OK] $cmd is installed." >&2
  else
    echo "[ERROR] $cmd is required." >&2
    exit -1
  fi
done

# Check if unshare can be run
unshare -r echo >/dev/null 2>/dev/null
if [[ $? != 0 ]]; then
  echo "[ERROR] Not enough privilegies to run unshare. Use: sudo echo 1 > /proc/sys/kernel/unprivileged_userns_clone" >&2
  exit -1
fi

alpine_url="http://dl-cdn.alpinelinux.org/alpine/v3.10/releases/$(uname -m)/"
latest_miniroot=$($curl - $alpine_url/latest-releases.yaml 2>/dev/null | grep 'file:' | grep miniroot | sed 's/ *file: *//g')

# Create temporary folder if the user does not provide a target folder
if [[ $user_root == "-" ]]; then
  echo "[CHROOT] Create temporary folder $root." >&2
  root=$(mktemp -d -t qflex-XXXXXXXXXXX)
else
  mkdir $user_root
  root=$user_root
fi

# Get commands with absolute path
unshare="$(get_location unshare) -muipUCrf"
chroot="$(get_location chroot) $root/ $(get_location env) -i PATH=/bin/:/sbin:/usr/bin/:/usr/sbin/:/usr/local/bin:/usr/local/sbin OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}"

# Download alpine
echo "[CHROOT] Download $alpine_url/$latest_miniroot." >&2
$curl $root/rootfs.tar.gz $alpine_url/$latest_miniroot

# Extract rootfs
echo "[CHROOT] Extract rootfs." >&2
tar xvf $root/rootfs.tar.gz -C $root >/dev/null

# Copy /etc/resolv.conf for internet access
echo "[CHROOT] Copy /etc/resolv.conf." >&2
cp -fv /etc/resolv.conf $root/etc/

# Clone qflex
echo "[CHROOT] Clone qFlex." >&2
git clone $SCRIPTPATH/../ $root/qflex

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
make -j\$OMP_NUM_THREADS
EOF

if [[ -z $no_tests || $no_tests != "1" ]]; then
cat >> $root/install_qflex.sh << EOF

# Change folder
cd tests

# Make tests
make -j\$OMP_NUM_THREADS

# Create script to run tests
echo "#!/bin/sh" > /qflex/tests/run_all.sh
COUNT=\$(ls -1 /qflex/tests/*.x | wc -l)
IDX=1
for file in /qflex/tests/*.x; do
  if [[ \$IDX -lt \$COUNT ]]; then
    if [[ \$IDX -gt 1 ]]; then
      echo -ne "\t"
    fi
    echo "\$file && \\ " 
  else
    echo -e "\t\$file"
  fi
  IDX=\$(echo \$IDX+1 | bc -l)
done >> /qflex/tests/run_all.sh
chmod 755 /qflex/tests/run_all.sh
EOF
fi
chmod 755 $root/install_qflex.sh

echo "[CHROOT] Install qFlex." >&2
$unshare $chroot /install_qflex.sh

if [[ -z $no_tests || $no_tests != "1" ]]; then
  echo "[CHROOT] Run tests." >&2
  $unshare $chroot /qflex/tests/run_all.sh
fi

echo "[CHROOT] Container in: $(realpath $root)" >&2
