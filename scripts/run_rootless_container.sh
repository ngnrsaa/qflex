#!/usr/bin/env bash

# See "docs/rootless-container.md" for more information about this file.

# Get Script path
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

print_help() {
  echo -e "\n\tUsage: $(basename $0) folder [-h]\n"                   >&2
  echo -e "\t\tRun qflex rootless container."                         >&2
  echo                                                                >&2
  echo -e "\tOptions:"                                                >&2
  echo -e "\t\t-h   : Print this help."                               >&2
  echo                                                                >&2
}

num_args=$#

# Check that args are given
if [[ $num_args == 0 ]]; then
  print_help
  exit -1
fi

# Parse options
for((idx=1; idx<=$num_args; ++idx)); do
  if [[ $idx == 1 ]]; then
    if [[ $idx == 1 && ${1:0:1} == "-" ]]; then
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

if [[ ! -d $user_root ]]; then
  echo "No such directory!" >&2
  print_help
  exit -1
fi

get_location() {
  if whereis --version >/dev/null 2>/dev/null; then
    location=$(whereis $1)
    location=($location)
    echo ${location[1]}
  else
    echo $1
  fi
}

for cmd in env chroot unshare; do
  if $(get_location $cmd) --version >/dev/null 2>/dev/null; then
    echo "[OK] $cmd is installed."
  else
    echo "[ERROR] $cmd is required."
    exit -1
  fi
done

# Check if unshare can be run
unshare -r echo >/dev/null 2>/dev/null
if [[ $? != 0 ]]; then
  echo "[ERROR] Not enough privilegies to run unshare. Use: sudo echo 1 > /proc/sys/kernel/unprivileged_userns_clone" >&2
  exit -1
fi

# Get commands with absolute path
unshare="$(get_location unshare) -muipUCrf"
chroot="$(get_location chroot) $user_root/ $(get_location env) -i PATH=/bin/:/sbin:/usr/bin/:/usr/sbin/:/usr/local/bin:/usr/local/sbin OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} LANG=${LANG:-en}"

$unshare $chroot /bin/sh
