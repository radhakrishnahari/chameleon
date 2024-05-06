#!/usr/bin/env bash
#
# @file run.sh
#
# @copyright 2020-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                      Univ. Bordeaux. All rights reserved.
#
# @version 1.2.0
# @author Florent Pruvost
# @date 2022-02-22
#
set -x

# Unset the binding environment of the CI for this specific case
unset STARPU_MPI_NOBIND
unset STARPU_WORKERS_NOBIND

# to avoid a lock during fetching chameleon branch in parallel
export XDG_CACHE_HOME=/tmp/guix-$$

# save guix commits
#guix describe --format=json > guix.json
guix time-machine -C ./tools/bench/guix-channels.scm -- describe --format=json > guix.json

# define env var depending on the node type
if [[ "$NODE" == "bora" ]]; then
  export SLURM_CONSTRAINTS="bora,omnipath"
  export CHAMELEON_BUILD_OPTIONS="-DCHAMELEON_USE_MPI=ON -DCMAKE_BUILD_TYPE=Release"
  export STARPU_HOSTNAME="bora"
elif [[ "$NODE" == "sirocco" ]]; then
  export SLURM_CONSTRAINTS="sirocco,v100"
  export CHAMELEON_BUILD_OPTIONS="-DCHAMELEON_USE_MPI=ON -DCHAMELEON_USE_CUDA=ON -DCMAKE_BUILD_TYPE=Release"
  export STARPU_HOSTNAME="sirocco"
  export LD_PRELOAD="/usr/lib64/libcuda.so"
else
  echo "$0: Please set the NODE environment variable to bora or sirocco."
  exit -1
fi

# define env var and guix rule to use depending on the mpi vendor
GUIX_ENV="chameleon"
export MPI_OPTIONS=""
if [[ "$MPI" == "openmpi" ]]; then
  export MPI_OPTIONS="--bind-to board"
  GUIX_ENV_MPI=""
  GUIX_ADHOC_MPI=""
  if [[ "$NODE" == "sirocco" ]]; then
    GUIX_ENV="chameleon-cuda"
  fi
elif [[ "$MPI" == "nmad" ]]; then
  export MPI_OPTIONS="-DPIOM_DEDICATED=1 -DPIOM_DEDICATED_WAIT=1 hwloc-bind --cpubind machine:0"
  GUIX_ENV_MPI="--with-input=openmpi=nmad"
  GUIX_ADHOC_MPI="which gzip zlib tar inetutils util-linux procps nmad"
else
  echo "$0: Please set the MPI environnement variable to openmpi or nmad."
  exit -1
fi
GUIX_ADHOC="coreutils gawk grep hwloc jube nss-certs openssh perl python python-click python-certifi python-elasticsearch python-gitpython python-matplotlib python-pandas python-seaborn r-ggplot2 r-plyr r-reshape2 sed slurm mkl@2019"
GUIX_RULE="-D $GUIX_ENV $GUIX_ENV_MPI $GUIX_ADHOC $GUIX_ADHOC_MPI"

# Submit jobs

#exec guix shell --pure \
exec guix time-machine -C ./tools/bench/guix-channels.scm -- shell --pure \
       --preserve=PLATFORM \
       --preserve=NODE \
       --preserve=LD_PRELOAD \
       --preserve=^CI \
       --preserve=proxy$ \
       --preserve=^SLURM \
       --preserve=^JUBE \
       --preserve=^MPI \
       --preserve=^STARPU \
       --preserve=^CHAMELEON \
       $GUIX_RULE \
       -- /bin/bash --norc ./tools/bench/plafrim/slurm.sh
err=$?

# clean tmp
rm -rf /tmp/guix-$$

# exit with error code from the guix command
exit $err
