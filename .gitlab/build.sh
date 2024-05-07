#!/usr/bin/env bash

fatal() {
    echo "$0: error occurred, exit"
    exit 1
}

set -x

SCAN=""

if [[ "$SYSTEM" != "windows" ]]; then
  if [[ "$SYSTEM" == "macosx" ]]; then
    # check starpu is already installed and install it if necessary
    DEP_INSTALLED=`brew ls --versions starpu | cut -d " " -f 2`
    if [[ -z "${DEP_INSTALLED}" ]]; then
      # dep is not installed, we have to install it
      brew install --build-from-source ./tools/homebrew/starpu.rb
    else
      # dep is already installed, check the version with our requirement
      DEP_REQUIRED=`brew info --json ./tools/homebrew/starpu.rb |jq -r '.[0].versions.stable'`
      if [[ "${DEP_INSTALLED}" != "${DEP_REQUIRED}" ]]; then
        # if the installed version is not the required one, re-install
        brew remove --force --ignore-dependencies starpu
        brew install --build-from-source ./tools/homebrew/starpu.rb
      fi
    fi
    # clang is used on macosx and it is not compatible with MORSE_ENABLE_COVERAGE=ON
    # to avoid the Accelerate framework and get Openblas we use BLA_PREFER_PKGCONFIG
    cmake -B build-${VERSION} -S . \
          -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=$PWD/install-${VERSION} \
          -DMORSE_ENABLE_COVERAGE=OFF -DBLA_PREFER_PKGCONFIG=ON || fatal
  else
    source .gitlab-ci-env.sh $CHAM_CI_ENV_ARG || fatal
    if [[ $CI_COMMIT_REF_NAME == $CI_DEFAULT_BRANCH ]]
    then
      SCAN="scan-build -plist --intercept-first --exclude CMakeFiles --analyze-headers -o analyzer_reports "
    fi
    eval '${SCAN}cmake -B build-${VERSION} -S . -C cmake_modules/gitlab-ci-initial-cache.cmake $BUILD_OPTIONS' || fatal
  fi
else
  # on windows the mpi_f08 interface is missing, see https://www.scivision.dev/windows-mpi-msys2/
  # do not use static libraries because executables are too large and the build
  # directory can reach more than 10Go
  cmake -GNinja -B build-${VERSION} -S . \
        -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=$PWD/install-${VERSION} \
        -DCHAMELEON_USE_MPI=OFF \
        || fatal
fi
eval '${SCAN}cmake --build build-${VERSION} -j 4' || fatal
cmake --install build-${VERSION} || fatal
rm -r install-${VERSION} || fatal
