#!/usr/bin/env bash

set -e
set -x

cd build-$VERSION
source ../.gitlab-ci-env.sh $CHAM_CI_ENV_ARG
ctest --output-on-failure --no-compress-output $TESTS_RESTRICTION --output-junit ../${LOGNAME}-junit.xml
if [[ "$SYSTEM" == "linux" ]]; then
  # clang is used on macosx and it is not compatible with MORSE_ENABLE_COVERAGE=ON
  # so that we can only make the coverage report on the linux runner with gcc
  cd ..
  lcov --directory build-$VERSION --capture --output-file ${LOGNAME}.lcov
fi
