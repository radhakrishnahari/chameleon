#!/bin/bash

# Performs an analysis of Chameleon source code:
# - we consider to be in Chameleon's source code root
# - we consider having the coverage file chameleon_coverage.xml in the root directory
# - we consider having cppcheck, rats, sonar-scanner programs available in the environment

# filter sources:
# - consider generated files in build
# - exclude base *z* files to avoid duplication
# - exclude cblas.h and lapacke-.h because not really part of chameleon and make cppcheck analysis too long
./tools/find_sources.sh

# Generate coverage xml output
INPUT_FILES=""
for name in $( ls -1 chameleon_*.lcov | grep -v simgrid)
do
    INPUT_FILES="$INPUT_FILES -a $name";
done
lcov $INPUT_FILES -o chameleon.lcov
lcov --summary chameleon.lcov
lcov_cobertura.py chameleon.lcov --output chameleon_coverage.xml

# Undefine this because not relevant in our configuration
export UNDEFINITIONS="-UCHAMELEON_USE_OPENCL -UWIN32 -UWIN64 -U_MSC_EXTENSIONS -U_MSC_VER -U__SUNPRO_C -U__SUNPRO_CC -U__sun -Usun -U__cplusplus"

# run cppcheck analysis
CPPCHECK_OPT=" -v -f --language=c --platform=unix64 --enable=all --xml --xml-version=2 --suppress=missingInclude ${UNDEFINITIONS}"
cppcheck $CPPCHECK_OPT --file-list=./filelist_none.txt 2> chameleon_cppcheck.xml
cppcheck $CPPCHECK_OPT -DPRECISION_s -UPRECISION_d -UPRECISION_c -UPRECISION_z --file-list=./filelist_s.txt 2>> chameleon_cppcheck.xml
cppcheck $CPPCHECK_OPT -DPRECISION_d -UPRECISION_s -UPRECISION_c -UPRECISION_z --file-list=./filelist_d.txt 2>> chameleon_cppcheck.xml
cppcheck $CPPCHECK_OPT -DPRECISION_c -UPRECISION_s -UPRECISION_d -UPRECISION_z --file-list=./filelist_c.txt 2>> chameleon_cppcheck.xml
cppcheck $CPPCHECK_OPT -DPRECISION_z -UPRECISION_s -UPRECISION_d -UPRECISION_c --file-list=./filelist_z.txt 2>> chameleon_cppcheck.xml

# create the sonarqube config file
cat > sonar-project.properties << EOF
sonar.host.url=https://sonarqube.inria.fr/sonarqube
sonar.login=$SONARQUBE_LOGIN

sonar.links.homepage=$CI_PROJECT_URL
sonar.links.scm=$CI_REPOSITORY_URL
sonar.links.ci=$CI_PROJECT_URL/pipelines
sonar.links.issue=$CI_PROJECT_URL/issues

sonar.projectKey=hiepacs:chameleon:gitlab:$CI_PROJECT_NAMESPACE:$CI_COMMIT_REF_NAME
sonar.projectDescription=Dense linear algebra subroutines for heterogeneous and distributed architectures
sonar.projectVersion=0.9

sonar.language=c
sonar.sources=build/runtime/openmp, build/runtime/parsec, build/runtime/quark, build/runtime/starpu, compute, control, coreblas, example, include, runtime, testing
sonar.inclusions=`cat filelist.txt | sed ':a;N;$!ba;s/\n/, /g'`
sonar.c.includeDirectories=$(echo | gcc -E -Wp,-v - 2>&1 | grep "^ " | tr '\n' ',').,$(find . -type f -name '*.h' | sed -r 's|/[^/]+$||' |sort |uniq | xargs echo | sed -e 's/ /,/g'),$PARSEC_DIR/include,$QUARK_DIR/include,$STARPU_DIR/include/starpu/1.3,$SIMGRID_DIR/include
sonar.sourceEncoding=UTF-8
sonar.c.errorRecoveryEnabled=true
sonar.c.gcc.charset=UTF-8
sonar.c.gcc.regex=(?<file>.*):(?<line>[0-9]+):[0-9]+:\\\x20warning:\\\x20(?<message>.*)\\\x20\\\[(?<id>.*)\\\]
sonar.c.gcc.reportPath=chameleon_build.log
sonar.c.coverage.reportPath=chameleon_coverage.xml
sonar.c.cppcheck.reportPath=chameleon_cppcheck.xml
sonar.c.clangsa.reportPath=build/analyzer_reports/*/*.plist
sonar.c.jsonCompilationDatabase=build/compile_commands.json, build/compile_commands_starpu_simgrid.json
EOF

# run sonar analysis + publish on sonarqube-dev
sonar-scanner -X > sonar.log
