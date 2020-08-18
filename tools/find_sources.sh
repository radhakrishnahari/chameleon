#!/bin/sh
#set -x

SRCDIR_TO_ANALYZE="build/runtime/openmp build/runtime/parsec build/runtime/quark build/runtime/starpu build compute control coreblas example include runtime testing"

echo $PWD
rm -f filelist.txt
for dir in ${SRCDIR_TO_ANALYZE}
do
    find $dir -name '*\.[ch]' >> filelist.txt
done

# Remove all CMakeFiles generated file
sed -i '/CMakeFiles/d' filelist.txt

# Remove installed files
sed -i '/build.*\/install.*/d' filelist.txt

# Remove original files used for precision generation
for file in `git grep "@precisions" | awk -F ":" '{ print $1 }'`
do
    sed -i "\:^$file.*:d" filelist.txt
done

# Remove external header files
for file in coreblas/include/coreblas/cblas.h coreblas/include/coreblas/lapacke.h coreblas/include/coreblas/lapacke_config.h coreblas/include/coreblas/lapacke_mangling.h
do
    sed -i "\:^$file.*:d" filelist.txt
done

rm -f filelist_*.txt
for name in $(cat filelist.txt)
do
    test=$(grep "@generated" $name | wc -l)
    if [ $test -gt 0 ]
    then
        prec=$(grep "@generated" $name | sed 's/^.*[scdz] -> \([sdcz]\).*$/\1/')
        echo $name >> filelist_${prec}.txt
    else
        echo $name >> filelist_none.txt
    fi
done
