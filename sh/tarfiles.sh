#!/bin/bash

# Srcipt to tar instalation files of imodel

date=` date +%F `
version=` date +%y.%m.%d `
echo "Today: " $date

sourcefiles="src/*.f90"

parfiles="par/*.par "

scripts="sh/*.sh "

others="Makefile \
README.*"

files="$sourcefiles $parfiles $scripts $others"

output="fvcs$version.tar.bz2"

tar cjfv $output $files

echo "File " $output " ready!"
echo

echo "-------------------------------------------"
