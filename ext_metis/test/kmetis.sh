#!/bin/bash
#
#  Get INCLUDE files.
#
cp ../Lib/*.h .
cp ../Lib/*.a .
#
gcc -c -g -I. kmetis.c 
if [ $? -ne 0 ]; then
  echo "Errors compiling kmetis.c."
  exit
fi
#
gcc -o kmetis kmetis.o libmetis_x86_64.a -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading the object files."
  exit
fi
#
