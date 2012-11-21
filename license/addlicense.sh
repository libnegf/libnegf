#!/bin/bash  
for x in $*; do  
	head -$LICENSELEN $x | diff libnegf.txt - || ( ( cat libnegf.txt; echo; cat $x) > /tmp/file;  
	mv /tmp/file $x )  
done  
