#!/bin/bash

echo "const char* `echo $1 | sed 's/\./_/g'` = " > ${1}.h
cat $1 | sed 's/$/\\n"/g' | sed 's/^/"/g' >> ${1}.h
echo '"";' >> ${1}.h
