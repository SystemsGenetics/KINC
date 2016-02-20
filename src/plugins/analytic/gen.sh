#!/bin/bash

nl='
'
header="";
code="";

for file in *
do
if [ -d $file ]
then
header="${header}${nl}#include \"${file}/${file}.h\""
code="${code}${nl}   if (strcmp(type.c_str(),\"${file}\")==0) ret = new ${file};"
fi
done

echo "$header" > plugin.h
echo >> plugin.h
echo "Analytic* KINCPlugins::new_analytic(const std::string& type)" >> plugin.h
echo "{" >> plugin.h
echo -n "   Analytic* ret = nullptr;" >> plugin.h
echo "$code" >> plugin.h
echo "   return ret;" >> plugin.h
echo "}" >> plugin.h

