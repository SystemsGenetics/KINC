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
code="${code}${nl}   if (strcmp(type.c_str(),\"${file}\")==0) ret = new ${file}(type,name);"
fi
done

echo "$header" > plugin.h
echo >> plugin.h
echo "DataPlugin* KINCPlugins::new_data(const std::string& type, const std::string& name)" >> plugin.h
echo "{" >> plugin.h
echo -n "   DataPlugin* ret = nullptr;" >> plugin.h
echo "$code" >> plugin.h
echo "   return ret;" >> plugin.h
echo "}" >> plugin.h

