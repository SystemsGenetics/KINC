#!/bin/bash

nl='
'
header="#include \"../plugins.h\"";
code="";

for file in *
do
if [ -d $file ]
then
header="${header}${nl}#include \"${file}/${file}.h\""
code="${code}${nl}   if (strcmp(type.c_str(),\"${file}\")==0) ret = new ${file}(type,name);"
fi
done

echo "$header" > plugin.cpp
echo >> plugin.cpp
echo "DataPlugin* KINCPlugins::new_data(const std::string& type, const std::string& name)" >> plugin.cpp
echo "{" >> plugin.cpp
echo -n "   DataPlugin* ret = nullptr;" >> plugin.cpp
echo "$code" >> plugin.cpp
echo "   return ret;" >> plugin.cpp
echo "}" >> plugin.cpp

