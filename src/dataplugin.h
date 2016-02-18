#ifndef DATAPLUGIN_H
#define DATAPLUGIN_H
#include <string>
#include "data.h"
#include "kincfile.h"



class DataPlugin : public Data, public KincFile
{
   inline DataPlugin(const std::string&,const std::string&);
   std::string type();
private:
   std::string _type;
};



inline DataPlugin::DataPlugin(const std::string& type,
                              const std::string& fileName):
   _type(type),
   KincFile(fileName)
{}



#endif
