#ifndef DATAPLUGIN_H
#define DATAPLUGIN_H
#include <string>
#include "data.h"
#include "kincfile.h"



class DataPlugin : public Data, public KincFile
{
public:
   inline std::string type();
protected:
   inline DataPlugin(const std::string&,const std::string&);
private:
   std::string _type;
};



inline std::string DataPlugin::type()
{
   return _type;
}



inline DataPlugin::DataPlugin(const std::string& type,
                              const std::string& fileName):
   _type(type),
   KincFile(fileName)
{}



#endif
