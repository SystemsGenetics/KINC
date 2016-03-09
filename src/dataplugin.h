#ifndef DATAPLUGIN_H
#define DATAPLUGIN_H
#include <string>
#include "data.h"
#include "kincfile.h"



class DataPlugin : public Data, public KincFile
{
public:
   using string = std::string;
   string type();
protected:
   DataPlugin(const string&,const string&);
private:
   string _type;
};



inline DataPlugin::string DataPlugin::type()
{
   return _type;
}



inline DataPlugin::DataPlugin(const string& type, const string& fileName):
   _type(type),
   KincFile(fileName)
{}



#endif
