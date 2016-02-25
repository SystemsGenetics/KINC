#ifndef PLUGINS_H
#define PLUGINS_H
#include <string>
#include "../analytic.h"
#include "../dataplugin.h"



namespace KINCPlugins
{
   Analytic* new_analytic(const std::string&);
   DataPlugin* new_data(const std::string&,const std::string&);
}



#endif
