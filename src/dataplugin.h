#ifndef DATAPLUGIN_H
#define DATAPLUGIN_H
#include "data.h"
#include "kincfile.h"



class DataPlugin : public Data, public KincFile
{};



#endif
