#ifndef CONSOLE_H
#define CONSOLE_H
#include <string>
#include <list>
#include "consolestream.h"



class Data;
class Analytic;



class Console
{
private:
   std::list<Data*> _data;
public:
   ConsoleStream out;
   ConsoleStream warn;
   ConsoleStream err;
   Console();
   void run(int,char*[]);
   bool add(Data*,std::string&);
   bool del(std::string&);
   Data* find(std::string&);
   Data* new_data(std::string&);
   Analytic* new_analytic(std::string&);
};



extern Console g_console;



#endif
