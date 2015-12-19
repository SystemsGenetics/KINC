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
   struct {
      int platform;
      int device;
   } _gpu;
   std::list<Data*> _data;
   bool gpu(std::list<std::string>&);
   bool process(std::string&);
   void command();
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
   void touch_output();
};



extern Console g_console;



#endif
