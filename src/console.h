#ifndef CONSOLE_H
#define CONSOLE_H
#define __CL_ENABLE_EXCEPTIONS
#include <string>
#include <list>
#include <CL/cl.hpp>
#include "consolestream.h"



class Data;
class Analytic;



class Console
{
public:
   struct UnitTest;
private:
   enum command_e {gpu,glist,ginfo,gset,gclear,quit,error};
   struct {
      int i;
      int j;
      cl::Device* device;
   } _gpu;
   std::list<Data*> _data;
   cl::Device* gpu_get_device(std::string&,int&,int&);
   bool gpu_process(command_e,std::list<std::string>&);
   bool gpu_decode(std::list<std::string>&);
   bool process(command_e,std::list<std::string>&);
   bool decode(std::list<std::string>&);
   bool parse(std::string&);
   void get_terminal_line_loop();
public:
   ConsoleStream out;
   ConsoleStream warn;
   ConsoleStream err;
   Console();
   ~Console();
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
