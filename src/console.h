/*
 * See GPL.txt for copyright information.
 *
 * Author: Joshua Burns
 *
 */
#ifndef CONSOLE_H
#define CONSOLE_H
#define __CL_ENABLE_EXCEPTIONS
#include <string>
#include <list>
#include <CL/cl.hpp>
#include "terminal.h"
#include "cldevice.h"



class Data;
class Analytic;



/*
 * NOTE: for new_data and new_analytic see plugin.cpp for source code and
 * comments.
 */
namespace KINCPlugins
{
   extern Data* new_kinc_data(std::string&);
   extern Analytic* new_kinc_analytic(std::string&);
}



/// @brief Takes over execution of program and controls data/analytic objects.
///
/// Designed to take over execution of the program and manage commands supplied
/// by the user and all data and analytic objects. The run() function will take
/// over execution and release when the user quits the program. Only one
/// instance of this class is allowed to exist for the entire program.
///
/// @pre Only one instance of this class can exist at any one time.
class Console
{
   struct UnitTest;
public:
   // ****************************** Basic Methods **************************
   Console(const Console&) = delete;
   Console(Console&&) = delete;
   Console& operator=(const Console&) = delete;
   Console& operator=(Console&&) = delete;
   Console(int argc, char* argv[], Terminal& tm);
   ~Console();
   // ****************************** Functions ******************************
   void run();
   //bool add(Data*,std::string&);
   //bool del(std::string&);
   //Data* find(std::string&);
private:
   // ****************************** Enumerations ***************************
   enum class Command {gpu,quit,error};
   enum class GpuCommand { list,info,set,clear,error };
   // ****************************** Functions ******************************
   void terminal_loop();
   bool parse(std::string&);
   bool decode(std::list<std::string>&);
   bool process(Command,std::list<std::string>&);
   bool gpu_decode(std::list<std::string>&);
   bool gpu_process(GpuCommand,std::list<std::string>&);
   void gpu_list();
   bool gpu_info(std::list<std::string>&);
   bool gpu_set(std::list<std::string>&);
   void gpu_clear();
   // ****************************** Static Variables ***********************
   static bool _lock;
   // ****************************** Variables ******************************
   Terminal& _tm;
   CLDevice* _device;
   //std::list<Data*> _data;
   bool _alive;
};



#endif
