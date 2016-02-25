/*
 * See GPL.txt for copyright information.
 *
 * Author: Joshua Burns
 *
 */
#ifndef CONSOLE_H
#define CONSOLE_H
#include <string>
#include <list>
#include "terminal.h"
#include "cldevlist.h"
#include "datamap.h"
#include "exception.h"
#include "getopts.h"



class CLDevice;



/// @brief Main program Console.
///
/// Designed to take over execution of the program and manage commands supplied
/// by the user and all data and analytic objects. The run() function will take
/// over execution and release when the user quits the program. Only one
/// instance of this class is allowed to exist for the entire program.
///
/// @warning Only one instance of this class can exist at any one time.
///
/// @author Josh Burns
/// @date 22 Jan 2016
class Console
{
public:
   // *
   // * EXCEPTIONS
   // *
   struct Exception;
   struct InvalidUse;
   // *
   // * BASIC METHODS
   // *
   Console(const Console&) = delete;
   Console(Console&&) = delete;
   Console& operator=(const Console&) = delete;
   Console& operator=(Console&&) = delete;
   Console(int,char*[],Terminal&,DataMap&);
   ~Console();
   // *
   // * FUNCTIONS
   // *
   void run();
   //bool add(Data*,std::string&);
   //bool del(std::string&);
   //Data* find(std::string&);
private:
   using string = std::string;
   struct CommandError
   {
      CommandError(const char* c, const std::string& m): cmd {c}, msg {m} {}
      const char* cmd;
      std::string msg;
   };
   struct CommandQuit {};
   // *
   // * FUNCTIONS
   // *
   void terminal_loop();
   void process(GetOpts&);
   void gpu_process(GetOpts&);
   void gpu_list();
   void gpu_info(GetOpts&);
   void gpu_set(GetOpts&);
   void gpu_clear();
   /*void data_open(slist&);
   void data_load(slist&);
   void data_dump(slist&);
   void data_query(slist&);
   void data_close(slist&);
   void data_list();
   void analytic(slist&);
   DataPlugin* find_data(const string&);
   void parse_data_options(DataPlugin*,slist&);
   void parse_analytic_inputs(aptr&,const string&);
   void parse_analytic_outputs(aptr&,const string&,dlist&);
   DataPlugin* parse_analytic_ndata(const string&,string&);
   void parse_analytic_options(aptr&,slist&);*/
   // *
   // * STATIC VARIABLES
   // *
   static bool _lock;
   // *
   // * VARIABLES
   // *
   /// Reference to program's main Terminal interface.
   Terminal& _tm;
   //
   DataMap& _dataMap;
   /// Pointer to OpenCL device that is set for computation acceleration. By
   /// or if clear command issues this is set to nullptr.
   CLDevice* _device;
   /// List of all possible OpenCL devices on program's computer.
   CLDevList _devList;
   //std::list<Data*> _data;
   /// Boolean variable that is used as state information stating if the console
   /// is still active and should continue accepting input from user. Once this
   /// is false the run() command will exit and return control.
   bool _alive;
};



/// Base exception class for any exception thrown from Console class.
struct Console::Exception : public ::Exception
{
   using ::Exception::Exception;
};

/// Exception thrown if console class is used in an invalid way.
struct Console::InvalidUse : public Console::Exception
{
   InvalidUse(const char* file, int line):
      Exception(file,line,"Terminal::InvalidUse")
   {}
};



#endif
