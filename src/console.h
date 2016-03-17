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
/// by the user for all data and analytic objects. Also manages and displays
/// information about all OpenCL devices attached to the program's machine.
///
/// @warning Only one instance of this class can exist at any one time.
///
/// @author Josh Burns
/// @date 17 Mar 2016
class Console
{
public:
   // *
   // * EXCEPTIONS
   // *
   struct Exception;
   struct InvalidUse;
   // *
   // * DECLERATIONS
   // *
   class CommandError;
   // *
   // * BASIC METHODS
   // *
   Console(int,char*[],Terminal&,DataMap&);
   ~Console();
   // *
   // * COPY METHODS
   // *
   Console(const Console&) = delete;
   Console& operator=(const Console&) = delete;
   // *
   // * MOVE METHODS
   // *
   Console(Console&&) = delete;
   Console& operator=(Console&&) = delete;
   // *
   // * FUNCTIONS
   // *
   void run();
private:
   using string = std::string;
   using hiter = History::Iterator;
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
   void data_open(GetOpts&);
   void data_close(GetOpts&);
   void data_select(GetOpts&);
   void data_clear();
   void data_list();
   void data_history(GetOpts&);
   void data_load(GetOpts&);
   void data_dump(GetOpts&);
   void data_query(GetOpts&);
   void analytic(GetOpts&);
   void seperate(const string&,const string&,string&,string&);
   void rec_history(hiter,hiter,int);
   void print_pad(int);
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



class Console::CommandError
{
public:
   using string = std::string;
   CommandError(const string&,const string&);
   void print(Terminal&);
private:
   string _who;
   string _msg;
};



inline Console::CommandError::CommandError(const string& who,
                                           const string& msg):
   _who(who),
   _msg(msg)
{}



inline void Console::CommandError::print(Terminal& tm)
{
   tm << _who << ": " << _msg << Terminal::endl;
}



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
