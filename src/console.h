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
/// @date 18 March 2016
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
   /// Used to determine if an instance of this class exists or not.
   static bool _lock;
   // *
   // * VARIABLES
   // *
   /// Reference to program's main Terminal interface.
   Terminal& _tm;
   /// Reference to program's data list manager.
   DataMap& _dataMap;
   /// Pointer to OpenCL device that is set for computation acceleration. Set to
   /// nullptr if no device is selected.
   CLDevice* _device;
   /// List of all possible OpenCL devices on program's machine.
   CLDevList _devList;
};



/// @brief Single console error.
///
/// This holds a single error that occured and is thrown from somewhere outside
/// of the console class which the console class is designed to catch and
/// report to the user.
///
/// @author Josh Burns
/// @date 18 March 2016
class Console::CommandError
{
public:
   using string = std::string;
   CommandError(const string&,const string&);
   void print(Terminal&);
private:
   /// Holds who threw the command error.
   string _who;
   /// Holds the message of the command error.
   string _msg;
};



/// Initializes a new command error and sets who threw it and its message.
///
/// @param who Identifies who threw the error message.
/// @param msg The message describing the error that occured.
inline Console::CommandError::CommandError(const string& who,
                                           const string& msg):
   _who(who),
   _msg(msg)
{}



/// Prints the error message to the terminal.
///
/// @param tm The program's terminal that will be printed to.
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
