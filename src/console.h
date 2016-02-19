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
   struct CommandError
   {
      CommandError(const char* c, const std::string& m): cmd {c}, msg {m} {}
      const char* cmd;
      std::string msg;
   };
   struct CommandQuit {};
   // *
   // * ENUMERATIONS
   // *
   ///
   /// Main commands.
   enum class Command
   {
      gpu, ///< This is an OpenCL command.
      open,
      load,
      dump,
      query,
      close,
      list,
      analytic,
      quit, ///< The quit command.
      error ///< Error at parsing command.
   };
   ///
   /// OpenCL subcommands.
   enum class GpuCommand
   {
      list, ///< The list subcommand.
      info, ///< The info subcommand.
      set, ///< The set subcommand.
      clear ///< The clear subcommand.
   };
   // *
   // * FUNCTIONS
   // *
   void terminal_loop();
   void parse(std::string&);
   void decode(std::list<std::string>&);
   void process(Command,std::list<std::string>&);
   void gpu_decode(std::list<std::string>&);
   void gpu_process(GpuCommand,std::list<std::string>&);
   void gpu_list();
   void gpu_info(std::list<std::string>&);
   void gpu_set(std::list<std::string>&);
   void gpu_clear();
   void data_open(std::list<std::string>&);
   void data_load(std::list<std::string>&);
   void data_dump(std::list<std::string>&);
   void data_query(std::list<std::string>&);
   void data_close(std::list<std::string>&);
   void data_list(std::list<std::string>&);
   void analytic(std::list<std::string>&);
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
