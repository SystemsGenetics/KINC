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
#include "exception.h"



class Data;
class Analytic;
class CLDevice;



/*
 * NOTE: for new_data and new_analytic see plugin.cpp for source code and
 * comments.
 */
namespace KINCPlugins
{
   extern Data* new_kinc_data(std::string&);
   extern Analytic* new_kinc_analytic(std::string&);
}



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
   Console(int argc, char* argv[], Terminal& tm);
   ~Console();
   // *
   // * FUNCTIONS
   // *
   void run();
   //bool add(Data*,std::string&);
   //bool del(std::string&);
   //Data* find(std::string&);
private:
   // *
   // * ENUMERATIONS
   // *
   ///
   /// Main commands.
   enum class Command
   {
      gpu, ///< This is an OpenCL command.
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
      clear, ///< The clear subcommand.
      error ///< Error at parsing OpenCL subcommand.
   };
   // *
   // * FUNCTIONS
   // *
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
   // *
   // * STATIC VARIABLES
   // *
   static bool _lock;
   // *
   // * VARIABLES
   // *
   /// Reference to program's main Terminal interface.
   Terminal& _tm;
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
