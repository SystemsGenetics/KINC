/*
 * See GPL.txt for copyright information.
 *
 * Author: Joshua Burns
 *
 */
#include <vector>
#include <sstream>
#include "console.h"
#include "exception.h"
#include "cldevice.h"



bool Console::_lock {false};



/// @brief Prepares console interface.
///
/// Initializes internal variables, locks console interface, and sets terminal
/// header.
///
/// @pre Console must not be locked.
Console::Console(int argc, char* argv[], Terminal& terminal):
   _alive {false},
   _tm {terminal},
   _device {nullptr}
{
   InvalidUse::assert(!_lock,__FILE__,__LINE__);
   _lock = true;
   _tm.header("KINC:> ");
}



/// @brief Cleans up after console interface.
///
/// Deletes OpenCL device if one is set and unsets lock to console.
Console::~Console()
{
   _lock = false;
   if (_device)
   {
      delete _device;
   }
}



/// @brief Main program starter.
///
/// Prints welcome message and goes directly to terminal loop.
void Console::run()
{
   _tm << "Welcome to KINC!\nVersion 0.0001 :)\n\n";
   terminal_loop();
}



/// @brief Main program terminal loop.
///
/// This is the main function loop of the program. This loop will continue to
/// grab one input from the Terminal interface until the quit command has been
/// processed. Each user line is processed by calling parse. The _alive state
/// of the program is set to false in a subfunction if quit is given.
void Console::terminal_loop()
{
   _alive = true;
   while (_alive)
   {
      std::string line;
      _tm >> line;
      _tm << Terminal::endl;
      parse(line);
      _tm << Terminal::endl;
   }
}



/// @brief Parses one line of user input into list of strings.
///
/// Takes one line of user input and seperates it into list of strings using
/// space or tab characters as delimiters between arguments. Once parsed into a
/// list passes onto function that decodes user input and processes the command.
///
/// @param line One line of user input.
///
/// @return True if the user command was processed successfully.
bool Console::parse(std::string& line)
{
   enum {_new,build} state = _new;
   std::list<std::string> list;
   char newBuf[2] = {" "};
   for (std::string::iterator i=line.begin();i!=line.end();i++)
   {
      switch (state)
      {
      case _new:
         if (*i!=' '&&*i!='\t')
         {
            newBuf[0] = *i;
            list.push_back(newBuf);
            state = build;
         }
         break;
      case build:
         if (*i!=' '&&*i!='\t')
         {
            list.back() += *i;
         }
         else
         {
            state = _new;
         }
         break;
      }
   }
   return decode(list);
}



/// @brief Decodes string command into enumerated type.
///
/// Takes the first argument in user command and decodes into specific
/// enumerated command type. Possibly modified command argument list and
/// enumerated command type is then passed to subfunction.
///
/// @param list List of user arguments from command.
///
/// @return True if command was successful.
bool Console::decode(std::list<std::string>& list)
{
   Command comm;
   if (list.size()>0)
   {
      if (list.front()=="gpu")
      {
         comm = Command::gpu;
         list.pop_front();
      }
      else if (list.front()=="quit")
      {
         comm = Command::quit;
         list.pop_front();
      }
      else
      {
         comm = Command::error;
         _tm << "Error: " << list.front() << " command/analytic not found."
             << Terminal::endl;
      }
   }
   else
   {
      comm = Command::error;
      _tm << "Error: empty command given." << Terminal::endl;
   }
   return process(comm,list);
}



/// @brief Processes enumerated command.
///
/// Processes decoded user command and routes to specific command.
///
/// @param comm Enumerated command to be processed.
/// @param list List of user arguments for command.
///
/// @return True if command was successful.
bool Console::process(Command comm, std::list<std::string>& list)
{
   bool ret;
   switch (comm)
   {
   case Command::gpu:
      ret = gpu_decode(list);
      break;
   case Command::quit:
      _alive = false;
      ret = true;
      break;
   case Command::error:
      ret = false;
   }
   return ret;
}



/// @brief decodes string OpenCL subcommand into enumerated type.
///
/// Takes the first argument of the OpenCL subcommand and decodes to proper
/// enumerated command. Possibly modified subcommand argument list and
/// enumerated command type is then passed to subfunction.
///
/// @param list List of arguments for OpenCL command.
///
/// @return True if the command was successful.
bool Console::gpu_decode(std::list<std::string>& list)
{
   GpuCommand comm;
   if (list.size()>0)
   {
      if (list.front()=="list")
      {
         comm = GpuCommand::list;
         list.pop_front();
      }
      else if (list.front()=="info")
      {
         comm = GpuCommand::info;
         list.pop_front();
      }
      else if (list.front()=="set")
      {
         comm = GpuCommand::set;
         list.pop_front();
      }
      else if (list.front()=="clear")
      {
         comm = GpuCommand::clear;
         list.pop_front();
      }
      else
      {
         comm = GpuCommand::error;
         _tm << "Error: " << list.front() << " GPU subcommand not found."
             << Terminal::endl;
      }
   }
   else
   {
      comm = GpuCommand::error;
      _tm << "Error: GPU subcommand required." << Terminal::endl;
   }
   return gpu_process(comm,list);
}



/// @brief Processes enumerated OpenCL subcommand.
///
/// Processes decoded OpenCL command and routes to specific command given.
///
/// @param comm The OpenCL command to be processed.
/// @param list List of arguments for this command.
///
/// @return True if the command was successful.
bool Console::gpu_process(GpuCommand comm, std::list<std::string>& list)
{
   bool ret;
   switch (comm)
   {
   case GpuCommand::list:
      gpu_list();
      ret = true;
      break;
   case GpuCommand::info:
      ret = gpu_info(list);
      break;
   case GpuCommand::set:
      ret = gpu_set(list);
      break;
   case GpuCommand::clear:
      gpu_clear();
      ret = true;
      break;
   case GpuCommand::error:
      ret = false;
      break;
   }
   return ret;
}



/// @brief List all OpenCL devices.
///
/// Executes command to list all available OpenCL devices.
void Console::gpu_list()
{
   for (auto d:_devList)
   {
      _tm << d.info(CLDevice::ident) << " ";
      _tm << d.info(CLDevice::name);
      if (_device&&d==*_device)
      {
         _tm << " ***";
      }
      _tm << Terminal::endl;
   }
   _tm << Terminal::flush;
}



/// @brief Prints info about specific OpenCL device.
///
/// Executes command to print basic info of of given OpenCL device.
///
/// @param list List of arguments for this command.
///
/// @return True if the command was successful.
bool Console::gpu_info(std::list<std::string>& list)
{
   bool ret = false;
   if (list.size()>0)
   {
      int p,d;
      char sep;
      std::istringstream str(list.front());
      if ((str >> p >> sep >> d)&&sep==':'&&_devList.exist(p,d))
      {
         CLDevice& dev {_devList.at(p,d)};
         _tm << "===== " << dev.info(CLDevice::name) << " ("
             << dev.info(CLDevice::type) << ") =====\n";
         _tm << "Online: " << dev.info(CLDevice::online) << ".\n";
         _tm << "Unified Memory: " << dev.info(CLDevice::unified_mem) << ".\n";
         _tm << dev.info(CLDevice::addr_space) << " bit address space.\n";
         _tm << dev.info(CLDevice::clock) << "Mhz max clock frequency.\n";
         _tm << dev.info(CLDevice::compute_units) << " compute unit(s), "
             << dev.info(CLDevice::work_size) << " work-item(s) per unit.\n";
         _tm << dev.info(CLDevice::global_mem) << " global memory, "
             << dev.info(CLDevice::local_mem) << " local memory."
             << Terminal::endl;
         ret = true;
      }
      else
      {
         _tm << "Error: cannot find OpenCL device " << list.front()
                    << Terminal::endl;
      }
   }
   else
   {
      _tm << "Error: the gpu info command requires one argument."
                 << Terminal::endl;
   }
   return ret;
}



/// @brief Sets OpenCL device for use in console.
///
/// Executes command that sets OpenCL device for analytic computation.
///
/// @param list List of arguments for this command.
///
/// @return True if the command was successful.
bool Console::gpu_set(std::list<std::string>& list)
{
   bool ret = false;
   if (list.size()>0)
   {
      int p,d;
      char sep;
      std::istringstream str(list.front());
      if ((str >> p >> sep >> d)&&sep==':'&&_devList.exist(p,d))
      {
         if (_device)
         {
            delete _device;
         }
         _device = new CLDevice {_devList.at(p,d)};
         ret = true;
      }
      else
      {
         _tm << "Error: cannot find OpenCL device " << list.front()
                    << Terminal::endl;
      }
   }
   else
   {
      _tm << "Error: the gpu set command requires one argument."
                 << Terminal::endl;
   }
   return ret;
}



/// @brief Clear any previously set OpenCL device.
///
/// Executes command that clears any OpenCL device set for computation.
void Console::gpu_clear()
{
   if (_device)
   {
      delete _device;
      _device = nullptr;
   }
}
