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



// Makes sure only one console object exists.
bool Console::_lock {false};



Console::Console(int argc, char* argv[], Terminal& terminal)
/*
 * Prints welcome message to terminal interface, set lock to console.
 *
 * argc: the argc argument from main().
 * argv: the argv argument from main().
 * terminal: reference to terminal interface for input and output to user.
 */
   :_alive(false),
   _tm(terminal),
   _device(nullptr)
{
   InvalidUse::assert(!_lock,__FILE__,__LINE__);
   _lock = true;
   _tm.header("KINC:> ");
   _tm << "Welcome to KINC!\nVersion 0.0001 :)\n\n";
}



Console::~Console()
/*
 * Deletes OpenCL device if one is set, unset lock to console.
 */
{
   InvalidUse::assert(_lock,__FILE__,__LINE__);
   _lock = false;
   if (_device)
   {
      delete _device;
   }
}



void Console::run()
/*
 * Goes directly to terminal loop.
 */
{
   terminal_loop();
}



void Console::terminal_loop()
/*
 * This is the main function loop of the program. This loop will continue to
 * grab one input from the Terminal interface until the quit command has been
 * processed. Each user line is processed by calling parse.
 */
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



bool Console::parse(std::string& line)
/*
 * Takes one line of user input and seperates it into list of strings using
 * space or tab characters as delimiters between arguments. Once parsed into a
 * list passes onto function that decodes user input and processes the command.
 *
 * line: one line of user input.
 *
 * Returns true if the user command was processed successfully.
 */
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



bool Console::decode(std::list<std::string>& list)
/*
 * Decodes user command from what user typed.
 *
 * list: list of user arguments.
 *
 * Returns true if command was successful.
 */
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



bool Console::process(Command comm, std::list<std::string>& list)
/*
 * Processes decoded user command and routes to specific command.
 *
 * comm: user command to be processed.
 * list: list of user arguments for command.
 *
 * Returns true if command was successful.
 */
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



bool Console::gpu_decode(std::list<std::string>& list)
/*
 * Decodes the specific OpenCL command the user typed and passes it to a
 * processing function.
 *
 * list: list of arguments for OpenCL command.
 *
 * Returns true if the command was successful.
 */
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



bool Console::gpu_process(GpuCommand comm, std::list<std::string>& list)
/*
 * Processes decoded OpenCL command and routes to specific command given.
 *
 * comm: the OpenCL command to be processed.
 * list: list of arguments for this command.
 *
 * Returns true if the command was successful.
 */
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



void Console::gpu_list()
/*
 * Executes command to list all available OpenCL devices.
 */
{
   for (auto i = CLDevice::begin();i!=CLDevice::end();i++)
   {
      CLDevice& dev = *i;
      _tm << dev.info(CLInfo::ident) << " ";
      _tm << dev.info(CLInfo::name);
      if (_device&&dev==*_device)
      {
         _tm << " ***";
      }
      _tm << Terminal::endl;
   }
   _tm << Terminal::flush;
}



bool Console::gpu_info(std::list<std::string>& list)
/*
 * Executes command to print basic info of of given OpenCL device.
 *
 * list : list of arguments for this command.
 *
 * Returns true if the command was successful.
 */
{
   bool ret = false;
   if (list.size()>0)
   {
      int p,d;
      char sep;
      std::istringstream str(list.front());
      if ((str >> p >> sep >> d)&&sep==':'&&CLDevice::exist(p,d))
      {
         CLDevice dev(p,d);
         _tm << "===== " << dev.info(CLInfo::name) << " ("
             << dev.info(CLInfo::type) << ") =====\n";
         _tm << "Online: " << dev.info(CLInfo::online) << ".\n";
         _tm << "Unified Memory: " << dev.info(CLInfo::unified_mem) << ".\n";
         _tm << dev.info(CLInfo::addr_space) << " bit address space.\n";
         _tm << dev.info(CLInfo::clock) << "Mhz max clock frequency.\n";
         _tm << dev.info(CLInfo::compute_units) << " compute unit(s), "
             << dev.info(CLInfo::work_size) << " work-item(s) per unit.\n";
         _tm << dev.info(CLInfo::global_mem) << " global memory, "
             << dev.info(CLInfo::local_mem) << " local memory."
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



bool Console::gpu_set(std::list<std::string>& list)
/*
 * Executes command that sets OpenCL device for analytic computation.
 *
 * list : list of arguments for this command.
 *
 * Returns true if the command was successful.
 */
{
   bool ret = false;
   if (list.size()>0)
   {
      int p,d;
      char sep;
      std::istringstream str(list.front());
      if ((str >> p >> sep >> d)&&sep==':'&&CLDevice::exist(p,d))
      {
         if (_device)
         {
            delete _device;
         }
         _device = new CLDevice(p,d);
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



void Console::gpu_clear()
/*
 * Executes command that clears any OpenCL device set for computation.
 */
{
   if (_device)
   {
      delete _device;
      _device = nullptr;
   }
}
