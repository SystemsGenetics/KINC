/*
 * See GPL.txt for copyright information.
 *
 * Author: Joshua Burns
 *
 */
#include <vector>
#include <sstream>
#include "console.h"
#include "cldevice.h"
#include "data.h"
#include "analytic.h"
#include "plugins/plugins.h"



bool Console::_lock {false};



/// @brief Prepares console interface.
///
/// Initializes internal variables, locks console interface, and sets terminal
/// header.
///
/// @pre Console must not be locked.
Console::Console(int argc, char* argv[], Terminal& terminal, DataMap& datamap):
   _alive {false},
   _tm {terminal},
   _dataMap {datamap},
   _device {nullptr}
{
   assert<InvalidUse>(!_lock,__FILE__,__LINE__);
   _lock = true;
   _tm.header("KINC:> ");
}



/// Deletes OpenCL device if one is set and unsets lock to console.
Console::~Console()
{
   _lock = false;
   if (_device)
   {
      delete _device;
   }
}



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
   bool alive = true;
   while (alive)
   {
      std::string line;
      _tm >> line;
      _tm << Terminal::endl;
      try
      {
         parse(line);
      }
      catch (CommandError e)
      {
         _tm << e.cmd << ": " << e.msg << Terminal::endl;
      }
      catch (CommandQuit)
      {
         alive = false;
      }
   }
}



/// @brief Parses one line of user input into list of strings.
///
/// Takes one line of user input and seperates it into list of strings using
/// space or tab characters as delimiters between arguments. Once parsed into a
/// list passes onto function that decodes user input and processes the command.
///
/// @param line One line of user input.
void Console::parse(std::string& line)
{
   enum {_new,build} state = _new;
   std::list<std::string> list;
   char newBuf[2] = {" "};
   for (auto i:line)
   {
      switch (state)
      {
      case _new:
         if (i!=' '&&i!='\t')
         {
            newBuf[0] = i;
            list.push_back(newBuf);
            state = build;
         }
         break;
      case build:
         if (i!=' '&&i!='\t')
         {
            list.back() += i;
         }
         else
         {
            state = _new;
         }
         break;
      }
   }
   decode(list);
}



/// @brief Decodes string command into enumerated type.
///
/// Takes the first argument in user command and decodes into specific
/// enumerated command type. Possibly modified command argument list and
/// enumerated command type is then passed to subfunction.
///
/// @param list List of user arguments from command.
void Console::decode(std::list<std::string>& list)
{
   Command comm = Command::error;
   if (list.size()>0)
   {
      if (list.front()=="gpu")
      {
         comm = Command::gpu;
         list.pop_front();
      }
      else if (list.front()=="open")
      {
         comm = Command::open;
         list.pop_front();
      }
      else if (list.front()=="load")
      {
         comm = Command::load;
         list.pop_front();
      }
      else if (list.front()=="dump")
      {
         comm = Command::dump;
         list.pop_front();
      }
      else if (list.front()=="query")
      {
         comm = Command::query;
         list.pop_front();
      }
      else if (list.front()=="close")
      {
         comm = Command::close;
         list.pop_front();
      }
      else if (list.front()=="list")
      {
         comm = Command::close;
         list.pop_front();
      }
      else if (list.front()=="quit")
      {
         comm = Command::quit;
         list.pop_front();
      }
      else
      {
         comm = Command::analytic;
      }
   }
   process(comm,list);
}



/// Processes decoded user command and routes to specific command.
///
/// @param comm Enumerated command to be processed.
/// @param list List of user arguments for command.
void Console::process(Command comm, std::list<std::string>& list)
{
   switch (comm)
   {
   case Command::gpu:
      gpu_decode(list);
      break;
   case Command::open:
   case Command::load:
   case Command::dump:
   case Command::query:
   case Command::close:
      data_main(comm,list);
      break;
   case Command::list:
      break;
   case Command::analytic:
      break;
   case Command::quit:
      throw CommandQuit();
   case Command::error:
      return;
   }
}



/// @brief decodes string OpenCL subcommand into enumerated type.
///
/// Takes the first argument of the OpenCL subcommand and decodes to proper
/// enumerated command. Possibly modified subcommand argument list and
/// enumerated command type is then passed to subfunction.
///
/// @param list List of arguments for OpenCL command.
void Console::gpu_decode(std::list<std::string>& list)
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
         throw CommandError("gpu","Unknown command.");
      }
   }
   else
   {
      throw CommandError("gpu","No command given.");
   }
   gpu_process(comm,list);
}



/// Processes decoded OpenCL command and routes to specific command given.
///
/// @param comm The OpenCL command to be processed.
/// @param list List of arguments for this command.
void Console::gpu_process(GpuCommand comm, std::list<std::string>& list)
{
   switch (comm)
   {
   case GpuCommand::list:
      gpu_list();
      break;
   case GpuCommand::info:
      gpu_info(list);
      break;
   case GpuCommand::set:
      gpu_set(list);
      break;
   case GpuCommand::clear:
      gpu_clear();
      break;
   }
}



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



/// Executes command to print basic info of of given OpenCL device.
///
/// @param list List of arguments for this command.
void Console::gpu_info(std::list<std::string>& list)
{
   if (list.size()==0)
   {
      throw CommandError("gpu info","command requires 1 argument.");
   }
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
   }
   else
   {
      std::ostringstream buffer;
      buffer << "cannot find OpenCL device \"" << list.front() << "\".";
      throw CommandError("gpu info",buffer.str());
   }
}



/// Executes command that sets OpenCL device for analytic computation.
///
/// @param list List of arguments for this command.
void Console::gpu_set(std::list<std::string>& list)
{
   if (list.size()==0)
   {
      throw CommandError("gpu info","command requires 1 argument.");
   }
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
      _tm << "OpenCL device set to \"" << list.front() << "\"."
          << Terminal::endl;
   }
   else
   {
      std::ostringstream buffer;
      buffer << "cannot find OpenCL device \"" << list.front() << "\".";
      throw CommandError("gpu info",buffer.str());
   }
}



/// Executes command that clears any OpenCL device set for computation.
void Console::gpu_clear()
{
   if (_device)
   {
      delete _device;
      _device = nullptr;
      _tm << "OpenCL device cleared." << Terminal::endl;
   }
   else
   {
      _tm << "no OpenCL devie set." << Terminal::endl;
   }
}



void Console::data_main(Command comm, std::list<std::string>& list)
{
   try
   {
      switch (comm)
      {
      case Command::open:
         data_open(list);
         break;
      case Command::load:
         data_load(list);
         break;
      case Command::dump:
         data_dump(list);
         break;
      case Command::query:
         data_query(list);
         break;
      case Command::close:
         data_close(list);
         break;
      }
   }
   catch (DataException e)
   {
      if (e.level()==DataException::Level::fatal)
      {
         _dataMap.del(e.who());
         if (e.who())
         {
            delete e.who();
         }
      }
      std::ostringstream buffer;
      buffer << "Data Exception Caught!" << std::endl;
      buffer << "Level: ";
      switch (e.level())
      {
      case DataException::Level::general:
         buffer << "General" << std::endl;
         break;
      case DataException::Level::caution:
         buffer << "Caution" << std::endl;
         break;
      case DataException::Level::warning:
         buffer << "Warning" << std::endl;
         break;
      case DataException::Level::severe:
         buffer << "Severe" << std::endl;
         break;
      case DataException::Level::fatal:
         buffer << "Fatal" << std::endl;
         break;
      }
      buffer << "Location: " << e.file() << ":" << e.line() << std::endl;
      buffer << "What: " << e.what();
      throw CommandError("data",buffer.str());
   }
}



void Console::data_open(std::list<std::string>& list)
{
   if (list.size()<2)
   {
      throw CommandError("open","command requires 2 arguments.");
   }
   std::string type = list.front();
   list.pop_front();
   std::string name = list.front();
   if (_dataMap.find(name)!=_dataMap.end())
   {
      std::ostringstream buffer;
      buffer << "data object with name \"" << name << "\" already exists.";
      throw CommandError("open",buffer.str());
   }
   DataPlugin* nd {nullptr};
   try
   {
      nd = KINCPlugins::new_data(type,name);
   }
   catch (DataException)
   {
      throw;
   }
   catch (...)
   {
      throw DataException(__FILE__,__LINE__,nd,"UNKNOWN",
                          DataException::Level::fatal);
   }
   if (!nd)
   {
      std::ostringstream buffer;
      buffer << "cannot find data type \"" << type << "\".";
      throw CommandError("open",buffer.str());
   }
   _dataMap.add(name,nd);
   _tm << "Added " << name << "(" << type << ") data object." << Terminal::endl;
}



void Console::data_load(std::list<std::string>& list)
{
   if (list.size()<2)
   {
      throw CommandError("load","command requires 2 arguments.");
   }
   std::string textFile = list.front();
   list.pop_front();
   DataPlugin* data = find_data(list.front());
   list.pop_front();
   try
   {
      parse_data_options(data,list);
      data->load(textFile,_tm);
   }
   catch (DataException)
   {
      throw;
   }
   catch (...)
   {
      throw DataException(__FILE__,__LINE__,data,"UNKNOWN",
                          DataException::Level::fatal);
   }
}



void Console::data_dump(std::list<std::string>& list)
{
   if (list.size()<2)
   {
      throw CommandError("dump","command requires 2 arguments.");
   }
   std::string textFile = list.front();
   list.pop_front();
   DataPlugin* data = find_data(list.front());
   list.pop_front();
   try
   {
      parse_data_options(data,list);
      data->dump(textFile,_tm);
   }
   catch (DataException)
   {
      throw;
   }
   catch (...)
   {
      throw DataException(__FILE__,__LINE__,data,"UNKNOWN",
                          DataException::Level::fatal);
   }
}



void Console::data_query(std::list<std::string>& list)
{
   if (list.size()<1)
   {
      throw CommandError("query","command requires 1 argument.");
   }
   DataPlugin* data = find_data(list.front());
   list.pop_front();
   try
   {
      parse_data_options(data,list);
      data->query(_tm);
   }
   catch (DataException)
   {
      throw;
   }
   catch (...)
   {
      throw DataException(__FILE__,__LINE__,data,"UNKNOWN",
                          DataException::Level::fatal);
   }
}



void Console::data_close(std::list<std::string>& list)
{
}



void Console::data_list(std::list<std::string>& list)
{
}



void Console::analytic(std::list<std::string>& list)
{
}



DataPlugin* Console::find_data(const std::string& name)
{
   auto i = _dataMap.find(name);
   if (i==_dataMap.end())
   {
      std::ostringstream buffer;
      buffer << "data object \"" << name << "\" not found.";
      throw CommandError("load",buffer.str());
   }
   return i->second;
}



void Console::parse_data_options(DataPlugin* data, std::list<std::string>& list)
{
   while (list.size()>0)
   {
      std::string key;
      std::string val;
      enum {front,back} state = front;
      for (auto i:list.front())
      {
         switch (state)
         {
         case front:
            if (i=='=')
            {
               state = back;
            }
            else
            {
               key += i;
            }
            break;
         case back:
            val += i;
            break;
         }
      }
      data->option(key,val);
   }
}
