/*
 * See GPL.txt for copyright information.
 *
 * Author: Joshua Burns
 *
 */
#include <vector>
#include <sstream>
#include <forward_list>
#include <memory>
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
      string line;
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
      catch (DataException e)
      {
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
      catch (AnalyticException e)
      {
         ;
      }
      catch (Exception e)
      {
         ;
      }
      catch (std::exception stde)
      {
         ;
      }
      catch (...)
      {
         ;
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
void Console::parse(string& line)
{
   enum {_new,build} state = _new;
   slist list;
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
void Console::decode(slist& list)
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
         comm = Command::list;
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
void Console::process(Command comm, slist& list)
{
   switch (comm)
   {
   case Command::gpu:
      gpu_decode(list);
      break;
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
   case Command::list:
      data_list();
      break;
   case Command::analytic:
      analytic(list);
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
void Console::gpu_decode(slist& list)
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
void Console::gpu_process(GpuCommand comm, slist& list)
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
void Console::gpu_info(slist& list)
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
void Console::gpu_set(slist& list)
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



void Console::data_open(slist& list)
{
   if (list.size()<2)
   {
      throw CommandError("open","command requires 2 arguments.");
   }
   string type = list.front();
   list.pop_front();
   if (_dataMap.exist(list.front()))
   {
      std::ostringstream buffer;
      buffer << "data object with name \"" << list.front()
             << "\" already exists.";
      throw CommandError("open",buffer.str());
   }
   std::unique_ptr<DataPlugin> nd(KINCPlugins::new_data(type,list.front()));
   if (!nd)
   {
      std::ostringstream buffer;
      buffer << "cannot find data type \"" << type << "\".";
      throw CommandError("open",buffer.str());
   }
   _dataMap.add(list.front(),nd.release());
   _tm << "Added " << list.front() << "(" << type << ") data object."
       << Terminal::endl;
}



void Console::data_load(slist& list)
{
   if (list.size()<2)
   {
      throw CommandError("load","command requires 2 arguments.");
   }
   string textFile = list.front();
   list.pop_front();
   DataPlugin* data = find_data(list.front());
   list.pop_front();
   try
   {
      parse_data_options(data,list);
      data->load(textFile,_tm);
   }
   catch (DataException e)
   {
      if (e.level()==DataException::Level::fatal)
      {
         _dataMap.del(data);
      }
      throw;
   }
   catch (...)
   {
      _dataMap.del(data);
      throw;
   }
}



void Console::data_dump(slist& list)
{
   if (list.size()<2)
   {
      throw CommandError("dump","command requires 2 arguments.");
   }
   string textFile = list.front();
   list.pop_front();
   DataPlugin* data = find_data(list.front());
   list.pop_front();
   try
   {
      parse_data_options(data,list);
      data->dump(textFile,_tm);
   }
   catch (DataException e)
   {
      if (e.level()==DataException::Level::fatal)
      {
         _dataMap.del(data);
      }
      throw;
   }
   catch (...)
   {
      _dataMap.del(data);
      throw;
   }
}



void Console::data_query(slist& list)
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
   catch (DataException e)
   {
      if (e.level()==DataException::Level::fatal)
      {
         _dataMap.del(data);
      }
      throw;
   }
   catch (...)
   {
      _dataMap.del(data);
      throw;
   }
}



void Console::data_close(slist& list)
{
   if (list.size()<1)
   {
      throw CommandError("close","command requires 1 argument.");
   }
   DataPlugin* data = _dataMap.find(list.front());
   if (!data)
   {
      std::ostringstream buffer;
      buffer << "data object \"" << list.front() << "\" not found.";
      throw CommandError("close",buffer.str());
   }
   else
   {
      _dataMap.del(data);
      _tm << "data object closed." << Terminal::endl;
   }
}



void Console::data_list()
{
   for (auto i=_dataMap.begin();i!=_dataMap.end();++i)
   {
      DataPlugin* data = i->second;
      _tm << i->first << " [" << data->type() << "]" << Terminal::endl;
   }
}



void Console::analytic(slist& list)
{
   if (list.size()<3)
   {
      throw CommandError("analytic","command requires 2 argument.");
   }
   string type = list.front();
   list.pop_front();
   aptr a(KINCPlugins::new_analytic(type));
   if (!a)
   {
      std::ostringstream buffer;
      buffer << "cannot find analytic type \"" << type << "\".";
      throw CommandError("open",buffer.str());
   }
   dlist nd;
   try
   {
      parse_analytic_inputs(a,list.front());
      list.pop_front();
      parse_analytic_outputs(a,list.front(),nd);
      parse_analytic_options(a,list);
      a->execute(_tm,&(_device->device()));
      for (auto i:nd)
      {
         _dataMap.add(i.first,i.second);
      }
   }
   catch (...)
   {
      for (auto i:nd)
      {
         if (i.second) delete i.second;
      }
   }
}



DataPlugin* Console::find_data(const string& name)
{
   DataPlugin* ret {_dataMap.find(name)};
   if (!ret)
   {
      std::ostringstream buffer;
      buffer << "data object \"" << name << "\" not found.";
      throw CommandError("find",buffer.str());
   }
   return ret;
}



void Console::parse_data_options(DataPlugin* data, slist& list)
{
   while (list.size()>0)
   {
      string key;
      string val;
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
      list.pop_front();
   }
}



void Console::parse_analytic_inputs(aptr& a, const string& arg)
{
   string name;
   auto i = arg.begin();
   while (true)
   {
      if (*i==','||i==arg.end())
      {
         if (name.size()>0)
         {
            a->input(find_data(name));
         }
         name.clear();
      }
      else
      {
         name += *i;
      }
      if (i==arg.end())
      {
         break;
      }
      ++i;
   }
}



void Console::parse_analytic_outputs(aptr& a, const string& arg, dlist& nd)
{
   string ndata;
   auto i = arg.begin();
   while (true)
   {
      if (*i==','||i==arg.end())
      {
         if (ndata.size()>0)
         {
            std::string name;
            DataPlugin* n {nullptr};
            try
            {
               n = parse_analytic_ndata(ndata,name);
               a->input(n);
            }
            catch(...)
            {
               if (n) delete n;
               throw;
            }
            nd.push_back({name,n});
         }
         ndata.clear();
      }
      else
      {
         ndata += *i;
      }
      if (i==arg.end())
      {
         break;
      }
      ++i;
   }
}



DataPlugin* Console::parse_analytic_ndata(const string& ndata, string& name)
{
   string type;
   enum {front,back} state = front;
   for (auto i:ndata)
   {
      switch (state)
      {
      case front:
         if (i==':')
         {
            state = back;
         }
         else
         {
            type += i;
         }
         break;
      case back:
         name += i;
         break;
      }
   }
   return KINCPlugins::new_data(type,name);
}



void Console::parse_analytic_options(aptr& a, slist& list)
{
   while (list.size()>0)
   {
      string key;
      string val;
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
      a->option(key,val);
      list.pop_front();
   }
}
