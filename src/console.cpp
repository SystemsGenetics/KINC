#include <iostream>
#include <vector>
#include <stdio.h>
#include "console.h"
#include "terminal.h"
#include "exception.h"



Console g_console;



cl::Device* Console::gpu_get_device(std::string& ident, int& i, int& j)
{
   cl::Device* ret = NULL;
   if (sscanf(ident.c_str(),"%d:%d",&i,&j)==2&&i>=0&&j>=0)
   {
      try
      {
         std::vector<cl::Platform> platforms;
         cl::Platform::get(&platforms);
         if (i<platforms.size())
         {
            std::vector<cl::Device> devices;
            platforms[i].getDevices(CL_DEVICE_TYPE_ALL,&devices);
            if (j<devices.size())
            {
               ret = new cl::Device(devices[j]);
            }
         }
      }
      catch (cl::Error e)
      {
         throw OpenCLError(__FILE__,__LINE__,e);
      }
   }
   return ret;
}



bool Console::gpu_process(command_e comm, std::list<std::string>& list)
{
   bool ret;
   switch (comm)
   {
   case glist:
      std::cout << "ADDR    NAME\n";
      try
      {
         std::vector<cl::Platform> platforms;
         cl::Platform::get(&platforms);
         for (int i=0;i<platforms.size();i++)
         {
            std::vector<cl::Device> devices;
            platforms[i].getDevices(CL_DEVICE_TYPE_ALL,&devices);
            for (int j=0;j<devices.size();j++)
            {
               std::cout << i << ":" << j;
               if (_gpu.i==i&&_gpu.j==j)
               {
                  std::cout << " *** ";
               }
               else
               {
                  std::cout << "     ";
               }
               std::cout << devices[j].getInfo<CL_DEVICE_NAME>() << "\n";
            }
         }
      }
      catch (cl::Error e)
      {
         throw OpenCLError(__FILE__,__LINE__,e);
      }
      std::cout << std::flush;
      ret = true;
      break;
   case ginfo:
      ret = false;
      if (list.size()>0)
      {
         int i;
         int j;
         cl::Device* ndev;
         if ((ndev = gpu_get_device(list.front(),i,j))!=NULL)
         {
            std::cout << ndev->getInfo<CL_DEVICE_NAME>() << std::endl;
            delete ndev;
            ret = true;
         }
         else
         {
            std::cout << "Error: cannot find OpenCL device " << list.front()
                      << std::endl;
         }
      }
      else
      {
         std::cout << "Error: the gpu info command requires one argument."
                   << std::endl;
      }
      break;
   case gset:
      ret = false;
      if (list.size()>0)
      {
         int i;
         int j;
         cl::Device* ndev;
         if ((ndev = gpu_get_device(list.front(),i,j))!=NULL)
         {
            _gpu.i = i;
            _gpu.j = j;
            _gpu.device = ndev;
            ret = true;
         }
         else
         {
            std::cout << "Error: cannot find OpenCL device " << list.front()
                      << std::endl;
         }
      }
      else
      {
         std::cout << "Error: the gpu set command requires one argument."
                   << std::endl;
      }
      break;
   case gclear:
      _gpu.i = -1;
      _gpu.j = -1;
      if (_gpu.device!=NULL)
      {
         delete _gpu.device;
      }
      _gpu.device = NULL;
      ret = true;
      break;
   default:
      ret = false;
      break;
   }
   return ret;
}



bool Console::gpu_decode(std::list<std::string>& list)
{
   command_e comm;
   if (list.size()>0)
   {
      if (list.front()=="list")
      {
         comm = glist;
         list.pop_front();
      }
      else if (list.front()=="info")
      {
         comm = ginfo;
         list.pop_front();
      }
      else if (list.front()=="set")
      {
         comm = gset;
         list.pop_front();
      }
      else if (list.front()=="clear")
      {
         comm = gclear;
         list.pop_front();
      }
      else
      {
         comm = error;
         std::cout << "Error: " << list.front()
                   << " GPU subcommand not found." << std::endl;
      }
   }
   else
   {
      comm = error;
      std::cout << "Error: GPU subcommand required." << std::endl;
   }
   return gpu_process(comm,list);
}



bool Console::process(command_e comm, std::list<std::string>& list)
{
   bool ret;
   switch (comm)
   {
   case gpu:
      ret = gpu_decode(list);
      break;
   case quit:
      ret = true;
      break;
   default:
      ret = false;
   }
   return ret;
}



bool Console::decode(std::list<std::string>& list)
{
   command_e comm;
   if (list.size()>0)
   {
      if (list.front()=="gpu")
      {
         comm = gpu;
         list.pop_front();
      }
      else if (list.front()=="quit")
      {
         comm = quit;
      }
      else
      {
         comm = error;
         std::cout << "Error: " << list.front()
                   << " command/analytic not found." << std::endl;
      }
   }
   else
   {
      comm = error;
      std::cout << "Error: empty command given." << std::endl;
   }
   return process(comm,list);
}



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



void Console::get_terminal_line_loop()
{
   Terminal cline;
   char input;
   bool alive = true;
   while (alive)
   {
      cline.clear("KINC:> ");
      while (!cline.read());
      std::string termCommand = cline.get();
      std::cout << "\n";
      parse(termCommand);
      if (termCommand=="quit")
      {
         alive = false;
      }
   }
}



Console::Console():
   out(ConsoleStream::out),
   warn(ConsoleStream::warning),
   err(ConsoleStream::error)
{
   _gpu.i = -1;
   _gpu.j = -1;
   _gpu.device = NULL;
}



Console::~Console()
{
   if (_gpu.device!=NULL)
   {
      delete _gpu.device;
   }
}



void Console::run(int argc, char* argv[])
{
   get_terminal_line_loop();
}
