#define __CL_ENABLE_EXCEPTIONS
#include <iostream>
#include <CL/cl.hpp>
#include <vector>
#include <stdio.h>
#include "console.h"
#include "terminal.h"
#include "exception.h"



Console g_console;



bool Console::gpu(std::list<std::string>& list)
{
   bool ret = false;
   if (list.size()>0)
   {
      if (list.front()=="list")
      {
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
                  if (_gpu.platform==i&&_gpu.device==j)
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
      }
      else if (list.front()=="set")
      {
         if (list.size()>1)
         {
            list.pop_front();
            int platform;
            int device;
            if (sscanf(list.front().c_str(),"%d:%d",&platform,&device)!=2||platform<0||device<0)
            {
               std::cout << "Error: argument given for OpenCL device ADDR is invalid." << std::endl;
            }
            else
            {
               try
               {
                  std::vector<cl::Platform> platforms;
                  cl::Platform::get(&platforms);
                  if (platform<platforms.size())
                  {
                     std::vector<cl::Device> devices;
                     platforms[platform].getDevices(CL_DEVICE_TYPE_ALL,&devices);
                     if (device<devices.size())
                     {
                        _gpu.platform = platform;
                        _gpu.device = device;
                        ret = true;
                     }
                  }
               }
               catch (cl::Error e)
               {
                  throw OpenCLError(__FILE__,__LINE__,e);
               }

               if (!ret)
               {
                  std::cout << "Error: cannot find OpenCL device " << platform << ":" << device
                            << "." << std::endl;
               }
            }
         }
         else
         {
            std::cout << "Error: the gpu set command requires one argument." << std::endl;
         }
      }
      else if (list.front()=="clear")
      {
         _gpu.platform = -1;
         _gpu.device = -1;
         ret = true;
      }
      else
      {
         std::cout << "Error: " << list.front() << " GPU subcommand not found." << std::endl;
      }
   }
   else
   {
      std::cout << "Error: GPU subcommand required." << std::endl;
   }
   return ret;
}



bool Console::process(std::string& line)
{
   bool ret = false;
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
   if (list.size()>0)
   {
      if (list.front()=="gpu")
      {
         list.pop_front();
         ret = gpu(list);
      }
      else if (list.front()=="quit")
      {
         ret = true;
      }
      else
      {
         std::cout << "Error: " << list.front() << " command/analytic not found." << std::endl;
      }
   }
   return ret;
}



void Console::command()
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
      process(termCommand);
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
   _gpu.platform = -1;
   _gpu.device = -1;
}



void Console::run(int argc, char* argv[])
{
   command();
}
