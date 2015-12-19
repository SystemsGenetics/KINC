#include "console.h"


const char* g_weclomeMsg = "Welcome to KINC!\n"
      "Version 0.0001 :)\n";


int main(int argc, char* argv[])
{
   g_console.out << g_weclomeMsg << "\n";
   g_console.out.flush();
   g_console.run(argc,argv);
   return 0;
}
