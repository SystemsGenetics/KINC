#include "console.h"



int main(int argc, char* argv[])
{
   g_console.out << "Hello! " << (short)2 << " two " << 3.14f << " pie\n";
   g_console.err.print("hello! ");
   g_console.err.print((short)2);
   g_console.err.print(" two ");
   g_console.err.print(3.14f);
   g_console.err.print(" pie\n");
   g_console.run(argc,argv);
   return 0;
}
