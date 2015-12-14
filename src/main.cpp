#include "console.h"



int main(int argc, char* argv[])
{
   g_console.out << "Hello! " << (short)2 << " two " << 3.14f << " pie\n";
   g_console.err.print("hello");
   g_console.err.print("foo\n");
   return 0;
}
