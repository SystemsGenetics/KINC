#include "console.h"

Console g_console();



Console::Console():
   out(ConsoleStream::out),
   warn(ConsoleStream::warning),
   err(ConsoleStream::error)
{}
