#ifndef LINUXTERM_H
#define LINUXTERM_H
#include <string>
#include <list>
#include "terminal.h"



/// @brief Implements Terminal in Linux OS environment.
///
/// Implements terminal interface for user in linux terminal environment. Uses
/// system calls to get the width of the terminal and VT100 and ANSI escape
/// codes to control the terminal in raw mode to create an interactive terminal
/// environment. See https://en.wikipedia.org/wiki/ANSI_escape_code for more
/// information about control characters used.
///
/// @pre Only one instance of this class can exist at any one time.
/// @pre stty_raw() must be called before instantiating this class.
/// @pre stty_cooked() must be called after destroying this class.
/// @pre The terminal this program is running in cannot change its width during
/// the lifetime of an instance.
class LinuxTerm : public Terminal
{
public:
   // ****************************** Deleted Methods ************************
   LinuxTerm(const LinuxTerm&) = delete;
   LinuxTerm& operator=(const LinuxTerm&) = delete;
   LinuxTerm& operator=(LinuxTerm&&) = delete;
   // ****************************** Static Functions ***********************
   static void stty_raw();
   static void stty_cooked();
   // ****************************** Constructors/Destructors ***************
   LinuxTerm();
   ~LinuxTerm();
   // ****************************** Functions ******************************
   void header(const std::string&);
   void set_ops(Ops);
   void precision(int);
   // ****************************** Operators ******************************
   LinuxTerm& operator<<(short);
   LinuxTerm& operator<<(unsigned short);
   LinuxTerm& operator<<(int);
   LinuxTerm& operator<<(unsigned int);
   LinuxTerm& operator<<(long);
   LinuxTerm& operator<<(unsigned long);
   LinuxTerm& operator<<(float);
   LinuxTerm& operator<<(double);
   LinuxTerm& operator<<(const char*);
   LinuxTerm& operator<<(const std::string&);
   void operator>>(std::string&);
private:
   // ****************************** Functions ******************************
   void reset_cursor(int);
   void reprint(bool);
   // ****************************** Constants ******************************
   constexpr static char _backspaceCh {'\x7f'};
   constexpr static char _escapeCh {'\x1b'};
   constexpr static char _arrowRightCh {'C'};
   constexpr static char _arrowLeftCh {'D'};
   const char* _cursorUpStr {"\x1b[A"};
   const char* _boldTextStr {"\x1b[1m"};
   const char* _normTextStr {"\x1b[0m"};
   // ****************************** Static Variables ***********************
   static bool _lock;
   static bool _cooked;
   // ****************************** Variables ******************************
   int _cols;
   int _chCount;
   std::string _header;
   std::list<char> _line;
   std::list<char>::iterator _i;
};



#endif
