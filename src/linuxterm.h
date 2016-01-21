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
   // ****************************** Basic Methods **************************
   LinuxTerm(const LinuxTerm&) = delete;
   LinuxTerm(LinuxTerm&&) = delete;
   LinuxTerm& operator=(const LinuxTerm&) = delete;
   LinuxTerm& operator=(LinuxTerm&&) = delete;
   LinuxTerm();
   ~LinuxTerm() override final;
   // ****************************** Static Functions ***********************
   static void stty_raw();
   static void stty_cooked();
   // ****************************** Functions ******************************
   void header(const std::string&) override final;
   void precision(int) override final;
   // ****************************** Operators ******************************
   LinuxTerm& operator<<(short) override final;
   LinuxTerm& operator<<(unsigned short) override final;
   LinuxTerm& operator<<(int) override final;
   LinuxTerm& operator<<(unsigned int) override final;
   LinuxTerm& operator<<(long) override final;
   LinuxTerm& operator<<(unsigned long) override final;
   LinuxTerm& operator<<(float) override final;
   LinuxTerm& operator<<(double) override final;
   LinuxTerm& operator<<(const char*) override final;
   LinuxTerm& operator<<(const std::string&) override final;
   void operator>>(std::string&) override final;
private:
   // ****************************** Functions ******************************
   void set_ops(Ops) override final;
   void reset_cursor(int);
   void reprint(bool);
   // ****************************** Constants ******************************
   static constexpr char _backspaceCh {'\x7f'};
   static constexpr char _escapeCh {'\x1b'};
   static constexpr char _arrowRightCh {'C'};
   static constexpr char _arrowLeftCh {'D'};
   static const char* _cursorUpStr;
   static const char* _boldTextStr;
   static const char* _normTextStr;
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
