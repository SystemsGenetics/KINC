#ifndef LINUXTERM_H
#define LINUXTERM_H
#include <string>
#include <list>
#include "terminal.h"



/// @brief Linux OS Terminal implementation.
///
/// Implements terminal interface for user in linux terminal environment. Uses
/// system calls to get the width of the terminal and VT100 and ANSI escape
/// codes to control the terminal in raw mode to create an interactive terminal
/// environment. See https://en.wikipedia.org/wiki/ANSI_escape_code for more
/// information about control characters used.
///
/// @pre Only one instance of this class can exist at any one time.
/// @pre stty_raw() must be called before instantiating this class.
/// @post stty_cooked() must be called after destroying this class.
/// @warning The terminal this program is running in cannot change its width
/// during the lifetime of an instance.
///
/// @author Josh Burns
/// @date 23 Jan 2016
class LinuxTerm : public Terminal
{
public:
   // *
   // * BASIC METHODS
   // *
   LinuxTerm(const LinuxTerm&) = delete;
   LinuxTerm(LinuxTerm&&) = delete;
   LinuxTerm& operator=(const LinuxTerm&) = delete;
   LinuxTerm& operator=(LinuxTerm&&) = delete;
   LinuxTerm();
   ~LinuxTerm() override final;
   // *
   // * STATIC FUNCTIONS
   // *
   static void stty_raw();
   static void stty_cooked();
   // *
   // * FUNCTIONS
   // *
   void header(const std::string&) override final;
   void precision(int) override final;
   // *
   // * OPERATORS
   // *
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
   // *
   // * FUNCTIONS
   // *
   void set_ops(Ops) override final;
   void reset_cursor(int);
   void reprint(bool);
   // *
   // * CONSTANTS
   // *
   static constexpr char _backspaceCh {'\x7f'};
   static constexpr char _escapeCh {'\x1b'};
   static constexpr char _arrowRightCh {'C'};
   static constexpr char _arrowLeftCh {'D'};
   static const char* _cursorUpStr;
   static const char* _boldTextStr;
   static const char* _normTextStr;
   // *
   // * STATIC VARIABLES
   // *
   static bool _lock;
   static bool _cooked;
   // *
   // * VARIABLES
   // *
   /// The width of the user terminal in characters.
   int _cols;
   /// The total number of characters that has been echoed to output from user
   /// input.
   int _chCount;
   /// The header that will be printed at the beginning of each user input line.
   std::string _header;
   /// Stores current line being read from the user terminal.
   std::list<char> _line;
   /// Stores where the cursor of the current input line is positioned within
   /// the temporary line of new input.
   std::list<char>::iterator _i;
};



#endif
