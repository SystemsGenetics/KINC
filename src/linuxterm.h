#ifndef LINUXTERM_H
#define LINUXTERM_H
#include <string>
#include <list>
#include "terminal.h"



class LinuxTerm : public Terminal
/*
 * Linux OS specific implemenation.
 *
 * CONDITIONS:
 * 1. Only one instance of this class can exist at anytime.
 * 2. stty_raw() must be called before constructor.
 * 3. stty_cooked() must be called after destructor.
 * 4. The terminal this program is running in cannot change its width during
 *    the lifetime of an instance.
 */
{
public:
   // ****************************** Functions ******************************
   static void stty_raw();
   static void stty_cooked();
   LinuxTerm();
   ~LinuxTerm();
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
   // ****************************** Variables ******************************
   static bool _lock;
   static bool _cooked;
   int _cols;
   int _chCount;
   std::string _header;
   std::list<char> _line;
   std::list<char>::iterator _i;
};



#endif
