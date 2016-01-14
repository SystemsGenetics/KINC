#ifndef LINUXTERM_H
#define LINUXTERM_H
#include <string>
#include <list>
#include "terminal.h"



class LinuxTerm : public Terminal
{
private:
   static bool _lock;
   static bool _cooked;
   int _cols;
   int _chCount;
   std::string _header;
   std::list<char> _line;
   std::list<char>::iterator _i;
   void reset_cursor(int);
   void calc_new_lines();
   void reprint(bool);
public:
   static void stty_raw();
   static void stty_cooked();
   LinuxTerm();
   ~LinuxTerm();
   void header(const std::string&);
   void set_ops(Ops);
   void precision(int);
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
};



#endif
