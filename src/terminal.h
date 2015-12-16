#ifndef TERMINAL_H
#define TERMINAL_H
#include <string>
#include <list>



class Terminal
{
private:
   int _cols;
   int _chCount;
   std::string _header;
   std::list<char> _line;
   std::list<char>::iterator _i;
   void reset_cursor(int);
   void calc_new_lines();
   void reprint(bool);
public:
   Terminal();
   void reset_line();
   std::string get();
   void clear(const char*);
   Terminal& operator<<(char);
   static void stty_raw();
   static void stty_cooked();
};



#endif
