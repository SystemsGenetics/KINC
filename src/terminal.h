#ifndef TERMINAL_H
#define TERMINAL_H
#include <string>
#include <list>



class Terminal
{
private:
   static bool __lock;
   int _cols;
   int _chCount;
   std::string _header;
   std::list<char> _line;
   std::list<char>::iterator _i;
   static void stty_raw();
   static void stty_cooked();
   void reset_cursor(int);
   void calc_new_lines();
   void reprint(bool);
public:
   Terminal();
   ~Terminal();
   void reset_line();
   std::string get();
   void clear(const char*);
   bool read();
   Terminal& operator<<(char);
};



#endif
