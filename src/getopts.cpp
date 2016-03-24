#include <regex>
#include "getopts.h"



/// @brief Decomposes user line.
///
/// Decomposes the raw user input line given into commands and options, adding
/// them to this new object's lists.
///
/// @param raw User input line.
///
/// @exception InvalidSyntax An error occured in the syntax of the user line
/// while attempting to decompose it.
GetOpts::GetOpts(const string& raw):
   _orig(raw)
{
   using regex = std::regex;
   using regtok = std::sregex_token_iterator;
   regex pat {"[\\s\\t]+"};
   regtok end {};
   for (regtok r {raw.begin(),raw.end(),pat,-1};r!=end;++r)
   {
      const string& i {*r};
      if (i.front()=='-')
      {
         string val;
         string key;
         auto x = i.find_first_not_of('-');
         auto y = i.find('=');
         if (y!=string::npos)
         {
            bool cond = y==i.rfind('=')&&y>x;
            assert<InvalidSyntax>(cond,__FILE__,__LINE__);
            val = i.substr(x,y-x);
            key = i.substr(++y);
         }
         else
         {
            val = i.substr(x);
         }
         _opts.push_back({std::move(val),std::move(key)});
      }
      else
      {
         _comms.push_back(i);
      }
   }
}



/// @brief Interrogate front command.
///
/// Interrogate the command at front of list for this object, comparing it to
/// initializer list of strings given seeing if it matches any of them.
///
/// @param commands List of values to match front of list command by.
/// @return 0 if no match is found in list of command values, else increment
/// into initializer list where match was found starting with 1 for beginning
/// of list.
int GetOpts::com_get(initlist commands)
{
   int ret = 0;
   int count = 0;
   if (!_comms.empty())
   {
      for (auto i:commands)
      {
         count++;
         if (i==_comms.front())
         {
            ret = count;
            break;
         }
      }
   }
   return ret;
}



/// Remove the command at front of list of commands for this object.
void GetOpts::com_pop()
{
   if (!_comms.empty())
   {
      _comms.pop_front();
   }
}



/// Interrogate list of options to see if the key of any option matches the
/// value given, optionally removing all matches from list of options if found.
///
/// @param opt Key value to match against.
/// @param del True to remove all matching options, else false.
/// @return True if one or more matches found, else false.
bool GetOpts::has_opt(const string& opt, bool del)
{
   bool ret {false};
   for (auto i = _opts.begin();i!=_opts.end();)
   {
      if (i->first==opt)
      {
         ret = true;
         if (del)
         {
            i = _opts.erase(i);
         }
         else
         {
            ++i;
         }
      }
      else
      {
         ++i;
      }
   }
   return ret;
}
