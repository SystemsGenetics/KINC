#include <regex>
#include "getopts.h"



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



int GetOpts::com_get(std::initializer_list<string> commands)
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
