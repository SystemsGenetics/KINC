#include "getopts.h"



GetOpts::GetOpts(const string& raw):
   _orig(raw)
{
   clist rawList = explode(raw);
   for (auto i:rawList)
   {
      if (i.front()=='-')
      {
         auto x = i.begin();
         auto y = i.begin();
         for (;x!=i.end()&&*x=='-';++x);
         for (;y!=i.end()&&*y!='=';++y);
         if (y!=i.end())
         {
            ++y;
         }
         string val(y,i.end());
         string key(x,y);
         _opts.push_back({key,val});
      }
      else
      {
         _comms.push_back(std::move(i));
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



inline GetOpts::clist GetOpts::explode(const string& raw)
{
   clist list;
   string buffer;
   auto i = raw.begin();
   while (true)
   {
      if (*i!=' '&&*i!='\t'&&i!=raw.end())
      {
         buffer += *i;
      }
      else if ((*i==' '||*i=='\t'||i==raw.end())&&!buffer.empty())
      {
         list.push_back(std::move(buffer));
      }
      if (i==raw.end())
      {
         break;
      }
      ++i;
   }
   return list;
}



GetOpts::Iterator& GetOpts::Iterator::operator>>(short& n)
{
   istring buffer(_i->second);
   bool cond = buffer >> n;
   assert<InvalidType>(cond,__FILE__,__LINE__);
   return *this;
}



GetOpts::Iterator& GetOpts::Iterator::operator>>(unsigned short& n)
{
   istring buffer(_i->second);
   bool cond = buffer >> n;
   assert<InvalidType>(cond,__FILE__,__LINE__);
   return *this;
}



GetOpts::Iterator& GetOpts::Iterator::operator>>(int& n)
{
   istring buffer(_i->second);
   bool cond = buffer >> n;
   assert<InvalidType>(cond,__FILE__,__LINE__);
   return *this;
}



GetOpts::Iterator& GetOpts::Iterator::operator>>(unsigned int& n)
{
   istring buffer(_i->second);
   bool cond = buffer >> n;
   assert<InvalidType>(cond,__FILE__,__LINE__);
   return *this;
}



GetOpts::Iterator& GetOpts::Iterator::operator>>(long& n)
{
   istring buffer(_i->second);
   bool cond = buffer >> n;
   assert<InvalidType>(cond,__FILE__,__LINE__);
   return *this;
}



GetOpts::Iterator& GetOpts::Iterator::operator>>(unsigned long& n)
{
   istring buffer(_i->second);
   bool cond = buffer >> n;
   assert<InvalidType>(cond,__FILE__,__LINE__);
   return *this;
}



GetOpts::Iterator& GetOpts::Iterator::operator>>(float& n)
{
   istring buffer(_i->second);
   bool cond = buffer >> n;
   assert<InvalidType>(cond,__FILE__,__LINE__);
   return *this;
}



GetOpts::Iterator& GetOpts::Iterator::operator>>(double& n)
{
   istring buffer(_i->second);
   bool cond = buffer >> n;
   assert<InvalidType>(cond,__FILE__,__LINE__);
   return *this;
}



GetOpts::Iterator& GetOpts::Iterator::operator>>(string& n)
{
   istring buffer(_i->second);
   bool cond = buffer >> n;
   assert<InvalidType>(cond,__FILE__,__LINE__);
   return *this;
}
