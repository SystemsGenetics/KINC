#include "ematrix.h"
#include <sstream>
#include <cmath>



ematrix::ematrix(const string& type, const string& file):
   DataPlugin(type,file),
   _mem(*KincFile::mem())
{
   if (KincFile::head()==FileMem::nullPtr)
   {
      _mem.allot(_hdr);
      KincFile::head(_hdr.addr());
      _hdr.sampleSize() = 0;
      _hdr.geneSize() = 0;
      _hdr.transform() = Transform::none;
      _hdr.samplePtr() = FileMem::nullPtr;
      _hdr.genePtr() = FileMem::nullPtr;
      _hdr.expPtr() = FileMem::nullPtr;
      _mem.sync(_hdr,FileSync::write);
   }
   else
   {
      _hdr = KincFile::head();
      _mem.sync(_hdr,FileSync::read);
   }
}



void ematrix::load(GetOpts &ops, Terminal &tm)
{
   assert<NotNewFile>(head()==FileMem::nullPtr,__FILE__,__LINE__);
   ifile f(ops.com_front());
   assert<CannotOpen>(f.is_open(),__FILE__,__LINE__);
   _mem.allot(_hdr);
   try
   {
      load_samples(tm,f);
      load_genes(tm,f,"NaN");
   }
   catch (...)
   {
      KincFile::clear();
      throw;
   }
   KincFile::head(_hdr.addr());
   _mem.sync(_hdr,FileSync::write);
   tm << "Loaded " << _hdr.geneSize() << " gene(s) with " << _hdr.sampleSize()
      << " sample(s) each.\n";
}



void ematrix::load_samples(Terminal& tm, ifile& f)
{
   while (f.peek()==' '||f.peek()=='\t'||f.peek()=='\n')
   {
      assert<InvalidFile>(f,__FILE__,__LINE__);
      f.get();
   }
   tm << "Loading sample headers...\n";
   string buf;
   std::getline(f,buf);
   std::istringstream ibuf(buf);
   string tmp;
   std::vector<string> vs;
   while (ibuf >> tmp)
   {
      vs.push_back(std::move(tmp));
   }
   assert<InvalidFile>(vs.size()>0,__FILE__,__LINE__);
   tm << "Writing sample headers...\n";
   sHdr shdr;
   _mem.allot(shdr,vs.size());
   _hdr.sampleSize() = vs.size();
   _hdr.samplePtr() = shdr.addr();
   for (int i=0;i<vs.size();++i)
   {
      FString tname(&_mem);
      tname = vs[i];
      shdr.name() = tname.addr();
      _mem.sync(shdr,FileSync::write,i);
   }
}



void ematrix::load_genes(Terminal& tm, ifile& f, string nan)
{
   tm << "Loading genes and expression data...\n";
   std::vector<string> gns;
   std::vector<std::vector<float>> exps;
   while (f)
   {
      string buf;
      std::getline(f,buf);
      if (!buf.empty())
      {
         std::istringstream ibuf(buf);
         string tmp;
         bool cond = ibuf >> tmp;
         assert<InvalidFile>(cond,__FILE__,__LINE__);
         gns.push_back(std::move(tmp));
         exps.push_back({});
         int count = 0;
         while (ibuf >> tmp)
         {
            if (tmp==nan)
            {
               exps.back().push_back(NAN);
            }
            else
            {
               try
               {
                  exps.back().push_back(stof(tmp));
               }
               catch (std::exception)
               {
                  throw InvalidFile(__FILE__,__LINE__);
               }
            }
            ++count;
         }
         assert<InvalidFile>(count==_hdr.sampleSize(),__FILE__,__LINE__);
      }
   }
   assert<InvalidFile>(gns.size()>0,__FILE__,__LINE__);
   tm << "Writing gene headers...\n";
   gHdr ghdr;
   _mem.allot(ghdr,gns.size());
   _hdr.geneSize() = gns.size();
   _hdr.genePtr() = ghdr.addr();
   for (int i=0;i<gns.size();++i)
   {
      FString tname(&_mem);
      tname = gns[i];
      ghdr.name() = tname.addr();
      _mem.sync(ghdr,FileSync::write,i);
   }
   tm << "Writing expression data...\n";
   Exp exp;
   _mem.allot(exp,exps.size()*exps.back().size());
   for (int x=0;x<exps.size();++x)
   {
      for (int y=0;y<exps[x].size();++y)
      {
         exp.val() = exps[x][y];
         _mem.sync(exp,FileSync::write,x*_hdr.sampleSize()+y);
      }
   }
}
