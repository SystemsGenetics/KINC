#include "ematrix.h"
#include <sstream>



ematrix::ematrix(const string& type, const string& file):
   DataPlugin(type,file),
   _mem(*KincFile::mem())
{
   if (KincFile::is_new())
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
   if (!load_samples(f)||!load_genes(f)||!load_data(f))
   {
      ;//something failed.
   }
   ;//successfully loaded.
}



bool ematrix::load_samples(ifile& f)
{
   while (f.peek()==' '||f.peek()=='\t'||f.peek()=='\n')
   {
      f.get();
   }
   string buf;
   std::getline(f,buf);
   std::istringstream ibuf(buf);
   string tmp;
   std::vector<string> vs;
   while (ibuf >> tmp)
   {
      vs.push_back(std::move(tmp));
   }
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
