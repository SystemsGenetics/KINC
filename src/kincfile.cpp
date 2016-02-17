#include "kincfile.h"



KincFile::KincFile(const std::string& fileName):
   _mem(nullptr),
   _hist(nullptr)
{
   _mem = new FileMem(fileName);
   if (_mem->size()==0)
   {
      _mem->allot(_header);
      _hist = new History(*_mem);
      _header.dataHead() = FileMem::nullPtr;
      _header.histHead() = _hist->addr();
      _mem->sync(_header,FileSync::write);
   }
   else
   {
      bool cond = _mem->size()>=_hdrSz;
      assert<InvalidFile>(cond,__FILE__,__LINE__);
      _header = _mem->head();
      _mem->sync(_header,FileSync::read);
      cond = _header.histHead()!=FileMem::nullPtr;
      assert<InvalidFile>(cond,__FILE__,__LINE__);
      cond = true;
      for (int i=0;i<_idSz;++i)
      {
         if (_header.idString()[i]!=_idString[i])
         {
            cond = false;
         }
      }
      assert<InvalidFile>(cond,__FILE__,__LINE__);
      _hist = new History(*_mem,_header.histHead());
   }
}



KincFile::~KincFile()
{
   if (_hist)
   {
      delete _hist;
   }
   if (_mem)
   {
      delete _mem;
   }
}



FileMem::Ptr KincFile::head() const
{
   return _header.dataHead();
}



void KincFile::head(FileMem::Ptr ptr)
{
   _header.dataHead() = ptr;
   _mem->sync(_header,FileSync::write);
}



History& KincFile::history()
{
   return *_hist;
}
