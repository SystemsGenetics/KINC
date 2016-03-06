#include <cstring>
#include "kincfile.h"



KincFile::KincFile(const std::string& fileName)
{
   _mem = fptr(new FileMem(fileName));
   if (_mem->size()==0)
   {
      _mem->allot(_header);
      _hist = hptr(new History(*_mem));
      _header.histHead() = _hist->addr();
      _header.dataHead() = FileMem::nullPtr;
      _header.ident() = FileMem::nullPtr;
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
      cond = strncmp(_header.idString(),_idString,_idSz)==0;
      assert<InvalidFile>(cond,__FILE__,__LINE__);
      _hist = hptr(new History(*_mem,_header.histHead()));
      _ident.addr(_header.ident());
      _new = false;
   }
}



bool KincFile::is_new()
{
   return _new;
}



History& KincFile::history()
{
   return *_hist;
}



KincFile::string KincFile::ident() const
{
   return *_ident;
}



void KincFile::ident(const string& id)
{
   try
   {
      _ident = id;
   }
   catch (FString::AlreadySet)
   {
      throw AlreadySet(__FILE__,__LINE__);
   }
   _header.ident() = _ident.addr();
   _mem->sync(_header,FileSync::write);
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
