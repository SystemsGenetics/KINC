#include "fstring.h"



FString::FString(FileMem* mem, FPtr ptr):
   _mem(mem),
   _hdr(ptr)
{
   bool cond = mem!=nullptr;
   assert<InvalidPtr>(cond,__FILE__,__LINE__);
   if (_hdr.addr()!=FileMem::nullPtr)
   {
      load();
   }
}



FString::FString(FString&& tmp):
   _mem(tmp._mem),
   _hdr(tmp._hdr),
   _str(std::move(tmp._str))
{
   tmp._hdr = FileMem::nullPtr;
}



FString& FString::operator=(FString&& tmp)
{
   _mem = tmp._mem;
   _hdr = tmp._hdr;
   _str = std::move(tmp._str);
   tmp._hdr = FileMem::nullPtr;
}



FString& FString::operator=(const string& nStr)
{
   bool cond = _hdr.addr()==FileMem::nullPtr;
   assert<AlreadySet>(cond,__FILE__,__LINE__);
   String fStr(nStr.size()+1);
   _mem->allot(_hdr);
   _mem->allot(fStr);
   _hdr.stripe() = FStringData::strip;
   _hdr.sSize() = nStr.size()+1;
   memcpy(fStr.c_str(),nStr.c_str(),nStr.size()+1);
   _mem->sync(_hdr,FileSync::write);
   _mem->sync(fStr,FileSync::write);
   _str = nStr;
}



void FString::addr(FPtr ptr)
{
   _hdr = ptr;
   _str.clear();
   if (_hdr.addr()!=FileMem::nullPtr)
   {
      load();
   }
}



inline void FString::load()
{
   _mem->sync(_hdr,FileSync::read);
   bool cond = _hdr.stripe()==FStringData::strip;
   assert<InvalidPtr>(cond,__FILE__,__LINE__);
   String fStr(_hdr.sSize(),_hdr.addr()+FStringData::hdrSz);
   _mem->sync(fStr,FileSync::read);
   _str = fStr.c_str();
}
