#include "fstring.h"



FString::FString(FileMem* mem, FPtr ptr):
   _mem(mem),
   _hdr(ptr)
{
   if (_hdr.addr()!=FileMem::nullPtr)
   {
      _mem->sync(_hdr,FileSync::read);
      bool cond = _hdr.stripe()==FStringData::strip;
      assert<InvalidPtr>(cond,__FILE__,__LINE__);
      String fStr(_hdr.sSize());
      fStr = _hdr.addr()+FStringData::hdrSz;
      _mem->sync(fStr,FileSync::read);
      _str = fStr.c_str();
   }
}



FString::FString(FString&& tmp):
   _mem(tmp._mem),
   _hdr(tmp._hdr)
{
   tmp._hdr = FileMem::nullPtr;
}



FString& FString::operator=(FString&& tmp)
{
   _mem = tmp._mem;
   _hdr = tmp._hdr;
   tmp._hdr = FileMem::nullPtr;
}



const FString::string& FString::operator*()
{
   return _str;
}



const FString::string* FString::operator->()
{
   return &_str;
}



FString& FString::operator=(const string& nStr)
{
   bool cond = _hdr.addr()==FileMem::nullPtr;
   assert<AlreadySet>(cond,__FILE__,__LINE__);
   String fStr(nStr.size());
   _mem->allot(_hdr);
   _mem->allot(fStr);
   _hdr.stripe() = FStringData::strip;
   _hdr.sSize() = nStr.size();
   memcpy(fStr.c_str(),nStr.c_str(),nStr.size());
   _mem->sync(_hdr,FileSync::write);
   _mem->sync(fStr,FileSync::write);
   _str = nStr;
}
