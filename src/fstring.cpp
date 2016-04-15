#include "fstring.h"



/// @brief Initialize file string object.
///
/// Initializes file string, making sure file memory object pointer is valid and
/// loading the file string from memory if location of string given is not
/// nullptr.
///
/// @param mem Pointer to file memory object that is used.
/// @param ptr File pointer location where file string is location, nullptr if
/// string object is not yet set.
///
/// @exception InvalidPtr The pointer given for the file memory object is
/// nullptr.
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



/// Move given file string.
///
/// @param tmp Object to take data from.
FString::FString(FString&& tmp):
   _mem(tmp._mem),
   _hdr(tmp._hdr),
   _str(std::move(tmp._str))
{
   tmp._hdr = FileMem::nullPtr;
}



/// Move given file string, overwriting what this object currently stores.
///
/// @param tmp Object to take data from.
FString& FString::operator=(FString&& tmp)
{
   _mem = tmp._mem;
   _hdr = tmp._hdr;
   _str = std::move(tmp._str);
   tmp._hdr = FileMem::nullPtr;
}



/// Set file string to given string.
///
/// @param nStr Value to set file string to in file memory.
/// @return Reference to this object.
///
/// @exception AlreadySet This file string object has already been set.
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



/// @brief Set new location for this object.
///
/// Set new location in file memory where there is an FString to read or set
/// the adress to nullptr so a new string can be allocated.
///
/// @param ptr New location of FString for nullptr if unset.
void FString::addr(FPtr ptr)
{
   _hdr = ptr;
   _str.clear();
   if (_hdr.addr()!=FileMem::nullPtr)
   {
      load();
   }
}



/// @brief Load file string from file.
///
/// Load value of file string from file memory object from location stored in
/// this object.
///
/// @exception InvalidPtr The file memory location is not a valid file string.
inline void FString::load()
{
   _mem->sync(_hdr,FileSync::read);
   bool cond = _hdr.stripe()==FStringData::strip;
   assert<InvalidPtr>(cond,__FILE__,__LINE__);
   String fStr(_hdr.sSize(),_hdr.addr()+FStringData::hdrSz);
   _mem->sync(fStr,FileSync::read);
   _str = fStr.c_str();
}
