#include <cstring>
#include "kincfile.h"



/// @brief Open KINC file object.
///
/// Open file memory object from file given. If the size of the file memory
/// object is zero then create a new KINC file from it, else attempt to load
/// KINC file from head location of file memory object.
///
/// @param fileName Location of file for memory object.
///
/// @exception InvalidFile The file being opened is not a valid KINC file.
KincFile::KincFile(const std::string& fileName):
   _mem(fileName),
   _ident(&_mem)
{
   if (_mem.size()==0)
   {
      create();
   }
   else
   {
      bool cond = _mem.size()>=_hdrSz;
      assert<InvalidFile>(cond,__FILE__,__LINE__);
      _hdr = _mem.head();
      _mem.sync(_hdr,FileSync::read);
      cond = _hdr.histHead()!=FileMem::nullPtr;
      assert<InvalidFile>(cond,__FILE__,__LINE__);
      cond = true;
      for (int i=0;i<_idSz;++i)
      {
         if (_hdr.idString()[i]!=_idString[i])
         {
            cond = false;
         }
      }
      assert<InvalidFile>(cond,__FILE__,__LINE__);
      _hist = hptr(new History(_mem,_hdr.histHead()));
      _ident.addr(_hdr.ident());
      _new = false;
   }
}



/// @brief Clear this object.
///
/// Clear all data in this object, including all data plugin memory, making this
/// object a new KINC file with no information.
void KincFile::clear()
{
   _hist.reset();
   _new = true;
   _hdr = FileMem::nullPtr;
   _ident.addr(FileMem::nullPtr);
   _mem.clear();
   create();
}



/// Set value for data plugin ident.
///
/// @param id Value for ident.
///
/// @exception AlreadySet The ident value for data plugin has already been set.
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
   _hdr.ident() = _ident.addr();
   _mem.sync(_hdr,FileSync::write);
}



/// Set value of file memory location for beginning of data plugin memory.
///
/// @param ptr Location of data plugin memory.
void KincFile::head(FileMem::Ptr ptr)
{
   _hdr.dataHead() = ptr;
   _mem.sync(_hdr,FileSync::write);
}



/// @brief Create new KINC file.
///
/// Create a new KINC file for this object, writing to this object's file memory
/// object.
void KincFile::create()
{
   _mem.allot(_hdr);
   _hist = hptr(new History(_mem));
   _hdr.histHead() = _hist->addr();
   _hdr.dataHead() = FileMem::nullPtr;
   _hdr.ident() = FileMem::nullPtr;
   memcpy(_hdr.idString(),_idString,_idSz);
   _mem.sync(_hdr,FileSync::write);
}
