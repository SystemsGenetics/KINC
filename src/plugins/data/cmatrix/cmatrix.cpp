#include "cmatrix.h"



cmatrix::cmatrix(const string& type, const string& file):
   DataPlugin(type,file),
   _mem(*KincFile::mem())
{
   if (KincFile::head()==FileMem::nullPtr)
   {
      _mem.allot(_hdr);
      KincFile::head(_hdr.addr());
      _hdr.geneSize() = 0;
      _hdr.sampleSize() = 0;
      _hdr.corrSize() = 0;
      _hdr.genePtr() = FileMem::nullPtr;
      _hdr.samplePtr() = FileMem::nullPtr;
      _hdr.corrPtr() = FileMem::nullPtr;
      _hdr.data() = FileMem::nullPtr;
      _mem.sync(_hdr,FileSync::write);
   }
   else
   {
      _hdr = KincFile::head();
      _mem.sync(_hdr,FileSync::read);
   }
}



void cmatrix::load(GetOpts&, Terminal &tm)
{
   tm << "This command is not supported by cmatrix data type.\n";
}



void cmatrix::dump(GetOpts&, Terminal &tm)
{
   tm << "This command is not supported by cmatrix data type.\n";
}



void cmatrix::query(GetOpts &ops, Terminal &tm)
{
   tm << "To be added.\n";
}



bool cmatrix::empty()
{
   return _hdr.geneSize()==0;
}



void cmatrix::set_gene_size(uint32_t size)
{
   assert<AlreadySet>(_hdr.genePtr()==FileMem::nullPtr,__FILE__,__LINE__);
   _hdr.geneSize() = size;
   GHdr ghdr;
   _mem.allot(ghdr,size);
   _hdr.genePtr() = ghdr.addr();
   ghdr.name() = FileMem::nullPtr;
   for (int i=0;i<size;++i)
   {
      _mem.sync(ghdr,FileSync::write,i);
   }
   _mem.sync(_hdr,FileSync::write);
}



void cmatrix::set_sample_size(uint32_t size)
{
   assert<AlreadySet>(_hdr.samplePtr()==FileMem::nullPtr,__FILE__,__LINE__);
   _hdr.sampleSize() = size;
   SHdr shdr;
   _mem.allot(shdr,size);
   _hdr.samplePtr() = shdr.addr();
   shdr.name() = FileMem::nullPtr;
   for (int i=0;i<size;++i)
   {
      _mem.sync(shdr,FileSync::write,i);
   }
   _mem.sync(_hdr,FileSync::write);
}



void cmatrix::set_correlation_size(uint32_t size)
{
   assert<AlreadySet>(_hdr.corrPtr()==FileMem::nullPtr,__FILE__,__LINE__);
   _hdr.corrSize() = size;
   CHdr chdr;
   _mem.allot(chdr,size);
   _hdr.corrPtr() = chdr.addr();
   chdr.name() = FileMem::nullPtr;
   for (int i=0;i<size;++i)
   {
      _mem.sync(chdr,FileSync::write,i);
   }
   _mem.sync(_hdr,FileSync::write);
}



void cmatrix::set_gene_name(uint32_t n, const string& name)
{
   assert<OutOfRange>(n<_hdr.geneSize(),__FILE__,__LINE__);
   FString newName(KincFile::mem());
   newName = name;
   GHdr ghdr(_hdr.genePtr());
   ghdr.name() = newName.addr();
   _mem.sync(ghdr,FileSync::write,n);
}



void cmatrix::set_sample_name(uint32_t n, const string& name)
{
   assert<OutOfRange>(n<_hdr.sampleSize(),__FILE__,__LINE__);
   FString newName(KincFile::mem());
   newName = name;
   SHdr shdr(_hdr.samplePtr());
   shdr.name() = newName.addr();
   _mem.sync(shdr,FileSync::write,n);
}



void cmatrix::set_correlation_name(uint32_t n, const string& name)
{
   assert<OutOfRange>(n<_hdr.corrSize(),__FILE__,__LINE__);
   FString newName(KincFile::mem());
   newName = name;
   CHdr chdr(_hdr.corrPtr());
   chdr.name() = newName.addr();
   _mem.sync(chdr,FileSync::write,n);
}



cmatrix::MIterator cmatrix::set_modes(uint32_t,uint32_t,uint8_t)
{
   ;
}
