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
   AccelCompEng::assert<AlreadySet>(_hdr.genePtr()==FileMem::nullPtr,__LINE__);
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
   AccelCompEng::assert<AlreadySet>(_hdr.samplePtr()==FileMem::nullPtr,__LINE__);
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
   AccelCompEng::assert<AlreadySet>(_hdr.corrPtr()==FileMem::nullPtr,__LINE__);
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
   AccelCompEng::assert<OutOfRange>(n<_hdr.geneSize(),__LINE__);
   FString newName(KincFile::mem());
   newName = name;
   GHdr ghdr(_hdr.genePtr());
   ghdr.name() = newName.addr();
   _mem.sync(ghdr,FileSync::write,n);
}



void cmatrix::set_sample_name(uint32_t n, const string& name)
{
   AccelCompEng::assert<OutOfRange>(n<_hdr.sampleSize(),__LINE__);
   FString newName(KincFile::mem());
   newName = name;
   SHdr shdr(_hdr.samplePtr());
   shdr.name() = newName.addr();
   _mem.sync(shdr,FileSync::write,n);
}



void cmatrix::set_correlation_name(uint32_t n, const string& name)
{
   AccelCompEng::assert<OutOfRange>(n<_hdr.corrSize(),__LINE__);
   FString newName(KincFile::mem());
   newName = name;
   CHdr chdr(_hdr.corrPtr());
   chdr.name() = newName.addr();
   _mem.sync(chdr,FileSync::write,n);
}



void cmatrix::create_data()
{
   AccelCompEng::assert<AlreadySet>(_hdr.data()==FileMem::nullPtr,__LINE__);
   Corr corr;
   _mem.allot(corr,_hdr.geneSize()*(_hdr.geneSize()-1)/2);
   _hdr.data() = corr.addr();
}



cmatrix::FPtr cmatrix::set_modes(uint32_t g1, uint32_t g2, uint8_t modes)
{
   AccelCompEng::assert<NotInitialized>(_hdr.data()!=FileMem::nullPtr,__LINE__);
   AccelCompEng::assert<InvalidGeneCorr>(g1!=g2,__LINE__);
   if (g1>g2)
   {
      uint32_t tmp {g1};
      g1 = g2;
      g2 = tmp;
   }
   FileMem::SizeT inc {(g1*(g1+1)/2)+g2};
   Corr corr(_hdr.data());
   Md md(_hdr.sampleSize(),_hdr.corrSize());
   _mem.allot(md,modes);
   corr.modePtr() = md.addr();
   corr.modeSize() = modes;
   _mem.sync(corr,FileSync::write,inc);
   return md.addr();
}



cmatrix::FPtr cmatrix::get_modes(uint32_t g1, uint32_t g2)
{
   AccelCompEng::assert<NotInitialized>(_hdr.data()!=FileMem::nullPtr,__LINE__);
   AccelCompEng::assert<InvalidGeneCorr>(g1!=g2,__LINE__);
   if (g1>g2)
   {
      uint32_t tmp {g1};
      g1 = g2;
      g2 = tmp;
   }
   FileMem::SizeT inc {(g1*(g1+1)/2)+g2};
   Corr corr(_hdr.data());
   _mem.sync(corr,FileSync::read,inc);
   return corr.modePtr();
}



void cmatrix::write_mode(FPtr ptr, uint8_t inc, const maskv& mask,
                         const floatv& corrs)
{
   AccelCompEng::assert<InvalidSize>(mask.size()==_hdr.sampleSize(),__LINE__);
   AccelCompEng::assert<InvalidSize>(corrs.size()==_hdr.corrSize(),__LINE__);
   Md md(_hdr.sampleSize(),_hdr.corrSize(),ptr);
   int x {0};
   for (auto i:mask)
   {
      md.mask(x++,i);
   }
   x = 0;
   for (auto i:corrs)
   {
      md.corr(x) = i;
   }
   _mem.sync(md,FileSync::write,inc);
}
