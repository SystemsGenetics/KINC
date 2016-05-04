#ifndef CMATRIX_HH
#define CMATRIX_HH
#include "../../../dataplugin.h"



namespace cmatrixData
{

struct Header : public FileMem::Static<44>
{
   using FPtr = FileMem::Ptr;
   using FileMem::Static<44>::Static;
   uint32_t& geneSize() { get<uint32_t>(0); }
   uint32_t& sampleSize() { get<uint32_t>(4); }
   uint32_t& corrSize() { get<uint32_t>(8); }
   FPtr& genePtr() { get<FPtr>(12); }
   FPtr& samplePtr() { get<FPtr>(20); }
   FPtr& corrPtr() { get<FPtr>(28); }
   FPtr& data() { get<FPtr>(36); }
};

struct GeneHdr : public FileMem::Static<8>
{
   using FPtr = FileMem::Ptr;
   using FileMem::Static<8>::Static;
   FPtr& name() { get<FPtr>(0); }
};

struct SampleHdr : public FileMem::Static<8>
{
   using FPtr = FileMem::Ptr;
   using FileMem::Static<8>::Static;
   FPtr& name() { get<FPtr>(0); }
};

struct CorrHdr : public FileMem::Static<8>
{
   using FPtr = FileMem::Ptr;
   using FileMem::Static<8>::Static;
   FPtr& name() { get<FPtr>(0); }
};

struct Correlation : public FileMem::Static<12>
{
   using FPtr = FileMem::Ptr;
   using FileMem::Static<12>::Static;
   uint32_t& modeSize() { get<uint32_t>(4); }
   FPtr& modePtr() { get<FPtr>(4); }
};

struct Mode : FileMem::Object
{
   using FPtr = FileMem::Ptr;
   Mode(int samples, int corrs, FPtr ptr = FileMem::nullPtr):
      Object((samples/8)+((samples%8)!=0)+4*corrs,ptr),
      _sBytes((samples/8)+((samples%8)!=0))
   {}
   uint8_t mask(int n)
   {
      int byte = n/8;
      n%=8;
      uint8_t tmp = get<uint8_t>(byte);
      return (tmp>>n)&0x1;
   }
   void mask(int n, uint8_t v)
   {
      int byte = n/8;
      n%=8;
      switch (v)
      {
      case 0:
         get<uint8_t>(byte)&(~(0x1<<n));
         break;
      case 1:
         get<uint8_t>(byte)|(0x1<<n);
         break;
      }
   }
   float& corr(int n) { get<float>(4*n+_sBytes); }
private:
   int _sBytes;
};

}



class cmatrix : public DataPlugin
{
public:
   cmatrix(const string& type, const string& file): DataPlugin(type,file) {}
   void load(GetOpts &ops, Terminal &tm) override final {}
   void dump(GetOpts &ops, Terminal &tm) override final {}
   void query(GetOpts &ops, Terminal &tm) override final {}
   bool empty() override final { return true; }
};



#endif
