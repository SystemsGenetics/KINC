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

struct Correlation : public FileMem::Static<9>
{
   using FPtr = FileMem::Ptr;
   using FileMem::Static<9>::Static;
   uint8_t& modeSize() { get<uint8_t>(0); }
   FPtr& modePtr() { get<FPtr>(1); }
};

struct Mode : public FileMem::Object
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
      if (v==0)
      {
         get<uint8_t>(byte)&(~(0x1<<n));
      }
      else
      {
         get<uint8_t>(byte)|(0x1<<n);
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
   // *
   // * DECLERATIONS
   // *
   using Hdr = cmatrixData::Header;
   using GHdr = cmatrixData::GeneHdr;
   using SHdr = cmatrixData::SampleHdr;
   using CHdr = cmatrixData::CorrHdr;
   using Corr = cmatrixData::Correlation;
   using Md = cmatrixData::Mode;
   using string = std::string;
   using floatv = std::vector<float>;
   using maskv = std::vector<uint8_t>;
   using FPtr = FileMem::Ptr;
   // *
   // * BASIC METHODS
   // *
   DATA_EXCEPTION(cmatrix,AlreadySet,fatal)
   DATA_EXCEPTION(cmatrix,OutOfRange,fatal)
   DATA_EXCEPTION(cmatrix,NotInitialized,fatal)
   DATA_EXCEPTION(cmatrix,InvalidGeneCorr,fatal)
   DATA_EXCEPTION(cmatrix,InvalidSize,fatal)
   // *
   // * BASIC METHODS
   // *
   cmatrix(const string& type, const string& file);
   // *
   // * VIRTUAL FUNCTIONS
   // *
   void load(GetOpts &ops, Terminal &tm) override final;
   void dump(GetOpts &ops, Terminal &tm) override final;
   void query(GetOpts &ops, Terminal &tm) override final;
   bool empty() override final;
   // *
   // * FUNCTIONS
   // *
   void set_gene_size(uint32_t);
   void set_sample_size(uint32_t);
   void set_correlation_size(uint32_t);
   void set_gene_name(uint32_t,const string&);
   void set_sample_name(uint32_t,const string&);
   void set_correlation_name(uint32_t,const string&);
   FPtr set_modes(uint32_t,uint32_t,uint8_t);
   FPtr get_modes(uint32_t,uint32_t);
   void write_mode(FPtr,uint8_t,const maskv&,const floatv&);
private:
   // *
   // * VARIABLES
   // *
   Hdr _hdr;
   FileMem& _mem;
};



#endif
