#ifndef CMATRIX_H
#define CMATRIX_H
#include <ace.h>
#include <vector>



namespace CMatrixData
{

using namespace AccelCompEng;

struct Header : public FileMem::Static<52>
{
   using FPtr = FileMem::Ptr;
   using FileMem::Static<52>::Static;
   uint32_t& geneSize() { get<uint32_t>(0); }
   uint32_t& sampleSize() { get<uint32_t>(4); }
   uint32_t& corrSize() { get<uint32_t>(8); }
   FPtr& genePtr() { get<FPtr>(12); }
   FPtr& samplePtr() { get<FPtr>(20); }
   FPtr& corrPtr() { get<FPtr>(28); }
   FPtr& modeData() { get<FPtr>(36); }
   FPtr& corrData() { get<FPtr>(44); }
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

struct Correlation : public FileMem::Static<4>
{
   using FileMem::Static<4>::Static;
   float& val() { get<float>(0); }
};

struct Correlations : public FileMem::Object
{
   using FPtr = FileMem::Ptr;
   Correlations(int size, FPtr ptr = FileMem::nullPtr):\
      FileMem::Object(size*4,ptr) {}
   float& val(int n) { get<float>(4*n); }
};

struct Mode : public FileMem::Object
{
   using FPtr = FileMem::Ptr;
   Mode(int size, FPtr ptr = FileMem::nullPtr): FileMem::Object(size,ptr) {}
   uint8_t& mask(int n) { get<uint8_t>(n); }
};

struct Modes : public FileMem::Object
{
   using FPtr = FileMem::Ptr;
   Modes(int size, int msize, FPtr ptr = FileMem::nullPtr):
      FileMem::Object(size*msize,ptr),_Size(size) {}
   uint8_t& mask(int mask, int n) { get<uint8_t>((mask*size)+n); }
private:
   int _Size;
};

}



class CMatrix : public AccelCompEng::DataPlugin
{
public:
   // *
   // * DECLERATIONS
   // *
   using string = std::string;
   using Hdr = CMatrixData::Header;
   using GHdr = CMatrixData::GeneHdr;
   using SHdr = CMatrixData::SampleHdr;
   using CHdr = CMatrixData::CorrHdr;
   using Corr = CMatrixData::Correlation;
   using Corrs = CMatrixData::Correlations;
   using Mode = CMatrixData::Mode;
   using Modes = CMatrixData::Modes;
   using GetOpts = AccelCompEng::GetOpts;
   using Terminal = AccelCompEng::Terminal;
   using FileMem = AccelCompEng::FileMem;
   using FileSync = AccelCompEng::FileSync;
   using File = AccelCompEng::File;
   using FString = AccelCompEng::FString;
   using FPtr = FileMem::Ptr;
   // *
   // * CLASSES
   // *
   class Mask;
   // *
   // * BASIC METHODS
   // *
   //ACE_EXCEPTION(cmatrix,AlreadySet)
   //ACE_EXCEPTION(cmatrix,OutOfRange)
   //ACE_EXCEPTION(cmatrix,NotInitialized)
   //ACE_EXCEPTION(cmatrix,InvalidGeneCorr)
   //ACE_EXCEPTION(cmatrix,InvalidSize)
   // *
   // * BASIC METHODS
   // *
   CMatrix(const string& type, const string& file);
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
   void initialize(uint32_t,uint32_t,uint32_t,uint8_t);
   void set_gene_name(uint32_t,const string&);
   void set_sample_name(uint32_t,const string&);
   void set_correlation_name(uint32_t,const string&);
private:
   // *
   // * VARIABLES
   // *
   Hdr _hdr;
};



#endif
