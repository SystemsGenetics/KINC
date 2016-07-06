#ifndef CMATRIX_H
#define CMATRIX_H
#include <ace.h>
#include <vector>



class CMatrix : public AccelCompEng::DataPlugin
{
public:
   // *
   // * DECLERATIONS
   // *
   ACE_DATA_HEADER()
   using string = std::string;
   // *
   // * CLASSES
   // *
   class GeneModes;
   class GeneCorrs;
   // *
   // * EXCEPTIONS
   // *
   ACE_EXCEPTION(CMatrix,AlreadyInitialized)
   ACE_EXCEPTION(CMatrix,InvalidSize)
   ACE_EXCEPTION(CMatrix,OutOfRange)
   ACE_EXCEPTION(CMatrix,NotInitialized)
   ACE_EXCEPTION(CMatrix,InvalidGeneCorr)
   ACE_EXCEPTION(CMatrix,ValueNotRead)
   // *
   // * BASIC METHODS
   // *
   CMatrix(const string& type, const string& file);
   // *
   // * VIRTUAL FUNCTIONS
   // *
   void load(GetOpts&,Terminal&) override final;
   void dump(GetOpts&,Terminal&) override final;
   void query(GetOpts&,Terminal&) override final;
   bool empty() override final;
   // *
   // * FUNCTIONS
   // *
   void initialize(uint32_t,uint32_t,uint32_t,uint8_t);
   void set_gene_name(uint32_t,const string&);
   void set_sample_name(uint32_t,const string&);
   void set_correlation_name(uint32_t,const string&);
   GeneModes get_modes(uint32_t,uint32_t);
   void get_modes(GeneModes&,uint32_t,uint32_t);
   GeneCorrs get_corrs(uint32_t,uint32_t);
   void get_corrs(GeneCorrs&,uint32_t,uint32_t);
private:
   // *
   // * DECLERATIONS
   // *
   using SizeT = AccelCompEng::FileMem::SizeT;
   // *
   // * FILEMEM TYPES
   // *
   ACE_FMEM_STATIC(Header,53)
      ACE_FMEM_VAL(gSize,uint32_t,0)
      ACE_FMEM_VAL(sSize,uint32_t,4)
      ACE_FMEM_VAL(cSize,uint32_t,8)
      ACE_FMEM_VAL(mdMax,uint8_t,12)
      ACE_FMEM_VAL(gPtr,FPtr,13)
      ACE_FMEM_VAL(sPtr,FPtr,21)
      ACE_FMEM_VAL(cPtr,FPtr,29)
      ACE_FMEM_VAL(mdData,FPtr,37)
      ACE_FMEM_VAL(crData,FPtr,45)
   ACE_FMEM_END()
   ACE_FMEM_STATIC(NmHead,8)
      ACE_FMEM_VAL(nPtr,FPtr,0)
   ACE_FMEM_END()
   // *
   // * FUNCTIONS
   // *
   void set_name(FPtr,uint32_t,uint32_t,const string&);
   SizeT diagonal(uint32_t,uint32_t);
   SizeT diag_size(uint32_t);
   // *
   // * VARIABLES
   // *
   Header _hdr;
};



class CMatrix::GeneModes
{
public:
   friend class CMatrix;
   bool mode(uint8_t,uint32_t);
   void mode(uint8_t,uint32_t,bool);
   void mode(uint8_t,bool);
   void mode(bool);
   void read();
   void write();
private:
   GeneModes(CMatrix*,FPtr);
   void addr(FPtr);
   FPtr static initialize(CMatrix&,uint32_t,uint32_t,uint8_t);
   // *
   // * FILEMEM TYPES
   // *
   ACE_FMEM_OBJECT(Modes)
      Modes(int mdSz, int mkSz, FPtr ptr = fNullPtr):
         Object(1+mdSz*mkSz,ptr),
         _mkSz(mkSz)
      {}
      ACE_FMEM_VAL(numMds,uint8_t,0)
      uint8_t& mask(int md, int mInc) { get<uint8_t>(1+(md*_mkSz)+mInc); }
      int _mkSz;
   ACE_FMEM_END()
   CMatrix* _p;
   Modes _data;
   bool _isRead {false};
};



class CMatrix::GeneCorrs
{
public:
   friend class CMatrix;
   float corr(uint8_t,uint32_t);
   void corr(uint8_t,uint32_t,float);
   void read();
   void write();
private:
   GeneCorrs(CMatrix*,FPtr);
   void addr(FPtr);
   FPtr static initialize(CMatrix&,uint32_t,uint32_t,uint8_t);
   // *
   // * FILEMEM TYPES
   // *
   ACE_FMEM_OBJECT(Corrs)
      Corrs(int mdSz, int crSz, FPtr ptr = fNullPtr):
         Object(mdSz*crSz*4,ptr),
         _mdSz(mdSz)
      {}
      float& val(int md, int Inc) { get<float>((md*_mdSz)+Inc); }
      int _mdSz;
   ACE_FMEM_END()
   CMatrix* _p;
   Corrs _data;
   bool _isRead {false};
};



#endif
