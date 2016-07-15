#ifndef CMATRIX_H
#define CMATRIX_H
#include <ace.h>
#include <vector>



class CMatrix : public AccelCompEng::DataPlugin
{
public:
   using string = std::string;
   ACE_DATA_HEADER()
   //ACE_EXCEPTION(CMatrix,AlreadyInitialized)
   //ACE_EXCEPTION(CMatrix,InvalidSize)
   ACE_EXCEPTION(CMatrix,OutOfRange)
   //ACE_EXCEPTION(CMatrix,NotInitialized)
   //ACE_EXCEPTION(CMatrix,InvalidGeneCorr)
   //ACE_EXCEPTION(CMatrix,ValueNotRead)
   //class GeneModes;
   //class GeneCorrs;
   class GPair;
   CMatrix(const string& type, const string& file);
   void load(GetOpts&,Terminal&) override final;
   void dump(GetOpts&,Terminal&) override final;
   void query(GetOpts&,Terminal&) override final;
   bool empty() override final;
   //void initialize(uint32_t,uint32_t,uint32_t,uint8_t);
   //void set_gene_name(uint32_t,const string&);
   //void set_sample_name(uint32_t,const string&);
   //void set_correlation_name(uint32_t,const string&);
   //GeneModes get_modes(uint32_t,uint32_t);
   //void get_modes(GeneModes&,uint32_t,uint32_t);
   //GeneCorrs get_corrs(uint32_t,uint32_t);
   //void get_corrs(GeneCorrs&,uint32_t,uint32_t);
   GPair begin();
   GPair end();
   GPair& at(int,int);
   GPair& ref(int,int);
private:
   ACE_FMEM_STATIC(Header,45)
      ACE_FMEM_VAL(gSize,int32_t,0)
      ACE_FMEM_VAL(sSize,int32_t,4)
      ACE_FMEM_VAL(cSize,int32_t,8)
      ACE_FMEM_VAL(mSize,int8_t,12)
      ACE_FMEM_VAL(gPtr,FPtr,13)
      ACE_FMEM_VAL(sPtr,FPtr,21)
      ACE_FMEM_VAL(cPtr,FPtr,29)
      ACE_FMEM_VAL(eData,FPtr,37)
   ACE_FMEM_END()
   ACE_FMEM_STATIC(NmHead,8)
      ACE_FMEM_VAL(nPtr,FPtr,0)
   ACE_FMEM_END()
   //void set_name(FPtr,uint32_t,uint32_t,const string&);
   static long long diagonal(int,int);
   static long long diag_size(int);
   Header _hdr;
   GPair* _iRelation {nullptr};
};



class CMatrix::GPair//
{
   friend class CMatrix;
   class Modes;
   class Corrs;
   ~GPair();
   void read();
   void write();
   Modes& modes();
   Correlations& correlations();
   void operator++();
   bool operator!=(const GPair&);
private:
   GPair(CMatrix*,int,int);
   void set(int,int);
   ACE_FMEM_OBJECT(Relation)
      Relation(int mSize, int sSize, int cSize, Fptr ptr = fNullPtr):
         Object((mSize*sSize)+(4*cSize*mSize),ptr),
         _mSize(mSize),
         _sSize(sSize),
         _cSize(cSize)
      {}
      int _mSize;
      int _sSize;
      int _cSize;
      int8_t& mVal(int mi, int i) { get<int8_t>((_sSize*mi)+i); }
      float& cVal(int mi, int ci) { get<float>((_mSize*_sSize)+(4*((_cSize*mi)+ci))); }
   ACE_FMEM_END()
   CMatrix* _p;
   Relation _data;
   int _x;
   int _y;
   Modes* _iModes {nullptr};
   Corrs* _iCorrs {nullptr};
};



class CMatrix::GPair::Modes//
{
   class Mode;
   ~Modes();
   Modes(const Modes&) = delete;
   Modes& operator=(const Modes&) = delete;
   Modes(Modes&&) = delete;
   Modes& operator=(Modes&&) = delete;
   Mode begin();
   Mode end();
   Mode& at(int);
   Mode& operator[](int);
private:
   Modes(GPair*);
   GPair* _p;
   Mode* _iMode {nullptr};
};



class CMatrix::GPair::Modes::Mode//
{
   class Iterator;
   Iterator begin();
   Iterator end();
   int8_t& at(int);
   int8_t& operator[](int);
   void operator++();
   bool operator!=(const Mode&);
private:
   Mode(Modes*,int);
   void set(int);
   Modes* _p;
   int _i;
};



class CMatrix::GPair::Modes::Mode::Iterator//
{
   int8_t& operator*();
   void operator++();
   bool operator!=(const Iterator&);
private:
   Iterator(Mode*,int);
   Mode* _p;
   int _i;
};



class CMatrix::GPair::Corrs
{
   class Corr;
   ~Corrs();
   Corrs(const Corrs&) = delete;
   Corrs& operator=(const Corrs&) = delete;
   Corrs(Corrs&&) = delete;
   Corrs& operator=(Corrs&&) = delete;
   Corr begin();
   Corr end();
   Corr& at(int);
   Corr& operator[](int);
private:
   Corrs(GPair*);
   GPair* _p;
   Corr* _iCorr {nullptr};
};



class CMatrix::GPair::Corrs::Corr
{
   class Iterator;
   Iterator begin();
   Iterator end();
   float& at(int);
   float& operator[](int);
   void operator++();
   bool operator!=(const Mode&);
private:
   Corr(Corrs*,int);
   void set(int);
   Corrs* _p;
   int _i;
};



class CMatrix::GPair::Corrs::Corr::Iterator
{
   float& operator*();
   void operator++();
   bool operator!=(const Iterator&);
private:
   Iterator(Corr*,int);
   Corr* _p;
   int _i;
};



////////////////////////////////////////////////////////////////////////////////////////
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
