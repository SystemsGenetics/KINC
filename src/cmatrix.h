#ifndef CMATRIX_H
#define CMATRIX_H
#include <ace.h>
#include <vector>


/*
class CMatrix : public AccelCompEng::DataPlugin
{
public:
   using string = std::string;
   ACE_DATA_HEADER()
   ACE_EXCEPTION(CMatrix,AlreadyInitialized)
   ACE_EXCEPTION(CMatrix,OutOfRange)
   ACE_EXCEPTION(CMatrix,NotInitialized)
   ACE_EXCEPTION(CMatrix,GreaterThanMax)
   ACE_EXCEPTION(CMatrix,AlreadySet)
   class GPair;
   CMatrix(const string& type, const string& file);
   void load(GetOpts&,Terminal&) override final;
   void dump(GetOpts&,Terminal&) override final;
   void query(GetOpts&,Terminal&) override final;
   bool empty() override final;
   static long long diagonal(int,int);
   static long long diag_size(int);
   void initialize(int,int,int,int,bool);
   int gSize() const;
   int sSize() const;
   int mSize() const;
   int cSize() const;
   string& gName(int);
   string& sName(int);
   string& cName(int);
   void write();
   GPair begin();
   GPair end();
   GPair& at(int,int);
   GPair& ref(int,int);
private:
   using svec = std::vector<string>;
   ACE_FMEM_STATIC(Header,46)
      ACE_FMEM_VAL(gSize,int32_t,0)
      ACE_FMEM_VAL(sSize,int32_t,4)
      ACE_FMEM_VAL(cSize,int32_t,8)
      ACE_FMEM_VAL(mSize,int8_t,12)
      ACE_FMEM_VAL(wr,int8_t,13)
      ACE_FMEM_VAL(gPtr,FPtr,14)
      ACE_FMEM_VAL(sPtr,FPtr,22)
      ACE_FMEM_VAL(cPtr,FPtr,30)
      ACE_FMEM_VAL(eData,FPtr,38)
   ACE_FMEM_END()
   ACE_FMEM_STATIC(NmHead,8)
      ACE_FMEM_VAL(nPtr,FPtr,0)
   ACE_FMEM_END()
   string& get_name(int,svec&,FPtr,int);
   void write_names(svec&,FPtr);
   Header _hdr;
   GPair* _iGPair {nullptr};
   svec _gNames;
   svec _sNames;
   svec _cNames;
};



class CMatrix::GPair
{
public:
   friend class CMatrix;
   class Modes;
   class Corrs;
   ~GPair();
   void read();
   void write();
   int x() const;
   int y() const;
   int size() const;
   void size(int);
   Modes& modes();
   Corrs& corrs();
   void operator++();
   bool operator!=(const GPair&);
private:
   GPair(CMatrix*,int,int);
   void set(int,int);
   static FPtr initialize(FileMem&,int,int,int,int);
   ACE_FMEM_OBJECT(Pair)
      Pair(int mSize, int sSize, int cSize, FPtr ptr = fNullPtr):
         Object(1+(mSize*sSize)+(4*cSize*mSize),ptr),
         _sSize(sSize),
         _cSize(cSize),
         _indent(1+(mSize*sSize))
      {}
      ACE_FMEM_VAL(mAmt,int8_t,0)
      int8_t& mVal(int mi, int i) { get<int8_t>(1+(_sSize*mi)+i); }
      float& cVal(int mi, int ci) { get<float>(_indent+(4*((_cSize*mi)+ci))); }
      int _sSize;
      int _cSize;
      int _indent;
   ACE_FMEM_END()
   CMatrix* _p;
   Pair _data;
   int _x;
   int _y;
   Modes* _iModes {nullptr};
   Corrs* _iCorrs {nullptr};
};



class CMatrix::GPair::Modes
{
public:
   friend class GPair;
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



class CMatrix::GPair::Modes::Mode
{
public:
   friend class Modes;
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



class CMatrix::GPair::Modes::Mode::Iterator
{
public:
   friend class Mode;
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
public:
   friend class GPair;
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
public:
   friend class Corrs;
   class Iterator;
   Iterator begin();
   Iterator end();
   float& at(int);
   float& operator[](int);
   void operator++();
   bool operator!=(const Corr&);
private:
   Corr(Corrs*,int);
   void set(int);
   Corrs* _p;
   int _i;
};



class CMatrix::GPair::Corrs::Corr::Iterator
{
public:
   friend class Corr;
   float& operator*();
   void operator++();
   bool operator!=(const Iterator&);
private:
   Iterator(Corr*,int);
   Corr* _p;
   int _i;
};

*/

#endif
