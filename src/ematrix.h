#ifndef EMATRIX_H
#define EMATRIX_H
#include <ace.h>
#include <fstream>
#include <vector>


/*
namespace ExprMatrixData {

using namespace AccelCompEng;


struct Header : public FileMem::Static<33>
{
   using FPtr = FileMem::Ptr;
   using FileMem::Static<33>::Static;
   uint32_t& sampleSize() { get<uint32_t>(0); }
   const uint32_t& sampleSize() const { get<uint32_t>(0); }
   uint32_t& geneSize() { get<uint32_t>(4); }
   const uint32_t& geneSize() const { get<uint32_t>(4); }
   uint8_t& transform() { get<uint8_t>(8); }
   FPtr& samplePtr() { get<FPtr>(9); }
   const FPtr& samplePtr() const { get<FPtr>(9); }
   FPtr& genePtr() { get<FPtr>(17); }
   const FPtr& genePtr() const { get<FPtr>(17); }
   FPtr& expPtr() { get<FPtr>(25); }
   const FPtr& expPtr() const { get<FPtr>(25); }
};

struct SampleHdr : public FileMem::Static<8>
{
   using FPtr = FileMem::Ptr;
   using FileMem::Static<8>::Static;
   FPtr& name() { get<FPtr>(0); }
};

struct GeneHdr : public FileMem::Static<8>
{
   using FPtr = FileMem::Ptr;
   using FileMem::Static<8>::Static;
   FPtr& name() { get<FPtr>(0); }
};

struct Expression : public FileMem::Static<4>
{
   using FileMem::Static<4>::Static;
   float& val() { get<float>(0); }
};

struct Expressions : public FileMem::Object
{
   using FPtr = FileMem::Ptr;
   Expressions(int amt, FPtr ptr = FileMem::nullPtr): Object(4*amt,ptr) {}
   float& val(int n) { get<float>(4*n); }
   const float& val(int n) const { get<float>(4*n); }
};

}



class EMatrix : public AccelCompEng::DataPlugin
{
public:
   ACE_DATA_HEADER()
   ACE_FMEM_STATIC(Hdr,33)
      using FPtr = FileMem::Ptr;
      ACE_FMEM_VAL(sampleSize,uint32_t,0)
      ACE_FMEM_VAL(geneSize,uint32_t,4)
      ACE_FMEM_VAL(transform,uint32_t,8)
      ACE_FMEM_VAL(samplePtr,FPtr,9)
      ACE_FMEM_VAL(genePtr,FPtr,17)
      ACE_FMEM_VAL(expPtr,FPtr,25)
   ACE_FMEM_END()
   // *
   // * ENUMERATIONS
   // *
   enum Transform { none=0,log,log2,log10 };
   // *
   // * DECLERATIONS
   // *
   //using Hdr = ExprMatrixData::Header;
   using gHdr = ExprMatrixData::GeneHdr;
   using sHdr = ExprMatrixData::SampleHdr;
   using Exp = ExprMatrixData::Expression;
   using Exps = ExprMatrixData::Expressions;
   using FileSync = AccelCompEng::FileSync;
   using KincFile = AccelCompEng::File;
   //using ace = AccelCompEng;
   using string = std::string;
   using ifile = std::ifstream;
   // *
   // * EXCEPTIONS
   // *
   ACE_EXCEPTION(EMatrix,NotNewFile)
   ACE_EXCEPTION(EMatrix,CannotOpen)
   ACE_EXCEPTION(EMatrix,InvalidFile)
   ACE_EXCEPTION(EMatrix,InvalidArg)
   ACE_EXCEPTION(EMatrix,OutOfRange)
   ACE_EXCEPTION(EMatrix,BufferNotLoaded)
   // *
   // * BASIC METHODS
   // *
   EMatrix(const string& type, const string& file);
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
   int sample_size() const;
   int gene_size() const;
   string sample_name(int) const;
   string gene_name(int) const;
   void load_buffer();
   void clear_buffer();
   const float* gene(int) const;
private:
   // *
   // * FUNCTIONS
   // *
   void lookup(GetOpts&,Terminal&);
   void load_samples(Terminal&,ifile&);
   void load_genes(Terminal&,ifile&,string,Transform);
   // *
   // * VARIABLES
   // *
   Hdr _hdr;
   FileMem& _mem;
   Exps* _data;
};
*/


class EMatrix : public AccelCompEng::DataPlugin
{
public:
   ACE_DATA_HEADER()
   ACE_EXCEPTION(EMatrix,OutOfRange)
   ACE_EXCEPTION(EMatrix,NoData)
   ACE_EXCEPTION(EMatrix,AlreadyInit)
   ACE_EXCEPTION(EMatrix,NotInit)
   ACE_EXCEPTION(EMatrix,AlreadySet)
   enum class Which {sample,gene};
   class Gene;
   class Mirror;
   using string = std::string;
   EMatrix(const string&,const string&);
   ~EMatrix();
   void load(GetOpts &ops, Terminal &tm) override final {}
   void dump(GetOpts &ops, Terminal &tm) override final {}
   void query(GetOpts &ops, Terminal &tm) override final {}
   bool empty() override final {}
   void initialize(int,int,bool);
   bool sHead() const;
   int gSize() const;
   int sSize() const;
   const string& gName(int) const;
   void gName(int,const string&);
   const string& sName(int) const;
   void sName(int,const string&);
   void write();
   Gene begin();
   Gene end();
   Gene& at(int);
   Gene& operator[](int);
private:
   enum Transform { none=0,log,log2,log10 };
   using svec = std::vector<string>;
   ACE_FMEM_STATIC(Header,33)
      ACE_FMEM_VAL(sSize,uint32_t,0)
      ACE_FMEM_VAL(gSize,uint32_t,4)
      ACE_FMEM_VAL(tr,uint32_t,8)
      ACE_FMEM_VAL(sPtr,FPtr,9)
      ACE_FMEM_VAL(gPtr,FPtr,17)
      ACE_FMEM_VAL(eData,FPtr,25)
   ACE_FMEM_END()
   ACE_FMEM_STATIC(NmHead,8)
      ACE_FMEM_VAL(nPtr,FPtr,0)
   ACE_FMEM_END()
   Header _hdr {fNullPtr};
   Gene* _iGene {nullptr};
   svec _gNames;
   svec _sNames;
};



class EMatrix::Gene
{
public:
   friend class EMatrix;
   class Iterator;
   void read();
   void write();
   Iterator begin();
   Iterator end();
   float& at(int);
   float& operator[](int);
   void operator++();
   bool operator!=(const Gene&);
private:
   Gene(EMatrix*,int);
   void set(int);
   ACE_FMEM_OBJECT(Expr)
      Expr(int cSize, FPtr ptr = fNullPtr): Object(4*cSize,ptr) {}
      float& val(int i) { get<float>(4*i); }
   ACE_FMEM_END()
   EMatrix* _p;
   Expr _head;
   int _i;
};



class EMatrix::Gene::Iterator
{
public:
   friend class Gene;
   void operator++();
   float& operator*();
   bool operator!=(const Iterator&);
private:
   Iterator(Gene*,int);
   Gene* _p;
   int _i;
};



class EMatrix::Mirror
{
public:
   friend class EMatrix;
   class Gene;
   void read();
   void write();
   Gene begin();
   Gene end();
   Gene& at(int);
   Gene& operator[](int);
private:
   Mirror(EMatrix*);
   ACE_FMEM_OBJECT(AllExps)
      AllExps(int gSize, int cSize, FPtr ptr = fNullPtr):
         Object(4*gSize*cSize,ptr),
         _cSize(cSize)
      {}
      float& val(int gi, int ci) { get<float>((gi*_cSize)+ci); }
      int _cSize;
   ACE_FMEM_END()
   AllExps _data;
   Gene* _iGene;
};



class EMatrix::Mirror::Gene
{
public:
   friend class Mirror;
   class Iterator;
   Iterator begin();
   Iterator end();
   float& at(int);
   float& operator[](int);
   void operator++();
   Gene& operator!=(const Gene&);
private:
   Gene(Mirror*,int);
   Mirror* _p;
   int _i;
};



class EMatrix::Mirror::Gene::Iterator
{
public:
   friend class Array;
   float& operator[](int);
   void operator++();
   float& operator*();
   Iterator& operator!=(const Iterator&);
private:
   Iterator(Gene*,int);
   Gene* _p;
   int _i;
};



#endif
