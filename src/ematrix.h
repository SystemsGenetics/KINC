#ifndef EMATRIX_H
#define EMATRIX_H
#include <ace.h>
#include <fstream>
#include <vector>



class EMatrix : public AccelCompEng::DataPlugin
{
public:
   ACE_DATA_HEADER()
   ACE_EXCEPTION(EMatrix,OutOfRange)
   ACE_EXCEPTION(EMatrix,AlreadyInitialized)
   ACE_EXCEPTION(EMatrix,NotInitialized)
   ACE_EXCEPTION(EMatrix,AlreadySet)
   ACE_EXCEPTION(EMatrix,InvalidArg)
   ACE_EXCEPTION(EMatrix,InvalidFile)
   ACE_EXCEPTION(EMatrix,NotNewFile)
   ACE_EXCEPTION(EMatrix,CannotOpen)
   enum Transform { none=0,log,log2,log10 };
   class Gene;
   class Mirror;
   using string = std::string;
   EMatrix(const string&,const string&);
   ~EMatrix();
   void load(GetOpts &ops, Terminal &tm) override final;
   void dump(GetOpts &ops, Terminal &tm) override final;
   void query(GetOpts &ops, Terminal &tm) override final;
   bool empty() override final;
   void initialize(int,int,bool);
   bool hasSampleHead() const;
   int gSize() const;
   int sSize() const;
   string& gName(int);
   string& sName(int);
   Transform transform() const;
   void transform(Transform);
   void write();
   Gene begin();
   Gene end();
   Gene find(int);
   Gene find(const string&);
   Gene& at(int);
   Gene& operator[](int);
private:
   using svec = std::vector<string>;
   ACE_FMEM_STATIC(Header,34)
      ACE_FMEM_VAL(sSize,uint32_t,0)
      ACE_FMEM_VAL(gSize,uint32_t,4)
      ACE_FMEM_VAL(tr,uint8_t,8)
      ACE_FMEM_VAL(wr,uint8_t,9)
      ACE_FMEM_VAL(sPtr,FPtr,10)
      ACE_FMEM_VAL(gPtr,FPtr,18)
      ACE_FMEM_VAL(eData,FPtr,26)
   ACE_FMEM_END()
   ACE_FMEM_STATIC(NmHead,8)
      ACE_FMEM_VAL(nPtr,FPtr,0)
   ACE_FMEM_END()
   void read_sizes(std::ifstream&,int);
   void read_header(std::ifstream&);
   void read_gene_expressions(std::ifstream&,const string&);
   void lookup(GetOpts&,Terminal&);
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
   using string = std::string;
   const string& name() const;
   void read();
   void write();
   Iterator begin();
   Iterator end();
   float& at(int);
   float& operator[](int);
   void operator++();
   bool operator!=(const Gene&);
   bool operator==(const Gene&);
private:
   Gene(EMatrix*,int);
   void set(int);
   static FPtr initialize(FileMem&,int,int);
   ACE_FMEM_OBJECT(Expr)
      Expr(int sSize, FPtr ptr = fNullPtr): Object(4*sSize,ptr) {}
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
   bool operator==(const Iterator&);
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
   Mirror(EMatrix&);
   ~Mirror();
   Mirror(const Mirror&) = delete;
   Mirror& operator=(const Mirror&) = delete;
   Mirror(Mirror&&) = delete;
   Mirror& operator=(Mirror&&) = delete;
   void read();
   void write();
   Gene begin();
   Gene end();
   Gene& at(int);
   Gene& operator[](int);
private:
   ACE_FMEM_OBJECT(AllExps)
      AllExps(int gSize, int sSize, FPtr ptr = fNullPtr):
         Object(4*gSize*sSize,ptr),
         _sSize(sSize)
      {}
      float& val(int gi, int si) { get<float>(4*((gi*_sSize)+si)); }
      int _sSize;
   ACE_FMEM_END()
   EMatrix* _p;
   AllExps _data;
   Gene* _iGene {nullptr};
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
   bool operator!=(const Gene&);
   bool operator==(const Gene&);
private:
   Gene(Mirror*,int);
   void set(int);
   Mirror* _p;
   int _i;
};



class EMatrix::Mirror::Gene::Iterator
{
public:
   friend class Gene;
   float& operator*();
   void operator++();
   bool operator!=(const Iterator&);
   bool operator==(const Iterator&);
private:
   Iterator(Gene*,int);
   Gene* _p;
   int _i;
};



#endif
