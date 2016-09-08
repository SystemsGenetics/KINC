#ifndef EMATRIX_H
#define EMATRIX_H
#include <ace.h>
#include <fstream>
#include <vector>



namespace Ace = AccelCompEng;



class EMatrix : public Ace::Data
{
public:
   struct OutOfRange : public Ace::Exception { using Ace::Exception::Exception; };
   struct NotInitialized : public Ace::Exception { using Ace::Exception::Exception; };
   struct AlreadySet : public Ace::Exception { using Ace::Exception::Exception; };
   struct InvalidArg : public Ace::Exception { using Ace::Exception::Exception; };
   struct InvalidFile : public Ace::Exception { using Ace::Exception::Exception; };
   struct NotNewFile : public Ace::Exception { using Ace::Exception::Exception; };
   struct CannotOpen : public Ace::Exception { using Ace::Exception::Exception; };
   struct NullMatrix : public Ace::Exception { using Ace::Exception::Exception; };
   struct InvalidSize : public Ace::Exception { using Ace::Exception::Exception; };
   enum Transform { none=0,log,log2,log10 };
   class Gene;
   class Mirror;
   EMatrix() = default;
   ~EMatrix();
   void init() override final;
   void load(Ace::GetOpts &ops, Ace::Terminal &tm) override final;
   void dump(Ace::GetOpts &ops, Ace::Terminal &tm) override final;
   void query(Ace::GetOpts &ops, Ace::Terminal &tm) override final;
   bool empty() override final;
   void initialize(std::vector<std::string>&& geneNames,std::vector<std::string>&& sampleNames,
                   Transform transform);
   int gene_size() const;
   int sample_size() const;
   const std::string& gene_name(int i) const;
   const std::string& sample_name(int i) const;
   Transform transform() const;
   Gene begin();
   Gene end();
   Gene find(const std::string& geneName);
   Gene& at(int i);
   Gene& operator[](int);
private:
   struct Header : public Ace::NVMemory::Node
   {
      Header():
         Node(sizeof(HData))
      {
         init_data<HData>();
      }
      Header(const std::shared_ptr<Ace::NVMemory>& mem):
         Node(sizeof(HData),mem)
      {
         init_data<HData>();
      }
      Header(const std::shared_ptr<Ace::NVMemory>& mem, int64_t ptr):
         Node(sizeof(HData),mem,ptr)
      {
         init_data<HData>();
      }
      struct __attribute__ ((__packed__)) HData
      {
         int32_t _sampleSize {0};
         int32_t _geneSize {0};
         uint8_t _transform {none};
         int64_t _samplePtr {fnullptr};
         int64_t _genePtr {fnullptr};
         int64_t _expData {fnullptr};
      };
      HData& data() { return get<HData>(); }
      const HData& data() const { return get<HData>(); }
      void null_data()
      {
         get<HData>()._sampleSize = 0;
         get<HData>()._geneSize = 0;
         get<HData>()._transform = none;
         get<HData>()._samplePtr = fnullptr;
         get<HData>()._genePtr = fnullptr;
         get<HData>()._expData = fnullptr;
      }
      void flip_endian()
      {
         flip(0,4);
         flip(4,4);
         flip(9,8);
         flip(17,8);
         flip(25,8);
      }
   };
   void read_sizes(std::ifstream& file, int& geneSize, int& sampleSize);
   //void read_header(std::ifstream&);
   //void read_gene_expressions(std::ifstream&,const string&);
   //void lookup(Ace::GetOpts&,Ace::Terminal&);
   Header _header;
   bool _isNew {true};
   Gene* _iGene {nullptr};
   std::vector<std::string> _geneNames;
   std::vector<std::string> _sampleNames;
   constexpr static int _strSize {1024};
};



class EMatrix::Gene : private Ace::NVMemory::Node
{
public:
   friend class EMatrix;
   class Iterator;
   const std::string& name() const;
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
   void null_data() override final {}
   void flip_endian() override final;
   static int64_t initialize(Ace::NVMemory&,int,int);
   EMatrix* _p;
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



class EMatrix::Mirror : private Ace::NVMemory::Node
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
   float& value(int gene, int i);
   Gene begin();
   Gene end();
   Gene& at(int);
   Gene& operator[](int);
private:
   EMatrix* _p;
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
