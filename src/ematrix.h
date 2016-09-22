#ifndef EMATRIX_H
#define EMATRIX_H
#include <ace.h>
#include <fstream>
#include <vector>



namespace Ace = AccelCompEng;



class EMatrix : public Ace::Data, private Ace::NVMemory::Node
{
public:
   struct OutOfRange : public Ace::Exception { using Ace::Exception::Exception; };
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
   EMatrix();
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
   Gene find(int);
   Gene find(const std::string& geneName);
private:
   struct __attribute__ ((__packed__)) Header
   {
      int32_t _sampleSize {0};
      int32_t _geneSize {0};
      uint8_t _transform {none};
      int64_t _samplePtr {fnullptr};
      int64_t _genePtr {fnullptr};
      int64_t _expData {fnullptr};
   };
   Header& data() { return get<Header>(); }
   const Header& data() const { return get<Header>(); }
   void read_headers(std::ifstream& file, int sampleSize, Transform transform, bool hasHeaders);
   void read_gene_expressions(std::ifstream& file, Ace::Terminal& tm, const std::string& nanStr);
   void lookup(Ace::GetOpts&,Ace::Terminal&);
   void skip_blanks(std::ifstream& file);
   bool is_blank_line(const std::string& line);
   void null_data() override final;
   void flip_endian() override final;
   Header _header;
   bool _isNew {true};
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
   Mirror(const Mirror&) = delete;
   Mirror& operator=(const Mirror&) = delete;
   Mirror(Mirror&&) = delete;
   Mirror& operator=(Mirror&&) = delete;
   void read();
   void write();
   float& value(int gene, int i);
   Gene begin();
   Gene end();
   Gene find(int);
private:
   void null_data() override final {}
   void flip_endian() override final;
   EMatrix* _p;
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
