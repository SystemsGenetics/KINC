#ifndef EMATRIX_H
#define EMATRIX_H
#include <ace.h>
#include <fstream>
#include <vector>



namespace Ace = AccelCompEng;



/// Main class for expression matrix.
class EMatrix : public Ace::Data, private Ace::NVMemory::Node
{
public:
   struct OutOfRange : public Ace::Exception { using Ace::Exception::Exception; };
   struct AlreadySet : public Ace::Exception { using Ace::Exception::Exception; };
   struct InvalidArg : public Ace::Exception { using Ace::Exception::Exception; };
   struct NotNewFile : public Ace::Exception { using Ace::Exception::Exception; };
   struct CannotOpen : public Ace::Exception { using Ace::Exception::Exception; };
   struct NullMatrix : public Ace::Exception { using Ace::Exception::Exception; };
   struct InvalidSize : public Ace::Exception { using Ace::Exception::Exception; };
   enum Transform { none=0,log,log2,log10 };
   class Gene;
   class Mirror;
   EMatrix();
   /// ACE initialize call.
   void init() override final;
   /// ACE load call.
   void load(Ace::GetOpts &ops, Ace::Terminal &tm) override final;
   /// ACE dump call.
   void dump(Ace::GetOpts &ops, Ace::Terminal &tm) override final;
   /// ACE query call.
   void query(Ace::GetOpts &ops, Ace::Terminal &tm) override final;
   /// ACE empty call.
   bool empty() override final;
   /// Initialize a new expression matrix for this object.
   /// @param geneNames List of gene names.
   /// @param sampleNames List of samples names, if there are no sample names provide a list of
   /// empty strings that totals the number of samples per gene.
   /// @param transform The transform, if any, to be done on samples.
   void initialize(std::vector<std::string>&& geneNames,std::vector<std::string>&& sampleNames,
                   Transform transform);
   /// Get number of genes in expression matrix.
   int gene_size() const;
   /// Get number of samples per gene in expression matrix.
   int sample_size() const;
   /// Get gene name.
   /// @param i The increment in the list of gene names to get the name from.
   /// @return Name of gene.
   const std::string& gene_name(int i) const;
   /// Get sample name.
   /// @param i The increment in the list of sample names to get the name from.
   /// @return Name of sample.
   const std::string& sample_name(int i) const;
   /// Get the transform of this expression matrix.
   Transform transform() const;
   /// Get first gene iterator.
   Gene begin();
   /// Get end of list iterator for genes.
   Gene end();
   /// Find gene iterator.
   /// @param i Increment of gene to find within gene list.
   /// @return Gene iterator if found, else end of list iterator if not found.
   Gene find(int i);
   /// Find gene iterator.
   /// @param name Name of gene to find.
   /// @return Gene iterator if found, else end of list iterator if not found.
   Gene find(const std::string& name);
private:
   struct InvalidFileGOverflow {};
   struct InvalidFileGUnderflow {};
   struct InvalidFileSUnderflow {};
   struct InvalidFileSInvalid {};
   struct InvalidFileBlank {};
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



/// Iterator class for expression matrix that represents a single gene, containing all the gene's
/// samples.
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



/// Iterator class for EMatrix::Gene that represents a single sample of a gene.
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



/// Maps entire expression matrix into system memory for quick access. A mirror of a file expression
/// matrix stored in system memory.
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



/// Iterator class for EMatrix::Mirror that represents a single gene, containing all the gene's
/// samples.
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



/// Iterator class for EMatrix::Mirror::Gene that represents a single sample of a gene.
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




class EMatrixHelpItem : public AbstractHelpItem
{
public:
   std::string getName() const override final
   {
      return "emx";
   }
   std::string getDescription() const override final
   {
      return "emx Data Object\n"
            "\n"
            "This data object holds expression matrix data.\n"
             "\n"
             "dump command is not implemented.\n"
             "\n"
             "load (--noheader|--samples|--transform|--missing)\n"
             "--noheader specifies the flat file has no sample names. --samples=n states the flat"
             " file has n samples per gene. --transform=log|log2|log10 will cause all sample values"
             " to be transformed by log, log2, or log10. --missing=T states the flat file uses the"
             " string token T for missing sample values.\n"
             "\n"
             "query outputs basic information.\n";
   }
};


#endif
