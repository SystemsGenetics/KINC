#ifndef CMATRIX_H
#define CMATRIX_H
#include <ace.h>
#include <vector>
#include "ematrix.h"



namespace Ace = AccelCompEng;



/// Data plugin class that holds correlation data for gene statistics.
class CMatrix : public Ace::Data, private Ace::NVMemory::Node
{
public:
   struct NullMatrix : public Ace::Exception { using Ace::Exception::Exception; };
   struct OutOfRange : public Ace::Exception { using Ace::Exception::Exception; };
   struct NotInitialized : public Ace::Exception { using Ace::Exception::Exception; };
   struct GreaterThanMax : public Ace::Exception { using Ace::Exception::Exception; };
   struct AlreadySet : public Ace::Exception { using Ace::Exception::Exception; };
   struct InvalidSize : public Ace::Exception { using Ace::Exception::Exception; };
   struct AllocateFail : public Ace::Exception { using Ace::Exception::Exception; };
   class GPair;
   CMatrix();
   ~CMatrix();
   /// ACE initialize call.
   void init() override final;
   /// ACE load call.
   void load(Ace::GetOpts&,Ace::Terminal&) override final;
   /// ACE dump call.
   void dump(Ace::GetOpts&,Ace::Terminal&) override final;
   /// ACE query call.
   void query(Ace::GetOpts&,Ace::Terminal&) override final;
   /// ACE empty call.
   bool empty() override final;
   /// Get one dimensional position of diagonal gene correlation.
   /// @param x Index of first gene.
   /// @param y Index of second gene.
   /// @return One dimensional position of gene correlation.
   static long long diagonal(int x, int y);
   /// Get total size of diagonal matrix.
   /// @param x Total number of genes.
   /// @return Total number of correlations between all genes.
   static long long diag_size(int x);
   /// Initialize a new correlation matrix for this object.
   /// @param geneNames List of all gene names.
   /// @param sampleNames List of all sample names.
   /// @param correlationNames List of all correlation names.
   /// @param maxModes Maximum number of modes each correlation can hold.
   void initialize(std::vector<std::string>&& geneNames, std::vector<std::string>&& sampleNames,
                   std::vector<std::string>&& correlationNames, int maxModes);
   /// Get total number of genes.
   int gene_size() const;
   /// Get total number of samples per gene.
   int sample_size() const;
   /// Get maximum allowed number of modes per gene.
   int max_modes() const;
   /// Get number of correlations per gene.
   int correlation_size() const;
   /// Get name of specific gene.
   /// @param i Index of gene to lookup.
   /// Name of gene.
   const std::string& gene_name(int i) const;
   /// Get name of specific sample column.
   /// @param i Index of sample column to lookup.
   /// Name of sample column.
   const std::string& sample_name(int i) const;
   /// Get name of correlation column.
   /// @param i Index of correlation column.
   /// Name of correlation column.
   const std::string& correlation_name(int i) const;
   /// First gene pair iterator.
   GPair begin();
   /// End of list gene pair iterator.
   GPair end();
   /// Find gene pair iterator.
   /// @param x Index of first gene.
   /// @param y Index of second gene.
   GPair find(int x, int y);
private:
   struct __attribute__ ((__packed__)) Header
   {
      int32_t _geneSize;
      int32_t _sampleSize;
      int32_t _correlationSize;
      int8_t _maxModes;
      int64_t _genePtr;
      int64_t _samplePtr;
      int64_t _correlationPtr;
      int64_t _expData;
   };
   Header& data() { return get<Header>(); }
   const Header& data() const { return get<Header>(); }
   void null_data() override final;
   void flip_endian() override final;
   bool _isNew {true};
   GPair* _iGPair {nullptr};
   std::vector<std::string> _geneNames;
   std::vector<std::string> _sampleNames;
   std::vector<std::string> _correlationNames;
   constexpr static int _strSize {1024};
};



/// Iterator that holds correlation data for a single gene pair.
class CMatrix::GPair: private Ace::NVMemory::Node
{
public:
   friend class CMatrix;
   class Modes
   {
   public:
      friend class GPair;
      class Mode;
      Modes(const Modes&);
      Modes(Modes&&) = default;
      Modes& operator=(const Modes&);
      Modes& operator=(Modes&&) = default;
      Mode begin();
      Mode end();
      Mode find(int);
   private:
      Modes(GPair*);
      GPair* _p;
      std::unique_ptr<Mode> _mode {nullptr};
   };
   class Corrs
   {
   public:
      friend class GPair;
      class Corr;
      Corrs(const Corrs&);
      Corrs(Corrs&&) = default;
      Corrs& operator=(const Corrs&);
      Corrs& operator=(Corrs&&) = default;
      Corr begin();
      Corr end();
      Corr find(int);
   private:
      Corrs(GPair*);
      GPair* _p;
      std::unique_ptr<Corr> _corr {nullptr};
   };
   GPair(const GPair&) = default;
   GPair(GPair&&) = default;
   GPair& operator=(const GPair&) = default;
   GPair& operator=(GPair&&) = default;
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
   int8_t& mode_val(int mi, int i);
   float& correlation_val(int mi, int ci);
   static int64_t initialize(Ace::NVMemory& mem, int geneSize, int maxModes, int sampleSize,
                             int correlationSize);
   static size_t calc_size(int maxModes, int sampleSize, int correlationSize);
   void null_data() override final {}
   void flip_endian() override final;
   CMatrix* _p;
   int _x;
   int _y;
   Modes _modes;
   Corrs _corrs;
};



/// Iterator that holds a single mode mask for a single mode.
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



/// Iterator for a single bit of a mode mask.
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



/// Iterator for a single correlation of a single mode of a gene pair.
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



/// Iterator for a single correlation value of a single correlation of a single mode of a single
/// gene pair.
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



class CMatrixHelpItem : public AbstractHelpItem
{
public:
   std::string getName() const override final
   {
      return "cmx";
   }
   std::string getDescription() const override final
   {
      return "cmx Data Object\n"
            "\n"
            "This data object holds correlation matrix data.\n"
             "\n"
             "dump [file] (--threshold)\n"
             "dumps gene correlation data that is above a given threshold to flat file [file]."
             " Default threshold is 99.0 unless otherwise specified by --threshold.\n"
             "\n"
             "load command is not implemented\n"
             "\n"
             "query outputs basic information.\n";
   }
};


#endif
