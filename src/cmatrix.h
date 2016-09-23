#ifndef CMATRIX_H
#define CMATRIX_H
#include <ace.h>
#include <vector>
#include "ematrix.h"



namespace Ace = AccelCompEng;



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
   void init() override final;
   void load(Ace::GetOpts&,Ace::Terminal&) override final;
   /// Function handler for user dump command.
   void dump(Ace::GetOpts&,Ace::Terminal&) override final;
   void query(Ace::GetOpts&,Ace::Terminal&) override final;
   bool empty() override final;
   static long long diagonal(int,int);
   static long long diag_size(int);
   void initialize(std::vector<std::string>&& geneNames, std::vector<std::string>&& sampleNames,
                   std::vector<std::string>&& correlationNames, int maxModes);
   int gene_size() const;
   int sample_size() const;
   int max_modes() const;
   int correlation_size() const;
   const std::string& gene_name(int) const;
   const std::string& sample_name(int) const;
   const std::string& correlation_name(int) const;
   GPair begin();
   GPair end();
   GPair find(int,int);
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



#endif
