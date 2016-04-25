#ifndef EMATRIX_H
#define EMATRIX_H
#include <fstream>
#include "../../../dataplugin.h"



namespace ematrixData {

struct Header : public FileMem::Static<33>
{
   using FPtr = FileMem::Ptr;
   using FileMem::Static<33>::Static;
   uint32_t& sampleSize() { get<uint32_t>(0); }
   uint32_t& geneSize() { get<uint32_t>(4); }
   uint8_t& transform() { get<uint8_t>(8); }
   FPtr& samplePtr() { get<FPtr>(9); }
   FPtr& genePtr() { get<FPtr>(17); }
   FPtr& expPtr() { get<FPtr>(25); }
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
   Expressions(int amt): Object(4*amt) {}
   float& val(int n) { get<float>(4*n); }
};

}



class ematrix : public DataPlugin
{
public:
   enum Transform { none=0,log,log2,log10 };
   using Hdr = ematrixData::Header;
   using gHdr = ematrixData::GeneHdr;
   using sHdr = ematrixData::SampleHdr;
   using Exp = ematrixData::Expression;
   using string = std::string;
   using ifile = std::ifstream;
   struct NotNewFile;
   struct CannotOpen;
   ematrix(const string& type, const string& file);
   void load(GetOpts &ops, Terminal &tm) override final;
   void dump(GetOpts &ops, Terminal &tm) override final {}
   void query(GetOpts &ops, Terminal &tm) override final {}
   bool empty() override final { return true; }
private:
   bool load_samples(ifile&);
   bool load_genes(ifile&);
   bool load_data(ifile&);
   bool load_blank();
   Hdr _hdr;
   FileMem& _mem;
};



struct ematrix::NotNewFile : public DataException
{
   NotNewFile(const char* file, int line):
      DataException(file,line,"ematrix::NotNewFile",Level::fatal)
   {}
};

struct ematrix::CannotOpen : public DataException
{
   CannotOpen(const char* file, int line):
      DataException(file,line,"ematrix::CannotOpen",Level::fatal)
   {}
};



#endif
