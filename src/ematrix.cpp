#include "ematrix.h"
#include <sstream>
#include <cmath>



EMatrix::EMatrix(const string& type, const string& file):
   DataPlugin(type,file),
   _mem(*KincFile::mem()),
   _data(nullptr)
{
   // Does the KINC file have a header already set. This case would be for
   // a new file.
   if (KincFile::head()==FileMem::nullPtr)
   {
      // Allocate space in the file (file memory) for the header.
      _mem.allot(_hdr);
      AccelCompEng::File::head(_hdr.addr());

      // Set defaults for the ematrix header.
      _hdr.sampleSize() = 0;
      _hdr.geneSize() = 0;
      _hdr.transform() = Transform::none;
      _hdr.samplePtr() = FileMem::nullPtr;
      _hdr.genePtr() = FileMem::nullPtr;
      _hdr.expPtr() = FileMem::nullPtr;

      // Write the header initially to file.
      _mem.sync(_hdr,FileSync::write);
   }
   // The KINC file already has a header if the file already exists and is
   // being opened.
   else
   {
      _hdr = KincFile::head();
      _mem.sync(_hdr,FileSync::read);
   }
}


/// Imports data from an expression matrix text file.
///
///
///
/// @param ops
///   The command-line arguments provided by the user.
/// @param tm
///
void EMatrix::load(GetOpts &ops, Terminal &tm)
{
   bool hasHeaders {true};
   int sampleSize {0};
   Transform tr {none};
   for (auto i = ops.begin();i!=ops.end();++i)
   {
      if (i.is_key("noheader"))
      {
         hasHeaders = false;
      }
      else if (i.is_key("samples"))
      {
         sampleSize = i.value<int>();
      }
      else if (i.is_key("transform"))
      {
         string type = i.value<string>();
         if (type==string("log"))
         {
            tr = log;
         }
         else if (type==string("log2"))
         {
            tr = log2;
         }
         else if (type==string("log10"))
         {
            tr = log10;
         }
         else
         {
            throw InvalidArg(__LINE__);
         }
      }
   }
   if (!hasHeaders&&sampleSize==0)
   {
      tm << "If noheader option is given number of samples must be given.\n";
      throw InvalidArg(__LINE__);
   }
   AccelCompEng::assert<NotNewFile>(head()==FileMem::nullPtr,__LINE__);
   ifile f(ops.com_front());
   AccelCompEng::assert<CannotOpen>(f.is_open(),__LINE__);
   _mem.allot(_hdr);
   try
   {
      if (hasHeaders)
      {
         load_samples(tm,f);
      }
      else
      {
         _hdr.sampleSize() = sampleSize;
         _hdr.samplePtr() = FileMem::nullPtr;
      }
      load_genes(tm,f,"NaN",tr);
   }
   catch (...)
   {
      KincFile::clear();
      throw;
   }
   KincFile::head(_hdr.addr());
   _hdr.transform() = tr;
   _mem.sync(_hdr,FileSync::write);
   tm << "Loaded " << _hdr.geneSize() << " gene(s) with " << _hdr.sampleSize()
      << " sample(s) each.\n";
}



void EMatrix::dump(GetOpts &ops, Terminal &tm)
{
   tm << "dump command not implemented for EMatrix type.\n";
}



void EMatrix::query(GetOpts &ops, Terminal &tm)
{
   enum {Basic=0,Lookup,Dump};
   switch (ops.com_get({"lookup","dump"}))
   {
   case Basic:
      if (empty())
      {
         tm << "No data loaded.\n";
      }
      else
      {
         tm << _hdr.geneSize() << " gene(s) and " << _hdr.sampleSize()
            << " sample(s).\n";
         tm << "Sample Transformation: ";
         switch (_hdr.transform())
         {
         case none:
            tm << "none.\n";
            break;
         case log:
            tm << "Natural Logarithm.\n";
            break;
         case log2:
            tm << "Logarithm Base 2.\n";
            break;
         case log10:
            tm << "Logarithm Base 10.\n";
            break;
         }
         if (_hdr.samplePtr()==FileMem::nullPtr)
         {
            tm << "Sample header data does not exist.\n";
         }
         else
         {
            tm << "Sample header data exists.\n";
         }
      }
      break;
   case Lookup:
      ops.com_pop();
      lookup(ops,tm);
      break;
   case Dump:
      break;
   }
}



int EMatrix::sample_size() const
{
   return _hdr.sampleSize();
}



int EMatrix::gene_size() const
{
   return _hdr.geneSize();
}



EMatrix::string EMatrix::sample_name(int n) const
{
   bool cond {_hdr.samplePtr()!=FileMem::nullPtr&&n>=0||n<_hdr.sampleSize()};
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
   sHdr shdr(_hdr.samplePtr());
   _mem.sync(shdr,FileSync::read,n);
   FString name(&_mem,shdr.name());
   return *name;
}



EMatrix::string EMatrix::gene_name(int n) const
{
   bool cond {n>=0||n<_hdr.geneSize()};
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
   gHdr ghdr(_hdr.genePtr());
   _mem.sync(ghdr,FileSync::read,n);
   FString name(&_mem,ghdr.name());
   return *name;
}



void EMatrix::load_buffer()
{
   if (!_data)
   {
      _data = new Exps(_hdr.geneSize()*_hdr.sampleSize(),_hdr.expPtr());
      _mem.sync(*_data,FileSync::read);
   }
}



void EMatrix::clear_buffer()
{
   if (_data)
   {
      delete _data;
      _data = nullptr;
   }
}



const float* EMatrix::gene(int n) const
{
   AccelCompEng::assert<BufferNotLoaded>(_data,__LINE__);
   AccelCompEng::assert<OutOfRange>(n<_hdr.geneSize(),__LINE__);
   return &(_data->val(n*_hdr.sampleSize()));
}



void EMatrix::lookup(GetOpts &ops, Terminal &tm)
{
   int i;
   try
   {
      i = {std::stoi(ops.com_front())};
   }
   catch (std::exception)
   {
      i = 0;
      gHdr ghdr {_hdr.genePtr()};
      while (i<_hdr.geneSize())
      {
         _mem.sync(ghdr,FileSync::read,i++);
         FString name {&_mem,ghdr.name()};
         if (*name==ops.com_front())
         {
            --i;
            break;
         }
      }
   }
   if (i<0||i>=_hdr.geneSize())
   {
      tm << "Invalid gene index.\n";
      return;
   }
   gHdr ghdr {_hdr.genePtr()};
   _mem.sync(ghdr,FileSync::read,i);
   FString name {&_mem,ghdr.name()};
   tm << i << ": " << *name << " : ";
   Exps exps(_hdr.sampleSize(),_hdr.expPtr());
   _mem.sync(exps,FileSync::read,i);
   int x = 0;
   for (;x<(_hdr.sampleSize()-1);++x)
   {
      tm << exps.val(x) << ", ";
   }
   tm << exps.val(x) << "\n";
}



bool EMatrix::empty()
{
   return _hdr.geneSize()==0;
}



void EMatrix::load_samples(Terminal& tm, ifile& f)
{
   while (f.peek()==' '||f.peek()=='\t'||f.peek()=='\n')
   {
      AccelCompEng::assert<InvalidFile>(f,__LINE__);
      f.get();
   }
   tm << "Loading sample headers...\n";
   string buf;
   std::getline(f,buf);
   std::istringstream ibuf(buf);
   string tmp;
   std::vector<string> vs;
   while (ibuf >> tmp)
   {
      vs.push_back(std::move(tmp));
   }
   AccelCompEng::assert<InvalidFile>(vs.size()>0,__LINE__);
   tm << "Writing sample headers...\n";
   sHdr shdr;
   _mem.allot(shdr,vs.size());
   _hdr.sampleSize() = vs.size();
   _hdr.samplePtr() = shdr.addr();
   for (int i = 0;i<vs.size();++i)
   {
      FString tname(&_mem);
      tname = vs[i];
      shdr.name() = tname.addr();
      _mem.sync(shdr,FileSync::write,i);
   }
}



void EMatrix::load_genes(Terminal& tm, ifile& f, string nan, Transform t)
{
   tm << "Loading genes and expression data...\n";
   std::vector<string> gns;
   std::vector<std::vector<float>> exps;
   while (f)
   {
      string buf;
      std::getline(f,buf);
      if (!buf.empty())
      {
         std::istringstream ibuf(buf);
         string tmp;
         bool cond = ibuf >> tmp;
         AccelCompEng::assert<InvalidFile>(cond,__LINE__);
         gns.push_back(std::move(tmp));
         exps.push_back({});
         int count = 0;
         while (ibuf >> tmp)
         {
            if (tmp==nan)
            {
               exps.back().push_back(NAN);
            }
            else
            {
               try
               {
                  switch (t)
                  {
                  case none:
                     exps.back().push_back(std::stof(tmp));
                     break;
                  case log:
                     exps.back().push_back(logf(std::stof(tmp)));
                     break;
                  case log2:
                     exps.back().push_back(log2f(std::stof(tmp)));
                     break;
                  case log10:
                     exps.back().push_back(log10f(std::stof(tmp)));
                     break;
                  }
               }
               catch (std::exception)
               {
                  throw InvalidFile(__LINE__);
               }
            }
            ++count;
         }
         AccelCompEng::assert<InvalidFile>(count==_hdr.sampleSize(),__LINE__);
      }
   }
   AccelCompEng::assert<InvalidFile>(gns.size()>0,__LINE__);
   tm << "Writing gene headers...\n";
   gHdr ghdr;
   _mem.allot(ghdr,gns.size());
   _hdr.geneSize() = gns.size();
   _hdr.genePtr() = ghdr.addr();
   for (int i = 0;i<gns.size();++i)
   {
      FString tname(&_mem);
      tname = gns[i];
      ghdr.name() = tname.addr();
      _mem.sync(ghdr,FileSync::write,i);
   }
   tm << "Writing expression data...\n";
   Exp exp;
   _mem.allot(exp,exps.size()*exps.back().size());
   _hdr.expPtr() = exp.addr();
   for (int x = 0;x<exps.size();++x)
   {
      for (int y = 0;y<exps[x].size();++y)
      {
         exp.val() = exps[x][y];
         _mem.sync(exp,FileSync::write,x*_hdr.sampleSize()+y);
      }
   }
}
