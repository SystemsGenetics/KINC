#include "ematrix.h"
#include <sstream>
#include <cmath>



EMatrix::EMatrix():
   Node(sizeof(Header))
{
   init_data<Header>();
}



void EMatrix::init()
{
   try
   {
      Node::mem(File::mem());
      if (File::head()==fnullptr)
      {
         allocate();
         File::head(addr());
         null_data();
         write();
      }
      else
      {
         addr(File::head());
         read();
         Ace::FString fstr(File::mem(),data()._genePtr);
         fstr.static_buffer(_strSize);
         _geneNames.push_back(fstr.str());
         for (int i=1;i<data()._geneSize;++i)
         {
            fstr.bump();
            _geneNames.push_back(fstr.str());
         }
         fstr.load(data()._samplePtr);
         _sampleNames.push_back(fstr.str());
         for (int i=1;i<data()._sampleSize;++i)
         {
            fstr.bump();
            _sampleNames.push_back(fstr.str());
         }
         _isNew = false;
      }
   }
   catch (...)
   {
      File::clear();
      _isNew = true;
      throw;
   }
}



void EMatrix::load(Ace::GetOpts &ops, Ace::Terminal &tm)
{
   static const char* f = __PRETTY_FUNCTION__;
   bool hasHeaders {true};
   int sampleSize {0};
   std::string nan {"NaN"};
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
         std::string type = i.value<std::string>();
         if (type==std::string("log"))
         {
            tr = log;
         }
         else if (type==std::string("log2"))
         {
            tr = log2;
         }
         else if (type==std::string("log10"))
         {
            tr = log10;
         }
         else
         {
            Ace::assert<InvalidArg>(false,f,__LINE__);
         }
      } else if (i.is_key("nosample"))
      {
         nan = i.value<std::string>();
      }
   }
   if (!hasHeaders&&sampleSize==0)
   {
      tm << "If noheader option is given number of samples must be given.\n";
      Ace::assert<InvalidArg>(false,f,__LINE__);
   }
   Ace::assert<NotNewFile>(File::head()==fnullptr,f,__LINE__);
   std::ifstream file(ops.com_front());
   Ace::assert<CannotOpen>(file.is_open(),f,__LINE__);
   try
   {
      tm << "Importing header information...\n";
      read_headers(file,sampleSize,tr,hasHeaders);
      file.clear();
      file.seekg(0,std::ios_base::beg);
      read_gene_expressions(file,tm,nan);
      _isNew = false;
   }
   catch (...)
   {
      File::clear();
      null_data();
      _isNew = true;
      throw;
   }
}



void EMatrix::dump(Ace::GetOpts &ops, Ace::Terminal &tm)
{
   tm << "dump command not implemented for EMatrix type.\n";
}



void EMatrix::query(Ace::GetOpts &ops, Ace::Terminal &tm)
{
   enum {Basic=0,Lookup};
   if (empty())
   {
      tm << "No data loaded.\n";
   }
   else
   {
      switch (ops.com_get({"lookup","dump"}))
      {
      case Basic:
         tm << gene_size() << " gene(s) and " << sample_size() << " sample(s).\n";
         tm << "Sample Transformation: ";
         switch (transform())
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
         break;
      case Lookup:
         ops.com_pop();
         lookup(ops,tm);
         break;
      }
   }
}



bool EMatrix::empty()
{
   return _isNew;
}



void EMatrix::initialize(std::vector<std::string>&& geneNames,
                         std::vector<std::string>&& sampleNames, Transform transform)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<InvalidSize>(geneNames.size()>0&&sampleNames.size()>0,f,__LINE__);
   if (File::head()==fnullptr)
   {
      addr(fnullptr);
      allocate();
      File::head(addr());
   }
   else
   {
      Ace::assert<AlreadySet>(_isNew,f,__LINE__);
   }
   try
   {
      data()._geneSize = geneNames.size();
      data()._sampleSize = sampleNames.size();
      Ace::FString fstr(File::mem());
      auto i = geneNames.begin();
      data()._genePtr = fstr.write(*i);
      ++i;
      while (i!=geneNames.end())
      {
         fstr.reset();
         fstr.write(*i);
         ++i;
      }
      fstr.reset();
      i = sampleNames.begin();
      data()._samplePtr = fstr.write(*i);
      ++i;
      while (i!=sampleNames.end())
      {
         fstr.reset();
         fstr.write(*i);
         ++i;
      }
      data()._expData = Gene::initialize(File::rmem(),geneNames.size(),sampleNames.size());
      data()._transform = transform;
      write();
      _geneNames = std::move(geneNames);
      _sampleNames = std::move(sampleNames);
      _isNew = false;
   }
   catch (...)
   {
      File::head(fnullptr);
      null_data();
      _isNew = true;
   }
}



int EMatrix::gene_size() const
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   return data()._geneSize;
}



int EMatrix::sample_size() const
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   return data()._sampleSize;
}



const std::string& EMatrix::gene_name(int i) const
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   Ace::assert<OutOfRange>(i>=0&&i<data()._geneSize,f,__LINE__);
   return _geneNames[i];
}



const std::string& EMatrix::sample_name(int i) const
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   Ace::assert<OutOfRange>(i>=0&&i<data()._sampleSize,f,__LINE__);
   return _sampleNames[i];
}



EMatrix::Transform EMatrix::transform() const
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   return (Transform)data()._transform;
}



EMatrix::Gene EMatrix::begin()
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   return Gene(this,0);
}



EMatrix::Gene EMatrix::end()
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   return Gene(this,data()._geneSize);
}



EMatrix::Gene EMatrix::find(int i)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<OutOfRange>( i >= 0 && i < data()._geneSize ,f,__LINE__);
   return Gene(this,i);
}



EMatrix::Gene EMatrix::find(const std::string& name)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   int i {0};
   for (;i<_geneNames.size();++i)
   {
      if (_geneNames[i]==name)
      {
         break;
      }
   }
   return Gene(this,i);
}



void EMatrix::read_headers(std::ifstream& file, int sampleSize, Transform transform,
                           bool hasHeaders)
{
   static const char* f = __PRETTY_FUNCTION__;
   std::vector<std::string> sampleNames;
   std::vector<std::string> geneNames;
   skip_blanks(file);
   if (hasHeaders)
   {
      std::string buf;
      std::getline(file,buf);
      std::istringstream ibuf(buf);
      std::string tmp;
      while (ibuf >> tmp)
      {
         sampleNames.push_back(tmp);
      }
   }
   else
   {
      sampleNames.resize(sampleSize);
   }
   while (file)
   {
      std::string buf;
      std::getline(file,buf);
      if (!is_blank_line(buf))
      {
         std::string tmp;
         std::istringstream ibuf(buf);
         ibuf >> tmp;
         geneNames.push_back(tmp);
      }
   }
   initialize(std::move(geneNames),std::move(sampleNames),transform);
}



void EMatrix::read_gene_expressions(std::ifstream& file, Ace::Terminal& tm,
                                    const std::string& nanStr)
{
   static const char* f = __PRETTY_FUNCTION__;
   int totalGSize {data()._geneSize};
   int totalDone {0};
   tm << "Importing gene samples[0%]..." << Ace::Terminal::flush;
   Mirror m(*this);
   skip_blanks(file);
   while (file)
   {
      std::string buf;
      std::getline(file,buf);
      if (!is_blank_line(buf))
      {
         break;
      }
   }
   auto g = m.begin();
   while (file)
   {
      std::string buf;
      std::getline(file,buf);
      if (!is_blank_line(buf))
      {
         Ace::assert<InvalidFile>(g!=m.end(),f,__LINE__);
         std::istringstream ibuf(buf);
         std::string tmp;
         Ace::assert<InvalidFile>(ibuf >> tmp,f,__LINE__);
         for (auto i = g.begin();i!=g.end();++i)
         {
            Ace::assert<InvalidFile>(ibuf >> tmp,f,__LINE__);
            if (tmp==nanStr)
            {
               *i = NAN;
            }
            else
            {
               try
               {
                  switch (transform())
                  {
                  case none:
                     *i = std::stof(tmp);
                     break;
                  case log:
                     *i = logf(std::stof(tmp));
                     break;
                  case log2:
                     *i = log2f(std::stof(tmp));
                     break;
                  case log10:
                     *i = log10f(std::stof(tmp));
                     break;
                  }
               }
               catch (std::exception)
               {
                  Ace::assert<InvalidFile>(false,f,__LINE__);
               }
            }
         }
         ++g;
         ++totalDone;
         int percentDone {totalDone*100/totalGSize};
         tm << "\rImporting gene samples[" << percentDone << "%]..." << Ace::Terminal::flush;
      }
   }
   tm << "\n";
   Ace::assert<InvalidFile>(g==m.end(),f,__LINE__);
   m.write();
}



void EMatrix::lookup(Ace::GetOpts &ops, Ace::Terminal &tm)
{
   Gene l {end()};
   try
   {
      int i = std::stoi(ops.com_front());
      if (i<0||i>=gene_size())
      {
         tm << "Gene index is out of range.\n";
         return;
      }
      l = find(i);
   }
   catch (std::exception)
   {
      l = find(ops.com_front());
   }
   if (l!=end())
   {
      l.read();
      tm << l.name() << ": ";
      for (auto i = l.begin();i!=l.end();++i)
      {
         tm << *i << " ";
      }
      tm << "\n";
   }
   else
   {
      tm << "No such gene found.\n";
   }
}



void EMatrix::skip_blanks(std::ifstream& file)
{
   static const char* f = __PRETTY_FUNCTION__;
   while (std::isspace(file.peek()))
   {
      Ace::assert<InvalidFile>(file.good(),f,__LINE__);
      file.get();
   }
}



bool EMatrix::is_blank_line(const std::string& line)
{
   bool ret {true};
   for (auto i = line.begin();i!=line.end();++i)
   {
      if (*i=='#')
      {
         break;
      }
      else if (!std::isspace(*i))
      {
         ret = false;
         break;
      }
   }
   return ret;
}



void EMatrix::null_data()
{
   data()._sampleSize = 0;
   data()._geneSize = 0;
   data()._transform = none;
   data()._samplePtr = fnullptr;
   data()._genePtr = fnullptr;
   data()._expData = fnullptr;
}



void EMatrix::flip_endian()
{
   flip(0,4);
   flip(4,4);
   flip(9,8);
   flip(17,8);
   flip(25,8);
}



const std::string& EMatrix::Gene::name() const
{
   return _p->gene_name(_i);
}



void EMatrix::Gene::read()
{
   Node::read(_i);
}



void EMatrix::Gene::write()
{
   Node::write(_i);
}



EMatrix::Gene::Iterator EMatrix::Gene::begin()
{
   return Iterator(this,0);
}



EMatrix::Gene::Iterator EMatrix::Gene::end()
{
   return Iterator(this,_p->data()._sampleSize);
}



float& EMatrix::Gene::at(int i)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<OutOfRange>(i>=0&&i<_p->data()._sampleSize,f,__LINE__);
   return (*this)[i];
}



float& EMatrix::Gene::operator[](int i)
{
   return pget<float>()[i];
}



void EMatrix::Gene::operator++()
{
   if (_i<(_p->data()._geneSize))
   {
      ++_i;
   }
}



bool EMatrix::Gene::operator!=(const Gene& cmp)
{
   return _p!=cmp._p||_i!=cmp._i;
}



bool EMatrix::Gene::operator==(const Gene& cmp)
{
   return _p==cmp._p&&_i==cmp._i;
}



EMatrix::Gene::Gene(EMatrix* p, int i):
   Node(sizeof(float)*(p->data()._sampleSize),p->File::mem(),p->data()._expData),
   _p(p),
   _i(i)
{
   init_data<float>(p->data()._sampleSize);
}



void EMatrix::Gene::flip_endian()
{
   for (int i = 0;i<(_p->data()._sampleSize);++i)
   {
      flip(i*4,4);
   }
}



int64_t EMatrix::Gene::initialize(Ace::NVMemory& mem, int geneSize, int sampleSize)
{
   size_t nSize {sizeof(float)*geneSize*sampleSize};
   if ( nSize > mem.available() )
   {
      mem.reserve(nSize-mem.available());
   }
   return mem.allocate(nSize);
}



void EMatrix::Gene::Iterator::operator++()
{
   if (_i<(_p->_p->data()._sampleSize))
   {
      ++_i;
   }
}



float& EMatrix::Gene::Iterator::operator*()
{
   return (_p->pget<float>())[_i];;
}



bool EMatrix::Gene::Iterator::operator!=(const Iterator& cmp)
{
   return _p!=cmp._p||_i!=cmp._i;
}



bool EMatrix::Gene::Iterator::operator==(const Iterator& cmp)
{
   return _p==cmp._p&&_i==cmp._i;
}



EMatrix::Gene::Iterator::Iterator(Gene* p, int i):
   _p(p),
   _i(i)
{}



EMatrix::Mirror::Mirror(EMatrix& p):
   Node(sizeof(float)*p.data()._geneSize*p.data()._sampleSize,p.File::mem(),p.data()._expData),
   _p(&p)
{
   static const char* f = __PRETTY_FUNCTION__;
   AccelCompEng::assert<NullMatrix>(!(p._isNew),f,__LINE__);
   init_data<float>(p.data()._geneSize*p.data()._sampleSize);
}



void EMatrix::Mirror::read()
{
   Node::read();
}



void EMatrix::Mirror::write()
{
   Node::write();
}



float& EMatrix::Mirror::value(int gene, int i)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<OutOfRange>(gene>=0&&gene<(_p->data()._geneSize)&&i>=0&&
                           i<(_p->data()._sampleSize),f,__LINE__);
   return pget<float>()[(gene*(_p->data()._sampleSize))+i];
}



EMatrix::Mirror::Gene EMatrix::Mirror::begin()
{
   return Gene(this,0);
}



EMatrix::Mirror::Gene EMatrix::Mirror::end()
{
   return Gene(this,_p->data()._geneSize);
}



EMatrix::Mirror::Gene EMatrix::Mirror::find(int i)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<OutOfRange>(i>=0&&i<(_p->data()._geneSize),f,__LINE__);
   return Gene(this,i);
}



void EMatrix::Mirror::flip_endian()
{
   for (int i = 0;i<((_p->data()._sampleSize)*(_p->data()._geneSize));++i)
   {
      flip(i*4,4);
   }
}



EMatrix::Mirror::Gene::Iterator EMatrix::Mirror::Gene::begin()
{
   return Iterator(this,0);
}



EMatrix::Mirror::Gene::Iterator EMatrix::Mirror::Gene::end()
{
   return Iterator(this,_p->_p->data()._sampleSize);
}



float& EMatrix::Mirror::Gene::at(int i)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<OutOfRange>(i>=0&&i<(_p->_p->data()._geneSize),f,__LINE__);
   return (*this)[i];
}



float& EMatrix::Mirror::Gene::operator[](int i)
{
   return _p->value(_i,i);
}



void EMatrix::Mirror::Gene::operator++()
{
   if (_i<(_p->_p->data()._geneSize))
   {
      ++_i;
   }
}



bool EMatrix::Mirror::Gene::operator!=(const Gene& cmp)
{
   return _p!=cmp._p||_i!=cmp._i;
}



bool EMatrix::Mirror::Gene::operator==(const Gene& cmp)
{
   return _p==cmp._p&&_i==cmp._i;
}



EMatrix::Mirror::Gene::Gene(Mirror* p, int i):
   _p(p),
   _i(i)
{}



float& EMatrix::Mirror::Gene::Iterator::operator*()
{
   return _p->_p->value(_p->_i,_i);
}



void EMatrix::Mirror::Gene::Iterator::operator++()
{
   if (_i<(_p->_p->_p->data()._sampleSize))
   {
      ++_i;
   }
}



bool EMatrix::Mirror::Gene::Iterator::operator!=(const Iterator& cmp)
{
   return _p!=cmp._p||_i!=cmp._i;
}



bool EMatrix::Mirror::Gene::Iterator::operator==(const Iterator& cmp)
{
   return _p==cmp._p&&_i==cmp._i;
}



EMatrix::Mirror::Gene::Iterator::Iterator(Gene* p, int i):
   _p(p),
   _i(i)
{}
