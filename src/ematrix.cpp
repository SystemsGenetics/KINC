#include "ematrix.h"
#include <sstream>
#include <cmath>



EMatrix::~EMatrix()
{
   if (_iGene)
   {
      delete _iGene;
   }
}



void EMatrix::init()
{
   try
   {
      _header.mem(Ace::File::mem());
      if (Ace::File::head()==fnullptr)
      {
         _header.allocate();
         Ace::File::head(_header.addr());
         _header.null_data();
         _header.write();
      }
      else
      {
         _header.addr(Ace::File::head());
         _header.read();
         Ace::FString fstr(Ace::File::mem(),_header.data()._genePtr);
         fstr.static_buffer(_strSize);
         _geneNames.push_back(fstr.str());
         for (int i=1;i<_header.data()._geneSize;++i)
         {
            fstr.bump();
            _geneNames.push_back(fstr.str());
         }
         if (_header.data()._samplePtr!=fnullptr)
         {
            fstr.load(_header.data()._samplePtr);
            _sampleNames.push_back(fstr.str());
            for (int i=1;i<_header.data()._sampleSize;++i)
            {
               fstr.bump();
               _sampleNames.push_back(fstr.str());
            }
         }
         else
         {
            _sampleNames.resize(_header.data()._sampleSize);
         }
         _isNew = false;
      }
   }
   catch (...)
   {
      File::clear();
      throw;
   }
}



void EMatrix::load(Ace::GetOpts &ops, Ace::Terminal &tm)
{
/*   static const char* f = __PRETTY_FUNCTION__;
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
   Ace::assert<NotNewFile>(File::head()==fNullPtr,f,__LINE__);
   std::ifstream f(ops.com_front());
   Ace::assert<CannotOpen>(f.is_open(),f,__LINE__);
   read_sizes(f,sampleSize);
   transform(tr);
   try
   {
      f.clear();
      f.seekg(0,std::ios_base::beg);
      while (f.peek()==' '||f.peek()=='\t'||std::iscntrl(f.peek()))
      {
         Ace::assert<InvalidFile>(f.good(),f,__LINE__);
         f.get();
      }
      if (hasHeaders)
      {
         read_header(f);
      }
      read_gene_expressions(f,"NaN");
      write();
      tm << gSize() << "\n";
      tm << sSize() << "\n";
   }
   catch (...)
   {
      File::clear();
      throw;
   }*/
}



void EMatrix::dump(Ace::GetOpts &ops, Ace::Terminal &tm)
{
   tm << "dump command not implemented for EMatrix type.\n";
}



void EMatrix::query(Ace::GetOpts &ops, Ace::Terminal &tm)
{
/*   enum {Basic=0,Lookup};
   if (empty())
   {
      tm << "No data loaded.\n";
   }
   else
   {
      switch (ops.com_get({"lookup","dump"}))
      {
      case Basic:
         tm << _hdr.gSize() << " gene(s) and " << _hdr.sSize()
            << " sample(s).\n";
         tm << "Sample Transformation: ";
         switch (_hdr.tr())
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
         if (hasSampleHead())
         {
            tm << "Sample header data exists.\n";
         }
         else
         {
            tm << "Sample header data does not exist.\n";
         }
         break;
      case Lookup:
         ops.com_pop();
         lookup(ops,tm);
         break;
      }
   }*/
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
      _header.addr(fnullptr);
      _header.allocate();
      Ace::File::head(_header.addr());
   }
   else
   {
      Ace::assert<AlreadySet>(_isNew,f,__LINE__);
   }
   try
   {
      _header.data()._geneSize = geneNames.size();
      _header.data()._sampleSize = sampleNames.size();
      Ace::FString fstr(Ace::File::mem());
      auto i = geneNames.begin();
      _header.data()._genePtr = fstr.write(*i);
      ++i;
      while (i!=geneNames.end())
      {
         fstr.reset();
         fstr.write(*i);
         ++i;
      }
      fstr.reset();
      i = sampleNames.begin();
      _header.data()._samplePtr = fstr.write(*i);
      ++i;
      while (i!=sampleNames.end())
      {
         fstr.reset();
         fstr.write(*i);
         ++i;
      }
      _header.data()._expData = Gene::initialize(Ace::File::rmem(),geneNames.size(),
                                                 sampleNames.size());
      _header.data()._transform = transform;
      _header.write();
      _geneNames = std::move(geneNames);
      _sampleNames = std::move(sampleNames);
      _isNew = false;
   }
   catch (...)
   {
      Ace::File::head(fnullptr);
      _isNew = true;
   }
}



int EMatrix::gene_size() const
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   return _header.data()._geneSize;
}



int EMatrix::sample_size() const
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   return _header.data()._sampleSize;
}



const std::string& EMatrix::gene_name(int i) const
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   Ace::assert<OutOfRange>(i>=0&&i<_header.data()._geneSize,f,__LINE__);
   return _geneNames[i];
}



const std::string& EMatrix::sample_name(int i) const
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   Ace::assert<OutOfRange>(i>=0&&i<_header.data()._sampleSize,f,__LINE__);
   return _sampleNames[i];
}



EMatrix::Transform EMatrix::transform() const
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   return (Transform)_header.data()._transform;
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
   return Gene(this,_header.data()._geneSize);
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



EMatrix::Gene& EMatrix::at(int i)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   Ace::assert<OutOfRange>(i>=0&&i<_header.data()._geneSize,f,__LINE__);
   return (*this)[i];
}



EMatrix::Gene& EMatrix::operator[](int i)
{
   if (!_iGene)
   {
      _iGene = new Gene(this,i);
   }
   else
   {
      _iGene->set(i);
   }
   return *_iGene;
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
   return Iterator(this,_p->_header.data()._sampleSize);
}



float& EMatrix::Gene::at(int i)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<OutOfRange>(i>=0&&i<_p->_header.data()._sampleSize,f,__LINE__);
   return (*this)[i];
}



float& EMatrix::Gene::operator[](int i)
{
   return pget<float>()[i];
}



void EMatrix::Gene::operator++()
{
   if (_i<(_p->_header.data()._geneSize))
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
   Node(sizeof(float)*(p->_header.data()._sampleSize),p->Ace::File::mem(),
        p->_header.data()._expData),
   _p(p),
   _i(i)
{
   init_data<float>(p->_header.data()._sampleSize);
   read();
}



void EMatrix::Gene::set(int i)
{
   bool change {_i!=i};
   _i = i;
   if (change)
   {
      read();
   }
}



void EMatrix::Gene::flip_endian()
{
   for (int i = 0;i<(_p->_header.data()._sampleSize);++i)
   {
      flip(i*4,4);
   }
}



int64_t EMatrix::Gene::initialize(Ace::NVMemory& mem, int geneSize, int sampleSize)
{
   return mem.allocate(sizeof(float)*geneSize*sampleSize);
}



void EMatrix::Gene::Iterator::operator++()
{
   if (_i<(_p->_p->_header.data()._sampleSize))
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


/////////////////
EMatrix::Mirror::Mirror(EMatrix& p):
   Node(sizeof(float)*p._header.data()._geneSize*p._header.data()._sampleSize,p.Ace::File::mem(),
        p._header.data()._expData),
   _p(&p)
{
   static const char* f = __PRETTY_FUNCTION__;
   AccelCompEng::assert<NullMatrix>(!(p._isNew),f,__LINE__);
   init_data<float>(p._header.data()._geneSize*p._header.data()._sampleSize);
}



EMatrix::Mirror::~Mirror()
{
   if (_iGene)
   {
      delete _iGene;
   }
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
   Ace::assert<OutOfRange>(gene>=0&&gene<(_p->_header.data()._geneSize)&&i>=0&&
                           i<(_p->_header.data()._sampleSize),f,__LINE__);
   return pget<float>()[(gene*(_p->_header.data()._sampleSize))+i];
}



EMatrix::Mirror::Gene EMatrix::Mirror::begin()
{
   return Gene(this,0);
}



EMatrix::Mirror::Gene EMatrix::Mirror::end()
{
   return Gene(this,_p->_header.data()._geneSize);
}



EMatrix::Mirror::Gene& EMatrix::Mirror::at(int i)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<OutOfRange>(i>=0&&i<(_p->_header.data()._geneSize),f,__LINE__);
   return (*this)[i];
}



EMatrix::Mirror::Gene& EMatrix::Mirror::operator[](int i)
{
   if (!_iGene)
   {
      _iGene = new Gene(this,i);
   }
   else
   {
      _iGene->set(i);
   }
   return *_iGene;
}




EMatrix::Mirror::Gene::Iterator EMatrix::Mirror::Gene::begin()
{
   return Iterator(this,0);
}



EMatrix::Mirror::Gene::Iterator EMatrix::Mirror::Gene::end()
{
   return Iterator(this,_p->_p->_header.data()._sampleSize);
}



float& EMatrix::Mirror::Gene::at(int i)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<OutOfRange>(i>=0&&i<(_p->_p->_header.data()._geneSize),f,__LINE__);
   return (*this)[i];
}



float& EMatrix::Mirror::Gene::operator[](int i)
{
   return _p->value(_i,i);
}



void EMatrix::Mirror::Gene::operator++()
{
   if (_i<(_p->_p->_header.data()._geneSize))
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



void EMatrix::Mirror::Gene::set(int i)
{
   _i = i;
}



float& EMatrix::Mirror::Gene::Iterator::operator*()
{
   return _p->_p->value(_p->_i,_i);
}



void EMatrix::Mirror::Gene::Iterator::operator++()
{
   if (_i<(_p->_p->_p->_header.data()._sampleSize))
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


/*
void EMatrix::read_sizes(std::ifstream& f, int sSize)
{
   bool hasHeaders {sSize==0};
   while (f.peek()==' '||f.peek()=='\t'||std::iscntrl(f.peek()))
   {
      AccelCompEng::assert<InvalidFile>(f.good(),__LINE__);
      f.get();
   }
   if (hasHeaders)
   {
      int count {1};
      enum {newline,chunk} state {chunk};
      while (std::isprint(f.peek()))
      {
         switch (state)
         {
         case newline:
            if (f.peek()!=' '&&f.peek()!='\t')
            {
               ++count;
               state = chunk;
            }
            break;
         case chunk:
            if (f.peek()==' '||f.peek()=='\t')
            {
               state = newline;
            }
            break;
         }
         f.get();
      }
      sSize = count;
      AccelCompEng::assert<InvalidFile>(f.good(),__LINE__);
   }
   int gSize {0};
   enum {newline,chunk} state {newline};
   while (f)
   {
      switch (state)
      {
      case newline:
         if (std::isprint(f.peek()))
         {
            state = chunk;
         }
         break;
      case chunk:
         if (std::iscntrl(f.peek()))
         {
            ++gSize;
            state = newline;
         }
         break;
      }
      f.get();
   }
   if (state==chunk)
   {
      ++gSize;
   }
   initialize(gSize,sSize,hasHeaders);
}



void EMatrix::read_header(std::ifstream& f)
{
   string buf;
   std::getline(f,buf);
   std::istringstream ibuf(buf);
   string tmp;
   int i {0};
   while (ibuf >> tmp)
   {
      sName(i++) = tmp;
   }
   bool cond {i==_hdr.sSize()};
   AccelCompEng::assert<InvalidFile>(cond,__LINE__);
}



void EMatrix::read_gene_expressions(std::ifstream& f, const string& nan)
{
   Mirror m(*this);
   auto g = m.begin();
   int gi {0};
   while (f)
   {
      string buf;
      std::getline(f,buf);
      if (!buf.empty())
      {
         if (g==m.end())
         {
            throw InvalidFile(__LINE__);
         }
         std::istringstream ibuf(buf);
         string tmp;
         if (!(ibuf >> tmp))
         {
            throw InvalidFile(__LINE__);
         }
         gName(gi++) = tmp;
         for (auto i = g.begin();i!=g.end();++i)
         {
            if (!(ibuf >> tmp))
            {
               throw InvalidFile(__LINE__);
            }
            if (tmp==nan)
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
                  throw InvalidFile(__LINE__);
               }
            }
         }
         ++g;
      }
   }
   if (g!=m.end())
   {
      throw InvalidFile(__LINE__);
   }
   m.write();
}



void EMatrix::lookup(GetOpts &ops, Terminal &tm)
{
   Gene l {end()};
   try
   {
      int i = std::stoi(ops.com_front());
      if (i<0||i>=gSize())
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
*/
