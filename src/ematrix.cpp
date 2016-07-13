#include "ematrix.h"
#include <sstream>
#include <cmath>



EMatrix::EMatrix(const string& type, const string& file):
   DataPlugin(type,file)
{
   try
   {
      // Does the KINC file have a header already set. This case would be for
      // a new file.
      if (File::head()==fNullPtr)
      {
         // Allocate space in the file (file memory) for the header.
         File::mem().allot(_hdr);
         File::head(_hdr.addr());
         // Set defaults for the ematrix header.
         _hdr.sSize() = 0;
         _hdr.gSize() = 0;
         _hdr.tr() = Transform::none;
         _hdr.wr() = 0;
         _hdr.sPtr() = fNullPtr;
         _hdr.gPtr() = fNullPtr;
         _hdr.eData() = fNullPtr;
         // Write the header initially to file.
         File::mem().sync(_hdr,FSync::write);
      }
      // The KINC file already has a header if the file already exists and is
      // being opened.
      else
      {
         _hdr.addr(File::head());
         File::mem().sync(_hdr,FSync::read);
         NmHead hd {_hdr.gPtr()};
         for (int i=0;i<_hdr.gSize();++i)
         {
            File::mem().sync(hd,FSync::read,i);
            FString name(&File::mem(),hd.nPtr());
            _gNames.push_back(*name);
         }
         if (_hdr.sPtr()!=fNullPtr)
         {
            hd.addr(_hdr.sPtr());
            for (int i=0;i<_hdr.sSize();++i)
            {
               File::mem().sync(hd,FSync::read,i);
               FString name(&File::mem(),hd.nPtr());
               _sNames.push_back(*name);
            }
         }
         else
         {
            _sNames.resize(_hdr.sSize());
         }
      }
   }
   catch (...)
   {
      File::clear();
      throw;
   }
}



EMatrix::~EMatrix()
{
   if (_iGene)
   {
      delete _iGene;
   }
}



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
   AccelCompEng::assert<NotNewFile>(File::head()==fNullPtr,__LINE__);
   std::ifstream f(ops.com_front());
   AccelCompEng::assert<CannotOpen>(f.is_open(),__LINE__);
   read_sizes(f,sampleSize);
   transform(tr);
   try
   {
      f.clear();
      f.seekg(0,std::ios_base::beg);
      while (f.peek()==' '||f.peek()=='\t'||std::iscntrl(f.peek()))
      {
         AccelCompEng::assert<InvalidFile>(f,__LINE__);
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
   }
}



void EMatrix::dump(GetOpts &ops, Terminal &tm)
{
   tm << "dump command not implemented for EMatrix type.\n";
}



void EMatrix::query(GetOpts &ops, Terminal &tm)
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
   }
}



bool EMatrix::empty()
{
   return _hdr.gSize()==0;
}



void EMatrix::initialize(int gSize, int sSize, bool sHeaders)
{
   if (File::head()==fNullPtr)
   {
      File::mem().allot(_hdr);
      File::head(_hdr.addr());
   }
   else
   {
      bool cond {_hdr.gPtr()==fNullPtr};
      AccelCompEng::assert<AlreadyInitialized>(cond,__LINE__);
   }
   NmHead hd;
   File::mem().allot(hd,gSize);
   _hdr.gPtr() = hd.addr();
   if (sHeaders)
   {
      File::mem().allot(hd,sSize);
      _hdr.sPtr() = hd.addr();
   }
   _hdr.gSize() = gSize;
   _hdr.sSize() = sSize;
   _hdr.eData() = Gene::initialize(File::mem(),gSize,sSize);
   _gNames.resize(gSize);
   _sNames.resize(sSize);
   File::mem().sync(_hdr,FSync::write);
   File::head(_hdr.addr());
}



bool EMatrix::hasSampleHead() const
{
   return _hdr.sPtr()!=fNullPtr;
}



int EMatrix::gSize() const
{
   return _hdr.gSize();
}



int EMatrix::sSize() const
{
   return _hdr.sSize();
}



const EMatrix::string& EMatrix::gName(int i) const
{
   bool cond {_hdr.gPtr()!=fNullPtr};
   AccelCompEng::assert<NoData>(cond,__LINE__);
   cond = {i>=0||i<_hdr.gSize()};
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
   return _gNames[i];
}



void EMatrix::gName(int i, const string& name)
{
   bool cond {_hdr.gPtr()!=fNullPtr};
   AccelCompEng::assert<NotInitialized>(cond,__LINE__);
   cond = {i>=0||i<_hdr.gSize()};
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
   _gNames[i] = name;
}



const EMatrix::string& EMatrix::sName(int i) const
{
   bool cond {_hdr.gPtr()!=fNullPtr};
   AccelCompEng::assert<NoData>(cond,__LINE__);
   cond = {i>=0||i<_hdr.gSize()};
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
   return _sNames[i];
}



void EMatrix::sName(int i, const string& name)
{
   bool cond {_hdr.sPtr()!=fNullPtr};
   AccelCompEng::assert<NotInitialized>(cond,__LINE__);
   cond = {i>=0||i<_hdr.sSize()};
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
   _sNames[i] = name;
}



EMatrix::Transform EMatrix::transform() const
{
   return (Transform)_hdr.tr();
}



void EMatrix::transform(Transform tr)
{
   _hdr.tr() = tr;
}



void EMatrix::write()
{
   bool cond {_hdr.gPtr()!=fNullPtr};
   AccelCompEng::assert<NotInitialized>(cond,__LINE__);
   NmHead hd {_hdr.gPtr()};
   cond = _hdr.wr()==0;
   AccelCompEng::assert<AlreadySet>(cond,__LINE__);
   _hdr.wr() = 1;
   for (int i=0;i<_gNames.size();++i)
   {
      FString tmp {&File::mem()};
      tmp = _gNames[i];
      hd.nPtr() = tmp.addr();
      File::mem().sync(hd,FSync::write,i);
   }
   if (_hdr.sPtr()!=fNullPtr)
   {
      hd.addr(_hdr.sPtr());
      for (int i=0;i<_sNames.size();++i)
      {
         FString tmp {&File::mem()};
         tmp = _sNames[i];
         hd.nPtr() = tmp.addr();
         File::mem().sync(hd,FSync::write,i);
      }
   }
   File::mem().sync(_hdr,FSync::write);
}



EMatrix::Gene EMatrix::begin()
{
   return Gene(this,0);
}



EMatrix::Gene EMatrix::end()
{
   return Gene(this,_hdr.gSize());
}



EMatrix::Gene EMatrix::find(int i)
{
   bool cond {_hdr.gPtr()!=fNullPtr};
   AccelCompEng::assert<NotInitialized>(cond,__LINE__);
   cond = i>=0&&i<_hdr.gSize();
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
   return Gene(this,i);
}



EMatrix::Gene EMatrix::find(const string& nm)
{
   int i {0};
   for (;i<_hdr.gSize();++i)
   {
      if (gName(i)==nm)
      {
         break;
      }
   }
   return Gene(this,i);
}



EMatrix::Gene& EMatrix::at(int i)
{
   bool cond {_hdr.gPtr()!=fNullPtr};
   AccelCompEng::assert<NotInitialized>(cond,__LINE__);
   cond = i>=0&&i<_hdr.gSize();
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
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



const EMatrix::Gene::string& EMatrix::Gene::name() const
{
   return _p->gName(_i);
}



void EMatrix::Gene::read()
{
   _p->File::mem().sync(_head,FSync::read,_i);
}



void EMatrix::Gene::write()
{
   _p->File::mem().sync(_head,FSync::write,_i);
}



EMatrix::Gene::Iterator EMatrix::Gene::begin()
{
   return Iterator(this,0);
}



EMatrix::Gene::Iterator EMatrix::Gene::end()
{
   return Iterator(this,_p->_hdr.sSize());
}



float& EMatrix::Gene::at(int i)
{
   bool cond {_p->_hdr.sPtr()!=fNullPtr};
   AccelCompEng::assert<NotInitialized>(cond,__LINE__);
   cond = i>=0&&i<_p->_hdr.sSize();
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
   return (*this)[i];
}



float& EMatrix::Gene::operator[](int i)
{
   return _head.val(i);
}



void EMatrix::Gene::operator++()
{
   if (_i<_p->_hdr.gSize())
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
   _p(p),
   _head(p->_hdr.sSize(),p->_hdr.eData()),
   _i(i)
{}



void EMatrix::Gene::set(int i)
{
   _i = i;
}



EMatrix::FPtr EMatrix::Gene::initialize(FileMem& mem, int gSize, int sSize)
{
   Expr newData(sSize);
   mem.allot(newData,gSize);
   return newData.addr();
}



void EMatrix::Gene::Iterator::operator++()
{
   if (_i<_p->_p->_hdr.sSize())
   {
      ++_i;
   }
}



float& EMatrix::Gene::Iterator::operator*()
{
   return _p->_head.val(_i);
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
   _p(&p),
   _data(p._hdr.gSize(),p._hdr.sSize(),p._hdr.eData())
{
   bool cond {p._hdr.eData()!=fNullPtr};
   AccelCompEng::assert<NoData>(cond,__LINE__);
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
   _p->File::mem().sync(_data,FSync::read);
}



void EMatrix::Mirror::write()
{
   _p->File::mem().sync(_data,FSync::write);
}



EMatrix::Mirror::Gene EMatrix::Mirror::begin()
{
   return Gene(this,0);
}



EMatrix::Mirror::Gene EMatrix::Mirror::end()
{
   return Gene(this,_p->_hdr.gSize());
}



EMatrix::Mirror::Gene& EMatrix::Mirror::at(int i)
{
   bool cond {i>=0&&i<_p->_hdr.gSize()};
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
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
   return Iterator(this,_p->_p->_hdr.sSize());
}



float& EMatrix::Mirror::Gene::at(int i)
{
   bool cond {i>=0&&i<(_p->_p->_hdr.gSize())};
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
   return (*this)[i];
}



float& EMatrix::Mirror::Gene::operator[](int i)
{
   return _p->_data.val(_i,i);
}



void EMatrix::Mirror::Gene::operator++()
{
   if (_i<(_p->_p->_hdr.gSize()))
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
   return _p->_p->_data.val(_p->_i,_i);
}



void EMatrix::Mirror::Gene::Iterator::operator++()
{
   if (_i<(_p->_p->_p->_hdr.sSize()))
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



void EMatrix::read_sizes(std::ifstream& f, int sSize)
{
   bool hasHeaders {sSize==0};
   while (f.peek()==' '||f.peek()=='\t'||std::iscntrl(f.peek()))
   {
      AccelCompEng::assert<InvalidFile>(f,__LINE__);
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
      AccelCompEng::assert<InvalidFile>(f,__LINE__);
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
      sName(i++,tmp);
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
         gName(gi++,tmp);
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
