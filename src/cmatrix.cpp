#include "cmatrix.h"



CMatrix::CMatrix():
   Node(sizeof(Header))
{
   init_data<Header>();
}



CMatrix::~CMatrix()
{
   if (_iGPair)
   {
      delete _iGPair;
   }
}



void CMatrix::init()
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
         fstr.load(data()._correlationPtr);
         _correlationNames.push_back(fstr.str());
         for (int i=1;i<data()._sampleSize;++i)
         {
            fstr.bump();
            _correlationNames.push_back(fstr.str());
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



void CMatrix::load(Ace::GetOpts&, Ace::Terminal& tm)
{
   tm << "Not yet implemented.\n";
}



void CMatrix::dump(Ace::GetOpts&, Ace::Terminal& tm)
{
   tm << "Not yet implemented.\n";
}



void CMatrix::query(Ace::GetOpts&, Ace::Terminal& tm)
{
   tm << "Not yet implemented.\n";
}



bool CMatrix::empty()
{
   return _isNew;
}



long long CMatrix::diagonal(int x, int y)
{
   return (x*(x-1)/2)+y;
}



long long CMatrix::diag_size(int x)
{
   return x*(x+1)/2;
}


void CMatrix::initialize(std::vector<std::string>&& geneNames,
                         std::vector<std::string>&& sampleNames,
                         std::vector<std::string>&& correlationNames, int maxModes)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<InvalidSize>(geneNames.size()>0&&sampleNames.size()>0&&
                            correlationNames.size()>0,f,__LINE__);
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
      data()._correlationSize = correlationNames.size();
      data()._maxModes = maxModes;
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
      fstr.reset();
      i = correlationNames.begin();
      data()._correlationPtr = fstr.write(*i);
      ++i;
      while (i!=correlationNames.end())
      {
         fstr.reset();
         fstr.write(*i);
         ++i;
      }
      data()._expData = GPair::initialize(File::rmem(),geneNames.size(),sampleNames.size(),
                                          correlationNames.size(),maxModes);
      write();
      _geneNames = std::move(geneNames);
      _sampleNames = std::move(sampleNames);
      _correlationNames = std::move(correlationNames);
      _isNew = false;
   }
   catch (...)
   {
      File::head(fnullptr);
      null_data();
      _isNew = true;
   }
}



int CMatrix::gene_size() const
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   return data()._geneSize;
}



int CMatrix::sample_size() const
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   return data()._sampleSize;
}



int CMatrix::max_modes() const
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   return data()._maxModes;
}



int CMatrix::correlation_size() const
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   return data()._correlationSize;
}



const std::string& CMatrix::gene_name(int i) const
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   return _geneNames[i];
}



const std::string& CMatrix::sample_name(int i) const
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   return _sampleNames[i];
}



const std::string& CMatrix::correlation_name(int i) const
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   return _correlationNames[i];
}


CMatrix::GPair CMatrix::begin()
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   return GPair(this,1,0);
}



CMatrix::GPair CMatrix::end()
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   return GPair(this,data()._geneSize,0);
}



CMatrix::GPair& CMatrix::at(int x, int y)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NullMatrix>(!_isNew,f,__LINE__);
   Ace::assert<OutOfRange>(x>=0&&y>=0&&y<data()._geneSize&&x<data()._geneSize&&x!=y,f,__LINE__);
   return ref(x,y);
}



CMatrix::GPair& CMatrix::ref(int x, int y)
{
   if ( x < y )
   {
      std::swap(x,y);
   }
   if (!_iGPair)
   {
      _iGPair = new GPair(this,x,y);
   }
   else
   {
      _iGPair->set(x,y);
   }
   return *_iGPair;
}



void CMatrix::null_data()
{
   data()._geneSize = 0;
   data()._sampleSize = 0;
   data()._correlationSize = 0;
   data()._maxModes = 0;
   data()._genePtr = fnullptr;
   data()._samplePtr = fnullptr;
   data()._correlationPtr = fnullptr;
   data()._expData = fnullptr;
}



void CMatrix::flip_endian()
{
   flip(0,4);
   flip(4,4);
   flip(8,4);
   flip(13,4);
   flip(21,4);
   flip(29,4);
   flip(37,4);
}



void CMatrix::GPair::read()
{
   Node::read(diagonal(_x,_y));
}



void CMatrix::GPair::write()
{
   Node::write(diagonal(_x,_y));
}



int CMatrix::GPair::x() const
{
   return _x;
}



int CMatrix::GPair::y() const
{
   return _y;
}



int CMatrix::GPair::size() const
{
   return *(reinterpret_cast<const int8_t*>(pget<char>()));
}



void CMatrix::GPair::size(int size)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<GreaterThanMax>(size>=0&&size<=(_p->data()._maxModes),f,__LINE__);
   *(reinterpret_cast<int8_t*>(pget<char>())) = size;
}



CMatrix::GPair::Modes& CMatrix::GPair::modes()
{
   return _modes;
}



CMatrix::GPair::Corrs& CMatrix::GPair::corrs()
{
   return _corrs;
}



void CMatrix::GPair::operator++()
{
   if (_x<(_p->data()._geneSize))
   {
      if ((_y+1)<_x)
      {
         _y++;
      }
      else
      {
         _y = 0;
         ++_x;
      }
   }
}



bool CMatrix::GPair::operator!=(const GPair& cmp)
{
   return _p!=cmp._p||_x!=cmp._x||_y!=cmp._y;
}



CMatrix::GPair::GPair(CMatrix* p, int x, int y):
   Node(calc_size(p->data()._maxModes,p->data()._sampleSize,p->data()._correlationSize),
        p->File::mem(),p->data()._expData),
   _p(p),
   _x(x),
   _y(y),
   _modes(this),
   _corrs(this)
{}



void CMatrix::GPair::set(int x, int y)
{
   bool change {_x!=x||_y!=y};
   _x = x;
   _y = y;
   if (change)
   {
      read();
   }
}



int8_t& CMatrix::GPair::mode_val(int mi, int i)
{
   size_t inc {sizeof(int8_t)*((mi*(_p->data()._sampleSize))+i+1)};
   return *reinterpret_cast<int8_t*>(&(pget<char>()[inc]));
}



float& CMatrix::GPair::correlation_val(int mi, int ci)
{
   size_t inc {sizeof(int8_t)*(((_p->data()._maxModes)*(_p->data()._sampleSize))+1)+
            sizeof(float)*((mi*(_p->data()._maxModes))+ci)};
   return *reinterpret_cast<float*>(&(pget<char>()[inc]));
}



int64_t CMatrix::GPair::initialize(Ace::NVMemory& mem, int geneSize, int maxModes, int sampleSize,
                                   int correlationSize)
{
   return mem.allocate(calc_size(maxModes,sampleSize,correlationSize)*diag_size(geneSize));
}



size_t CMatrix::GPair::calc_size(int maxModes, int sampleSize, int correlationSize)
{
   return sizeof(int8_t)+
          (sizeof(int8_t)*maxModes*sampleSize)+
          (sizeof(float)*maxModes*correlationSize);
}



void CMatrix::GPair::flip_endian()
{
   size_t inc {sizeof(int8_t)*(((_p->data()._maxModes)*(_p->data()._sampleSize))+1)};
   for (int i=0;i<((_p->data()._maxModes)*(_p->data()._correlationSize));++i)
   {
      flip(inc+(sizeof(float)*i),4);
   }
}



CMatrix::GPair::Modes::Modes(const Modes& copy):
   _p(copy._p)
{
   if (copy._mode.get())
   {
      _mode.reset(new Mode(this,copy._mode->_i));
   }
}



CMatrix::GPair::Modes& CMatrix::GPair::Modes::operator=(const Modes& copy)
{
   _p = copy._p;
   if (copy._mode.get())
   {
      _mode.reset(new Mode(this,copy._mode->_i));
   }
}



CMatrix::GPair::Modes::Mode CMatrix::GPair::Modes::begin()
{
   return Mode(this,0);
}



CMatrix::GPair::Modes::Mode CMatrix::GPair::Modes::end()
{
   return Mode(this,_p->size());
}



CMatrix::GPair::Modes::Mode& CMatrix::GPair::Modes::at(int i)
{
   static const char* f = __PRETTY_FUNCTION__;
   AccelCompEng::assert<OutOfRange>(i>=0&&i<(_p->size()),f,__LINE__);
   return (*this)[i];
}



CMatrix::GPair::Modes::Mode& CMatrix::GPair::Modes::operator[](int i)
{
   if (!_mode.get())
   {
      _mode.reset(new Mode(this,i));
   }
   else
   {
      _mode->set(i);
   }
   return *_mode;
}



CMatrix::GPair::Modes::Modes(GPair* p):
   _p(p)
{}



CMatrix::GPair::Modes::Mode::Iterator CMatrix::GPair::Modes::Mode::begin()
{
   return Iterator(this,0);
}



CMatrix::GPair::Modes::Mode::Iterator CMatrix::GPair::Modes::Mode::end()
{
   return Iterator(this,_p->_p->_p->data()._sampleSize);
}



int8_t& CMatrix::GPair::Modes::Mode::at(int i)
{
   static const char* f = __PRETTY_FUNCTION__;
   AccelCompEng::assert<OutOfRange>(i>=0&&i<(_p->_p->_p->data()._sampleSize),f,__LINE__);
   return (*this)[i];
}



int8_t& CMatrix::GPair::Modes::Mode::operator[](int i)
{
   return _p->_p->mode_val(_i,i);
}



void CMatrix::GPair::Modes::Mode::operator++()
{
   if (_i<(_p->_p->size()))
   {
      ++_i;
   }
}



bool CMatrix::GPair::Modes::Mode::operator!=(const Mode& cmp)
{
   return _p!=cmp._p||_i!=cmp._i;
}



CMatrix::GPair::Modes::Mode::Mode(Modes* p, int i):
   _p(p),
   _i(i)
{}



void CMatrix::GPair::Modes::Mode::set(int i)
{
   _i = i;
}



int8_t& CMatrix::GPair::Modes::Mode::Iterator::operator*()
{
   return _p->_p->_p->mode_val(_p->_i,_i);
}



void CMatrix::GPair::Modes::Mode::Iterator::operator++()
{
   if (_i<(_p->_p->_p->_p->data()._sampleSize))
   {
      ++_i;
   }
}



bool CMatrix::GPair::Modes::Mode::Iterator::operator!=(const Iterator& cmp)
{
   return _p!=cmp._p||_i!=cmp._i;
}



CMatrix::GPair::Modes::Mode::Iterator::Iterator(Mode* p, int i):
   _p(p),
   _i(i)
{}



CMatrix::GPair::Corrs::Corrs(const Corrs& copy):
   _p(copy._p)
{
   if (copy._corr.get())
   {
      _corr.reset(new Corr(this,copy._corr->_i));
   }
}



CMatrix::GPair::Corrs& CMatrix::GPair::Corrs::operator=(const Corrs& copy)
{
   _p = copy._p;
   if (copy._corr.get())
   {
      _corr.reset(new Corr(this,copy._corr->_i));
   }
}



CMatrix::GPair::Corrs::Corr CMatrix::GPair::Corrs::begin()
{
   return Corr(this,0);
}



CMatrix::GPair::Corrs::Corr CMatrix::GPair::Corrs::end()
{
   return Corr(this,_p->size());
}



CMatrix::GPair::Corrs::Corr& CMatrix::GPair::Corrs::at(int i)
{
   static const char* f = __PRETTY_FUNCTION__;
   AccelCompEng::assert<OutOfRange>(i>=0&&i<(_p->size()),f,__LINE__);
   return (*this)[i];
}



CMatrix::GPair::Corrs::Corr& CMatrix::GPair::Corrs::operator[](int i)
{
   if (!_corr.get())
   {
      _corr.reset(new Corr(this,i));
   }
   else
   {
      _corr->set(i);
   }
   return *_corr;
}



CMatrix::GPair::Corrs::Corrs(GPair* p):
   _p(p)
{}



CMatrix::GPair::Corrs::Corr::Iterator CMatrix::GPair::Corrs::Corr::begin()
{
   return Iterator(this,0);
}



CMatrix::GPair::Corrs::Corr::Iterator CMatrix::GPair::Corrs::Corr::end()
{
   return Iterator(this,_p->_p->_p->data()._correlationSize);
}



float& CMatrix::GPair::Corrs::Corr::at(int i)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<OutOfRange>(i>=0&&i<(_p->_p->_p->data()._correlationSize),f,__LINE__);
   return (*this)[i];
}



float& CMatrix::GPair::Corrs::Corr::operator[](int i)
{
   return _p->_p->correlation_val(_i,i);
}



void CMatrix::GPair::Corrs::Corr::operator++()
{
   if (_i<(_p->_p->size()))
   {
      ++_i;
   }
}



bool CMatrix::GPair::Corrs::Corr::operator!=(const Corr& cmp)
{
   return _p!=cmp._p||_i!=cmp._i;
}



CMatrix::GPair::Corrs::Corr::Corr(Corrs* p, int i):
   _p(p),
   _i(i)
{}



void CMatrix::GPair::Corrs::Corr::set(int i)
{
   _i = i;
}



float& CMatrix::GPair::Corrs::Corr::Iterator::operator*()
{
   return _p->_p->_p->correlation_val(_p->_i,_i);
}



void CMatrix::GPair::Corrs::Corr::Iterator::operator++()
{
   if (_i<(_p->_p->_p->_p->data()._correlationSize))
   {
      ++_i;
   }
}



bool CMatrix::GPair::Corrs::Corr::Iterator::operator!=(const Iterator& cmp)
{
   return _p!=cmp._p||_i!=cmp._i;
}



CMatrix::GPair::Corrs::Corr::Iterator::Iterator(Corr* p, int i):
   _p(p),
   _i(i)
{}
