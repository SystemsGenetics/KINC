#include "cmatrix.h"


/*
CMatrix::CMatrix(const string& type, const string& file):
   DataPlugin(type,file)
{
   if (File::head()==FileMem::nullPtr)
   {
      File::mem().allot(_hdr);
      File::head(_hdr.addr());
      _hdr.gSize() = 0;
      _hdr.sSize() = 0;
      _hdr.cSize() = 0;
      _hdr.mSize() = 0;
      _hdr.wr() = 0;
      _hdr.gPtr() = fNullPtr;
      _hdr.sPtr() = fNullPtr;
      _hdr.cPtr() = fNullPtr;
      _hdr.eData() = fNullPtr;
      File::mem().sync(_hdr,FSync::write);
   }
   else
   {
      _hdr = File::head();
      File::mem().sync(_hdr,FSync::read);
      NmHead nm {_hdr.gPtr()};
      for (int i = 0;i<_hdr.gSize();++i)
      {
         File::mem().sync(nm,FSync::read,i);
         FString tmp(&File::mem(),nm.nPtr());
         _gNames.push_back(*tmp);
      }
      nm.addr(_hdr.cPtr());
      for (int i = 0;i<_hdr.cSize();++i)
      {
         File::mem().sync(nm,FSync::read,i);
         FString tmp(&File::mem(),nm.nPtr());
         _cNames.push_back(*tmp);
      }
      if (_hdr.sPtr()!=fNullPtr)
      {
         nm.addr(_hdr.sPtr());
         for (int i = 0;i<_hdr.sSize();++i)
         {
            File::mem().sync(nm,FSync::read,i);
            FString tmp(&File::mem(),nm.nPtr());
            _sNames.push_back(*tmp);
         }
      }
      else
      {
         _sNames.resize(_hdr.sSize());
      }
   }
}



void CMatrix::load(GetOpts&, Terminal& tm)
{
   tm << "Not yet implemented.\n";
}



void CMatrix::dump(GetOpts&,Terminal& tm)
{
   tm << "Not yet implemented.\n";
}



void CMatrix::query(GetOpts&,Terminal& tm)
{
   tm << "Not yet implemented.\n";
}



bool CMatrix::empty()
{
   return _hdr.gPtr()==fNullPtr;
}



void CMatrix::initialize(int gSize, int sSize, int mSize, int cSize, bool sHdrs)
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
   File::mem().allot(hd,cSize);
   _hdr.cPtr() = hd.addr();
   if (sHdrs)
   {
      File::mem().allot(hd,sSize);
      _hdr.sPtr() = hd.addr();
   }
   else
   {
      _hdr.sPtr() = fNullPtr;
   }
   _hdr.wr() = 0;
   _hdr.gSize() = gSize;
   _hdr.sSize() = sSize;
   _hdr.mSize() = mSize;
   _hdr.cSize() = cSize;
   _hdr.eData() = GPair::initialize(File::mem(),gSize,sSize,mSize,cSize);
   _gNames.resize(gSize);
   _sNames.resize(sSize);
   _cNames.resize(cSize);
   File::mem().sync(_hdr,FSync::write);
}



int CMatrix::gSize() const
{
   return _hdr.gSize();
}



int CMatrix::sSize() const
{
   return _hdr.sSize();
}



int CMatrix::mSize() const
{
   return _hdr.mSize();
}



int CMatrix::cSize() const
{
   return _hdr.cSize();
}



CMatrix::string& CMatrix::gName(int i)
{
   return get_name(i,_gNames,_hdr.gPtr(),_hdr.gSize());
}



CMatrix::string& CMatrix::sName(int i)
{
   return get_name(i,_sNames,_hdr.sPtr(),_hdr.sSize());
}



CMatrix::string& CMatrix::cName(int i)
{
   return get_name(i,_cNames,_hdr.cPtr(),_hdr.cSize());
}



void CMatrix::write()
{
   bool cond {_hdr.gPtr()!=fNullPtr};
   AccelCompEng::assert<NotInitialized>(cond,__LINE__);
   cond = _hdr.wr()==0;
   AccelCompEng::assert<AlreadySet>(cond,__LINE__);
   _hdr.wr() = 1;
   write_names(_gNames,_hdr.gPtr());
   write_names(_cNames,_hdr.cPtr());
   if (_hdr.sPtr()!=fNullPtr)
   {
      write_names(_sNames,_hdr.sPtr());
   }
   File::mem().sync(_hdr,FSync::write);
}



CMatrix::GPair CMatrix::begin()
{
   return GPair(this,1,0);
}



CMatrix::GPair CMatrix::end()
{
   return GPair(this,_hdr.gSize(),0);
}



CMatrix::GPair& CMatrix::at(int x, int y)
{
   bool cond {_hdr.gPtr()!=fNullPtr};
   AccelCompEng::assert<NotInitialized>(cond,__LINE__);
   cond = x>=0&&y>=0&&y<x&&x<_hdr.gSize();
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
   return ref(x,y);
}



CMatrix::GPair& CMatrix::ref(int x, int y)
{
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



CMatrix::string& CMatrix::get_name(int i, svec& names, FPtr ptr, int size)
{
   bool cond {ptr!=fNullPtr};
   AccelCompEng::assert<NotInitialized>(cond,__LINE__);
   cond = {i>=0||i<size};
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
   return names[i];
}



void CMatrix::write_names(svec& names, FPtr ptr)
{
   NmHead hd {ptr};
   for (int i=0;i<names.size();++i)
   {
      FString tmp {&File::mem()};
      tmp = names[i];
      hd.nPtr() = tmp.addr();
      File::mem().sync(hd,FSync::write,i);
   }
}



long long CMatrix::diagonal(int x, int y)
{
   return (x*(x-1)/2)+y;
}



long long CMatrix::diag_size(int x)
{
   return x*(x+1)/2;
}



CMatrix::GPair::~GPair()
{
   if (_iModes)
   {
      delete _iModes;
   }
   if (_iCorrs)
   {
      delete _iCorrs;
   }
}



void CMatrix::GPair::read()
{
   _p->File::mem().sync(_data,FSync::read,diagonal(_x,_y));
}



void CMatrix::GPair::write()
{
   _p->File::mem().sync(_data,FSync::write,diagonal(_x,_y));
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
   return _data.mAmt();
}



void CMatrix::GPair::size(int size)
{
   bool cond {size>=0&&size<=(_p->_hdr.mSize())};
   AccelCompEng::assert<GreaterThanMax>(cond,__LINE__);
   _data.mAmt() = size;
}



CMatrix::GPair::Modes& CMatrix::GPair::modes()
{
   if (!_iModes)
   {
      _iModes = new Modes(this);
   }
   return *_iModes;
}



CMatrix::GPair::Corrs& CMatrix::GPair::corrs()
{
   if (!_iCorrs)
   {
      _iCorrs = new Corrs(this);
   }
   return *_iCorrs;
}



void CMatrix::GPair::operator++()
{
   if (_x<(_p->_hdr.gSize()))
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
   _p(p),
   _data(p->_hdr.mSize(),p->_hdr.sSize(),p->_hdr.cSize(),p->_hdr.eData()),
   _x(x),
   _y(y)
{}



void CMatrix::GPair::set(int x, int y)
{
   _x = x;
   _y = y;
}



CMatrix::FPtr CMatrix::GPair::initialize(FileMem& mem, int gSize, int sSize, int mSize, int cSize)
{
   Pair ndat(mSize,sSize,cSize);
   mem.allot(ndat,diag_size(gSize));
   return ndat.addr();
}



CMatrix::GPair::Modes::~Modes()
{
   if (_iMode)
   {
      delete _iMode;
   }
}



CMatrix::GPair::Modes::Mode CMatrix::GPair::Modes::begin()
{
   return Mode(this,0);
}



CMatrix::GPair::Modes::Mode CMatrix::GPair::Modes::end()
{
   return Mode(this,_p->_data.mAmt());
}



CMatrix::GPair::Modes::Mode& CMatrix::GPair::Modes::at(int i)
{
   bool cond {i>=0&&i<(_p->_data.mAmt())};
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
   return (*this)[i];
}



CMatrix::GPair::Modes::Mode& CMatrix::GPair::Modes::operator[](int i)
{
   if (!_iMode)
   {
      _iMode = new Mode(this,i);
   }
   else
   {
      _iMode->set(i);
   }
   return *_iMode;
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
   return Iterator(this,_p->_p->_p->_hdr.sSize());
}



int8_t& CMatrix::GPair::Modes::Mode::at(int i)
{
   bool cond {i>=0&&i<(_p->_p->_p->_hdr.sSize())};
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
   return (*this)[i];
}



int8_t& CMatrix::GPair::Modes::Mode::operator[](int i)
{
   return _p->_p->_data.mVal(_i,i);
}



void CMatrix::GPair::Modes::Mode::operator++()
{
   if (_i<(_p->_p->_data.mAmt()))
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
   return _p->_p->_p->_data.mVal(_p->_i,_i);
}



void CMatrix::GPair::Modes::Mode::Iterator::operator++()
{
   if (_i<(_p->_p->_p->_p->_hdr.sSize()))
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



CMatrix::GPair::Corrs::~Corrs()
{
   if (_iCorr)
   {
      delete _iCorr;
   }
}



CMatrix::GPair::Corrs::Corr CMatrix::GPair::Corrs::begin()
{
   return Corr(this,0);
}



CMatrix::GPair::Corrs::Corr CMatrix::GPair::Corrs::end()
{
   return Corr(this,_p->_data.mAmt());
}



CMatrix::GPair::Corrs::Corr& CMatrix::GPair::Corrs::at(int i)
{
   bool cond {i>=0&&i<(_p->_data.mAmt())};
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
   return (*this)[i];
}



CMatrix::GPair::Corrs::Corr& CMatrix::GPair::Corrs::operator[](int i)
{
   if (!_iCorr)
   {
      _iCorr = new Corr(this,i);
   }
   else
   {
      _iCorr->set(i);
   }
   return *_iCorr;
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
   return Iterator(this,_p->_p->_p->_hdr.cSize());
}



float& CMatrix::GPair::Corrs::Corr::at(int i)
{
   bool cond {i>=0&&i<(_p->_p->_p->_hdr.cSize())};
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
   return (*this)[i];
}



float& CMatrix::GPair::Corrs::Corr::operator[](int i)
{
   return _p->_p->_data.cVal(_i,i);
}



void CMatrix::GPair::Corrs::Corr::operator++()
{
   if (_i<(_p->_p->_data.mAmt()))
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
   return _p->_p->_p->_data.cVal(_p->_i,_i);
}



void CMatrix::GPair::Corrs::Corr::Iterator::operator++()
{
   if (_i<(_p->_p->_p->_p->_hdr.cSize()))
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






///////////////////////////////////////////////////////////////////////////////////
/*


void CMatrix::load(GetOpts&, Terminal &tm)
{
   tm << "Not yet implemented.\n";
}



void CMatrix::dump(GetOpts&, Terminal &tm)
{
   tm << "Not yet implemented.\n";
}



void CMatrix::query(GetOpts&, Terminal &tm)
{
   tm << "Not yet implemented.\n";
}



bool CMatrix::empty()
{
   return _hdr.gSize()==0;
}



void CMatrix::initialize(uint32_t gSize, uint32_t sSize, uint32_t cSize,
                         uint8_t mdMax)
{
   bool cond {_hdr.gSize()==0};
   AccelCompEng::assert<AlreadyInitialized>(cond,__LINE__);
   cond = gSize>0&&sSize>0&&cSize>0&&mdMax>0;
   AccelCompEng::assert<InvalidSize>(cond,__LINE__);
   _hdr.gSize() = gSize;
   _hdr.sSize() = sSize;
   _hdr.cSize() = cSize;
   _hdr.mdMax() = mdMax;
   NmHead nh;
   File::mem().allot(nh,gSize);
   _hdr.gPtr() = nh.addr();
   File::mem().allot(nh,sSize);
   _hdr.sPtr() = nh.addr();
   File::mem().allot(nh,cSize);
   _hdr.cPtr() = nh.addr();
   _hdr.mdData() = GeneModes::initialize(*this,gSize,sSize,mdMax);
   _hdr.crData() = GeneCorrs::initialize(*this,gSize,cSize,mdMax);
   File::mem().sync(_hdr,FSync::write);
}



void CMatrix::set_gene_name(uint32_t i, const string& name)
{
   set_name(_hdr.gPtr(),i,_hdr.gSize(),name);
}



void CMatrix::set_sample_name(uint32_t i, const string& name)
{
   set_name(_hdr.sPtr(),i,_hdr.sSize(),name);
}



void CMatrix::set_correlation_name(uint32_t i, const string& name)
{
   set_name(_hdr.cPtr(),i,_hdr.cSize(),name);
}



CMatrix::GeneModes CMatrix::get_modes(uint32_t g1, uint32_t g2)
{
   bool cond {_hdr.mdData()!=fNullPtr};
   AccelCompEng::assert<NotInitialized>(cond,__LINE__);
   return GeneModes(this,_hdr.mdData()+diagonal(g1,g2));
}



void CMatrix::get_modes(GeneModes& gm, uint32_t g1, uint32_t g2)
{
   bool cond {_hdr.mdData()!=fNullPtr};
   AccelCompEng::assert<NotInitialized>(cond,__LINE__);
   gm.addr(_hdr.mdData()+diagonal(g1,g2));
}



CMatrix::GeneCorrs CMatrix::get_corrs(uint32_t g1, uint32_t g2)
{
   bool cond {_hdr.crData()!=fNullPtr};
   AccelCompEng::assert<NotInitialized>(cond,__LINE__);
   return GeneCorrs(this,_hdr.crData()+diagonal(g1,g2));
}



void CMatrix::get_corrs(GeneCorrs& gc, uint32_t g1, uint32_t g2)
{
   bool cond {_hdr.crData()!=fNullPtr};
   AccelCompEng::assert<NotInitialized>(cond,__LINE__);
   gc.addr(_hdr.crData()+diagonal(g1,g2));
}



void CMatrix::set_name(FPtr ptr, uint32_t i, uint32_t size, const string& name)
{
   bool cond {size>0};
   AccelCompEng::assert<NotInitialized>(cond,__LINE__);
   cond = i<size;
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
   NmHead nm {ptr};
   FString str {&File::mem()};
   str = name;
   nm.nPtr() = str.addr();
   File::mem().sync(nm,FSync::write,i);
}



CMatrix::SizeT CMatrix::diagonal(uint32_t g1, uint32_t g2)
{
   bool cond {g1!=g2};
   AccelCompEng::assert<InvalidGeneCorr>(cond,__LINE__);
   cond = g1<_hdr.gSize()&&g2<_hdr.gSize();
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
   if (g1<g2)
   {
      uint32_t tmp {g1};
      g1 = g2;
      g2 = tmp;
   }
   return (g1*(g1+1)/2)+g2;
}



CMatrix::SizeT CMatrix::diag_size(uint32_t gNum)
{
   return gNum*(gNum-1)/2;
}



bool CMatrix::GeneModes::mode(uint8_t md, uint32_t i)
{
   bool cond {md<_p->_hdr.mdMax()&&i<md<_p->_hdr.sSize()};
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
   AccelCompEng::assert<ValueNotRead>(_isRead,__LINE__);
   return _data.mask(md,i);
}



void CMatrix::GeneModes::mode(uint8_t md, uint32_t i, bool val)
{
   bool cond {md<_p->_hdr.mdMax()&&i<md<_p->_hdr.sSize()};
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
   _data.mask(md,i) = val;
}



void CMatrix::GeneModes::mode(uint8_t md, bool val)
{
   bool cond {md<_p->_hdr.mdMax()};
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
   for (int i=0;i<_p->_hdr.sSize();++i)
   {
      _data.mask(md,i) = val;
   }
}



void CMatrix::GeneModes::mode(bool val)
{
   for (int md=0;md<_p->_hdr.mdMax();++md)
   {
      for (int i=0;i<_p->_hdr.sSize();++i)
      {
         _data.mask(md,i) = val;
      }
   }
}



void CMatrix::GeneModes::read()
{
   _p->File::mem().sync(_data,FSync::read);
   _isRead = true;
}



void CMatrix::GeneModes::write()
{
   _p->File::mem().sync(_data,FSync::write);
}



CMatrix::GeneModes::GeneModes(CMatrix* p, FPtr ptr):
   _p(p),
   _data(p->_hdr.mdMax(),p->_hdr.sSize(),ptr)
{}



void CMatrix::GeneModes::addr(FPtr ptr)
{
   _data.addr(ptr);
   _isRead = false;
}



CMatrix::FPtr CMatrix::GeneModes::initialize(CMatrix& p, uint32_t gSize,
                                             uint32_t sSize,
                                             uint8_t mdMax)
{
   Modes md(mdMax,sSize);
   p.File::mem().allot(md,p.diag_size(gSize));
   return md.addr();
}



float CMatrix::GeneCorrs::corr(uint8_t md, uint32_t i)
{
   bool cond {md<_p->_hdr.mdMax()&&i<md<_p->_hdr.cSize()};
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
   AccelCompEng::assert<ValueNotRead>(_isRead,__LINE__);
   return _data.val(md,i);
}



void CMatrix::GeneCorrs::corr(uint8_t md, uint32_t i, float val)
{
   bool cond {md<_p->_hdr.mdMax()&&i<md<_p->_hdr.cSize()};
   AccelCompEng::assert<OutOfRange>(cond,__LINE__);
   _data.val(md,i) = val;
}



void CMatrix::GeneCorrs::read()
{
   _p->File::mem().sync(_data,FSync::read);
   _isRead = true;
}



void CMatrix::GeneCorrs::write()
{
   _p->File::mem().sync(_data,FSync::write);
}



CMatrix::GeneCorrs::GeneCorrs(CMatrix* p, FPtr ptr):
   _p(p),
   _data(p->_hdr.mdMax(),p->_hdr.cSize(),ptr)
{}



void CMatrix::GeneCorrs::addr(FPtr ptr)
{
   _data.addr(ptr);
   _isRead = false;
}



CMatrix::FPtr CMatrix::GeneCorrs::initialize(CMatrix& p, uint32_t gSize,
                                             uint32_t cSize,
                                             uint8_t mdMax)
{
   Corrs cr(mdMax,cSize);
   p.File::mem().allot(cr,p.diag_size(gSize));
   return cr.addr();
}
*/
