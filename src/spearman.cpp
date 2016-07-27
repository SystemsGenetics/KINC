#include "spearman.h"
#include "spearman.cl.h"
#include <gsl/gsl_statistics.h>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <chrono>



void Spearman::input(DataPlugin* input)
{
   AccelCompEng::assert<TooManyInputs>(!_in,__LINE__);
   bool cond {input->type()==string("emx")};
   AccelCompEng::assert<InvalidDataType>(cond,__LINE__);
   _in = dynamic_cast<EMatrix*>(input);
}



void Spearman::output(DataPlugin* output)
{
   AccelCompEng::assert<TooManyOutputs>(!_out,__LINE__);
   bool cond {output->type()==string("cmx")};
   AccelCompEng::assert<InvalidDataType>(cond,__LINE__);
   _out = dynamic_cast<CMatrix*>(output);
}



void Spearman::execute_cl(GetOpts& ops, Terminal& tm)
{
   using namespace std::chrono;
   auto t1 = system_clock::now();
   int blSize {8192};
   int smSize {4};
   for (auto i = ops.begin();i!=ops.end();++i)
   {
      if (i.is_key("slots"))
      {
         smSize = i.value<int>();
      }
      else if (i.is_key("bsize"))
      {
         blSize = i.value<int>();
      }
   }
   tm << "Loading kernel program into OpenCL device...\n";
   CLProgram::add_source(spearman_cl);
   if (!CLProgram::compile(""))
   {
      tm << CLProgram::log();
      return;
   }
   auto kern = CLProgram::mkernel("spearman");
   AccelCompEng::assert<NoDataInput>(_in,__LINE__);
   AccelCompEng::assert<NoDataOutput>(_out,__LINE__);
   int gSize {_in->gSize()};
   int sSize {_in->sSize()};
   tm << "Loading expression data into OpenCL device...\n";
   auto expList = CLContext::buffer<cl_float>(sSize*gSize);
   int inc {0};
   for (auto g = _in->begin();g!=_in->end();++g)
   {
      for (auto i = g.begin();i!=g.end();++i)
      {
         expList[inc++] = *i;
      }
   }
   CLEvent ev = CLCommandQueue::write_buffer(expList);
   ev.wait();
   tm << "Formating and copying header information to output file...\n";
   _out->initialize(gSize,sSize,1,1,_in->hasSampleHead());
   for (int i = 0;i<_in->gSize();++i)
   {
      _out->gName(i) = _in->gName(i);
   }
   if (_in->hasSampleHead())
   {
      for (int i = 0;i<_in->sSize();++i)
      {
         _out->sName(i) =
               _in->sName(i);
      }
   }
   _out->sName(0) = "spearman";
   _out->write();
   tm << "Calculating spearman values and saving to output file[0%]...";
   int bSize {pow2_ceil(sSize)};
   bool cond {bSize>0};
   AccelCompEng::assert<TooManySamples>(cond,__LINE__);
   int wSize {pow2_floor(kern.get_wg_size())};
   int chunk {1};
   if ((2*wSize)<bSize)
   {
      chunk = bSize/(2*wSize);
   }
   calculate(tm,kern,expList,sSize,wSize,chunk,blSize,smSize);
   auto t2 = system_clock::now();
   int s = duration_cast<seconds>(t2-t1).count();
   int m = s/60;
   int h = m/60;
   s %= 60;
   tm << "Finished in";
   if (h>0)
   {
      tm << " " << h << " hour(s)";
   }
   if (m>0)
   {
      tm << " " << m << " minute(s)";
   }
   if (s>0)
   {
      tm << " " << s << " second(s)";
   }
   tm << ".\n";
}



void Spearman::execute_pn(GetOpts&, Terminal& tm)
{
   using namespace std::chrono;
   auto t1 = system_clock::now();
   AccelCompEng::assert<NoDataInput>(_in,__LINE__);
   AccelCompEng::assert<NoDataOutput>(_out,__LINE__);
   int gSize {_in->gSize()};
   int sSize {_in->sSize()};
   tm << "Formating and copying header information to output file...\n";
   _out->initialize(gSize,sSize,1,1,_in->hasSampleHead());
   for (int i = 0;i<_in->gSize();++i)
   {
      _out->gName(i) = _in->gName(i);
   }
   if (_in->hasSampleHead())
   {
      for (int i = 0;i<_in->sSize();++i)
      {
         _out->sName(i) =
               _in->sName(i);
      }
   }
   _out->sName(0) = "spearman";
   _out->write();
   tm << "Calculating spearman values and saving to output file[0%]...";
   auto i = _out->begin();
   i.size(1);
   for (auto m = i.modes().at(0).begin();m!=i.modes().at(0).end();++m)
   {
      *m = 1;
   }
   double total  = CMatrix::diag_size(gSize);
   double count {0};
   double a[gSize];
   double b[gSize];
   double work[2*gSize];
   int delay {0};
   for (;i!=_out->end();++i)
   {
      for (int x = 0;x<sSize;++x)
      {
         a[x] = _in->at(i.x()).at(x);
         b[x] = _in->at(i.y()).at(x);
      }
      i.corrs().at(0).at(0) = gsl_stats_spearman(a,1,b,1,gSize,work);
      ++count;
      ++delay;
      if (delay==16384)
      {
         tm << "\rCalculating spearman values and saving to output file["
            << (int)((count/total*100)+0.5) << "%]..." << Terminal::flush;
         delay = 0;
      }
   }
   auto t2 = system_clock::now();
   int s = duration_cast<seconds>(t2-t1).count();
   int m = s/60;
   int h = m/60;
   s %= 60;
   tm << "Finished in";
   if (h>0)
   {
      tm << " " << h << " hour(s)";
   }
   if (m>0)
   {
      tm << " " << m << " minute(s)";
   }
   if (s>0)
   {
      tm << " " << s << " second(s)";
   }
   tm << ".\n";
}



int Spearman::pow2_ceil(int i)
{
   int ret = 1;
   while (ret<i)
   {
      ret<<=1;
      if (ret<=0)
      {
         break;
      }
   }
   return ret;
}



int Spearman::pow2_floor(int i)
{
   int ret = pow2_ceil(i);
   if (ret>i)
   {
      ret>>=1;
   }
   return ret;
}



void Spearman::calculate(Terminal& tm, CLKernel& kern, elist& expList, int size, int wSize,
                         int chunk, int blSize, int smSize)
{
   enum class State {start,in,exec,out,end};
   double total = CMatrix::diag_size(_in->gSize());
   kern.set_arg(0,size);
   kern.set_arg(1,chunk);
   kern.set_arg(3,&expList);
   kern.set_arg(5,sizeof(cl_float)*wSize*chunk*2);
   kern.set_arg(6,sizeof(cl_float)*wSize*chunk*2);
   kern.set_arg(7,sizeof(cl_int)*wSize*chunk*2);
   kern.set_arg(8,sizeof(cl_long)*wSize*chunk*2);
   kern.set_swarm_dims(1);
   kern.set_swarm_size(0,blSize*wSize,wSize);
   struct
   {
      State st;
      int x;
      int y;
      CLEvent ev;
      AccelCompEng::CLBuffer<int> ld;
      AccelCompEng::CLBuffer<cl_float> ans;
   } state[smSize];
   for (int i = 0;i<smSize;++i)
   {
      state[i].st = State::start;
      state[i].x = 0;
      state[i].y = 0;
      state[i].ev = CLEvent();
      state[i].ld = CLContext::buffer<int>(2*blSize);
      state[i].ans = CLContext::buffer<cl_float>(blSize);
   }
   int alive {smSize};
   int si {0};
   auto i = _out->begin();
   double done {0};
   while (alive>0)
   {
      switch (state[si].st)
      {
      case State::start:
         if (i!=_out->end())
         {
            state[si].x = i.x();
            state[si].y = i.y();
            int count {0};
            while (i!=_out->end()&&count<blSize)
            {
               state[si].ld[2*count] = i.x();
               state[si].ld[(2*count)+1] = i.y();
               ++count;
               ++i;
            }
            while (count<blSize)
            {
               state[si].ld[2*count] = 0;
               state[si].ld[(2*count)+1] = 0;
               ++count;
            }
            state[si].ev = CLCommandQueue::write_buffer(state[si].ld);
            state[si].st = State::in;
         }
         else
         {
            --alive;
            state[si].st = State::end;
         }
         break;
      case State::in:
         if (state[si].ev.is_done())
         {
            kern.set_arg(2,&(state[si].ld));
            kern.set_arg(4,&(state[si].ans));
            state[si].ev = CLCommandQueue::add_swarm(kern);
            state[si].st = State::exec;
         }
         break;
      case State::exec:
         if (state[si].ev.is_done())
         {
            state[si].ev = CLCommandQueue::read_buffer(state[si].ans);
            state[si].st = State::out;
         }
         break;
      case State::out:
         if (state[si].ev.is_done())
         {
            auto wi = _out->ref(state[si].x,state[si].y);
            wi.size(1);
            for (auto m = wi.modes().at(0).begin();m!=wi.modes().at(0).end();++m)
            {
               *m = 1;
            }
            int count {0};
            while (wi!=_out->end()&&count<blSize)
            {
               wi.corrs().at(0).at(0) = state[si].ans[count];
               ++count;
               wi.write();
               ++wi;
               ++done;
            }
            state[si].st = State::start;
         }
         break;
      case State::end:
         break;
      }
      ++si;
      if (si>=smSize)
      {
         si = 0;
      }
      tm << "\rCalculating spearman values and saving to output file["
         << (int)((done/total*100)+0.5) << "%]..." << Terminal::flush;
      usleep(100);
   }
   tm << "\n";
}
