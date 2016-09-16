#include "spearman.h"
#include "spearman.cl.h"
#include <gsl/gsl_statistics.h>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <chrono>
#include <unistd.h>



void Spearman::input(Ace::Data* input)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<TooManyInputs>(!_in,f,__LINE__);
   Ace::assert<InvalidDataType>(input->type()==string("emx"),f,__LINE__);
   _in = dynamic_cast<EMatrix*>(input);
}



void Spearman::output(Ace::Data* output)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<TooManyOutputs>(!_out,f,__LINE__);
   Ace::assert<InvalidDataType>(output->type()==string("cmx"),f,__LINE__);
   _out = dynamic_cast<CMatrix*>(output);
}



void Spearman::execute_cl(Ace::GetOpts& ops, Ace::Terminal& tm)
{
   static const char* f = __PRETTY_FUNCTION__;
   using namespace std::chrono;
   auto t1 = system_clock::now();
   int blSize {8192};
   int smSize {4};
   int minSize {30};
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
   Ace::assert<NoDataInput>(_in,f,__LINE__);
   Ace::assert<NoDataOutput>(_out,f,__LINE__);
   int gSize {_in->gene_size()};
   int sSize {_in->sample_size()};
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
   auto ev = CLCommandQueue::write_buffer(expList);
   ev.wait();
   tm << "Formating and copying header information to output file...\n";
   std::vector<std::string> geneNames;
   for (int i = 0;i<_in->gene_size();++i)
   {
      geneNames.push_back(_in->gene_name(i));
   }
   std::vector<std::string> sampleNames;
   for (int i = 0;i<_in->sample_size();++i)
   {
      sampleNames.push_back(_in->sample_name(i));
   }
   std::vector<std::string> correlations {"spearman"};
   _out->initialize(std::move(geneNames),std::move(sampleNames),std::move(correlations),1);
   tm << "Calculating spearman values and saving to output file[0%]...";
   int bSize {pow2_ceil(sSize)};
   Ace::assert<TooManySamples>(bSize>0,f,__LINE__);
   int wSize {pow2_floor(kern.get_wg_size())};
   int chunk {1};
   if ((2*wSize)<bSize)
   {
      chunk = bSize/(2*wSize);
   }
   calculate(tm,kern,expList,sSize,wSize,chunk,blSize,smSize,minSize);
   auto t2 = system_clock::now();
   int s = duration_cast<seconds>(t2-t1).count();
   int m = s/60;
   int h = m/60;
   m %= 60;
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



void Spearman::execute_pn(Ace::GetOpts&, Ace::Terminal& tm)
{
   static const char* f = __PRETTY_FUNCTION__;
   using namespace std::chrono;
   auto t1 = system_clock::now();
   int minSize {30};
   Ace::assert<NoDataInput>(_in,f,__LINE__);
   Ace::assert<NoDataOutput>(_out,f,__LINE__);
   int gSize {_in->gene_size()};
   int sSize {_in->sample_size()};
   tm << "Formating and copying header information to output file...\n";
   std::vector<std::string> geneNames;
   for (int i = 0;i<_in->gene_size();++i)
   {
      geneNames.push_back(_in->gene_name(i));
   }
   std::vector<std::string> sampleNames;
   for (int i = 0;i<_in->sample_size();++i)
   {
      sampleNames.push_back(_in->sample_name(i));
   }
   std::vector<std::string> correlations {"spearman"};
   _out->initialize(std::move(geneNames),std::move(sampleNames),std::move(correlations),1);
   tm << "Calculating spearman values and saving to output file[0%]...";
   auto i = _out->begin();
   i.size(1);
   for (auto m = i.modes().at(0).begin();m!=i.modes().at(0).end();++m)
   {
      *m = 1;
   }
   double total  = CMatrix::diag_size(gSize);
   double count {0};
   double a[sSize];
   double b[sSize];
   double work[2*sSize];
   int delay {0};
   for (;i!=_out->end();++i)
   {
      int size {0};
      for (int x = 0;x<sSize;++x)
      {
         if (!std::isnan(_in->at(i.x()).at(x))&&!std::isnan(_in->at(i.y()).at(x)))
         {
            a[size] = _in->at(i.x()).at(x);
            b[size] = _in->at(i.y()).at(x);
            ++size;
         }
      }
      if (size<minSize)
      {
         i.corrs().at(0).at(0) = NAN;
      }
      else
      {
         i.corrs().at(0).at(0) = gsl_stats_spearman(a,1,b,1,size,work);
      }
      ++count;
      ++delay;
      if (delay==16384)
      {
         tm << "\rCalculating spearman values and saving to output file["
            << (int)((count/total*100)+0.5) << "%]..." << Ace::Terminal::flush;
         delay = 0;
      }
   }
   auto t2 = system_clock::now();
   int s = duration_cast<seconds>(t2-t1).count();
   int m = s/60;
   int h = m/60;
   m %= 60;
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



void Spearman::calculate(Ace::Terminal& tm, Ace::CLKernel& kern, elist& expList, int size,
                         int wSize, int chunk, int blSize, int smSize, int minSize)
{
   enum class State {start,in,exec,out,end};
   double total = CMatrix::diag_size(_in->gene_size());
   int bufferSize {wSize*chunk*2};
   kern.set_arg(0,size);
   kern.set_arg(1,chunk);
   kern.set_arg(2,minSize);
   kern.set_arg(4,&expList);
   kern.set_arg(6,sizeof(cl_float)*bufferSize);
   kern.set_arg(7,sizeof(cl_float)*bufferSize);
   kern.set_arg(8,sizeof(cl_int)*bufferSize);
   kern.set_arg(9,sizeof(cl_int)*bufferSize);
   kern.set_arg(10,sizeof(cl_long)*bufferSize);
   kern.set_arg(11,sizeof(cl_float)*bufferSize);
   kern.set_arg(12,sizeof(cl_float)*bufferSize);
   kern.set_arg(13,sizeof(cl_int)*bufferSize);
   kern.set_arg(14,sizeof(cl_int)*bufferSize);
   kern.set_arg(15,sizeof(cl_int)*bufferSize);
   kern.set_arg(16,sizeof(cl_int)*bufferSize);
   kern.set_swarm_dims(1);
   kern.set_swarm_size(0,blSize*wSize,wSize);
   struct
   {
      State st;
      int x;
      int y;
      Ace::CLEvent ev;
      AccelCompEng::CLBuffer<int> ld;
      AccelCompEng::CLBuffer<cl_float> ans;
   } state[smSize];
   for (int i = 0;i<smSize;++i)
   {
      state[i].st = State::start;
      state[i].x = 0;
      state[i].y = 0;
      state[i].ev = Ace::CLEvent();
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
         << (int)((done/total*100)+0.5) << "%]..." << Ace::Terminal::flush;
      usleep(100);
   }
   tm << "\n";
}
