#include "spearman.h"
#include "spearman.cl.h"
#include "spearman2.cl.h"
#include <gsl/gsl_statistics.h>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <chrono>
#include <unistd.h>


/**
 * @brief Implements ACE's Analytic::input.
 */
void Spearman::input(Ace::Data* input)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<TooManyInputs>(!_in,f,__LINE__);
   Ace::assert<InvalidDataType>(input->type()==std::string("emx"),f,__LINE__);
   _in = dynamic_cast<EMatrix*>(input);
}


/**
 * @brief Implements ACE's Analytic::output.
 */
void Spearman::output(Ace::Data* output)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<TooManyOutputs>(!_out,f,__LINE__);
   Ace::assert<InvalidDataType>(output->type()==std::string("cmx"),f,__LINE__);
   _out = dynamic_cast<CMatrix*>(output);
}



void Spearman::theMaddening()
{
   // Compile opencl kernel code and load into memory.
   std::cout << "Loading kernel program into OpenCL device...\n";
   CLProgram::add_source(spearman2_cl);
   if (!CLProgram::compile(""))
   {
      std::cout << CLProgram::log();
      return;
   }
   auto kern = CLProgram::mkernel("calculateSpearmanBlock");
   auto expList = CLContext::buffer<cl_float>(28);
   int i=0;
   expList[i++] = 106;
   expList[i++] = 86;
   expList[i++] = 100;
   expList[i++] = NAN;
   expList[i++] = 101;
   expList[i++] = 99;
   expList[i++] = 103;
   expList[i++] = 97;
   expList[i++] = 113;
   expList[i++] = 1000000;
   expList[i++] = 112;
   expList[i++] = 1000000;
   expList[i++] = 1000000;
   expList[i++] = 110;
   //
   expList[i++] = 7;
   expList[i++] = 0;
   expList[i++] = 27;
   expList[i++] = 100000;
   expList[i++] = 50;
   expList[i++] = 28;
   expList[i++] = 29;
   expList[i++] = 20;
   expList[i++] = 12;
   expList[i++] = NAN;
   expList[i++] = 6;
   expList[i++] = NAN;
   expList[i++] = NAN;
   expList[i++] = 17;
   std::cout << "Loading kernel program into OpenCL device...\n";
   auto ev = CLCommandQueue::write_buffer(expList);
   ev.wait();
   kern.set_arg(0,14);
   kern.set_arg(1,16);
   kern.set_arg(2,1);
   kern.set_arg(4,&expList);
   auto buf1 = CLContext::buffer<cl_float>(2*16);
   kern.set_arg(5,&buf1);
   auto buf2 = CLContext::buffer<cl_int>(16);
   kern.set_arg(6,&buf2);
   auto ld = CLContext::buffer<int>(2);
   ld[0] = 0;
   ld[1] = 14;
   ev = CLCommandQueue::write_buffer(ld);
   ev.wait();
   auto ans = CLContext::buffer<cl_float>(1);
   kern.set_arg(3,&ld);
   kern.set_arg(7,&ans);
   kern.set_swarm_dims(1);
   kern.set_swarm_size(0,1,1);
   std::cout << "Loading kernel program into OpenCL device...\n";
   ev = CLCommandQueue::add_swarm(kern);
   ev.wait();
   ev = CLCommandQueue::read_buffer(buf1);
   ev.wait();
   ev = CLCommandQueue::read_buffer(buf2);
   ev.wait();
   for (int i = 0; i < 16 ;++i)
   {
      std::cout << buf1[i] << " ";
   }
   std::cout << "\n";
   for (int i = 16; i < 32 ;++i)
   {
      std::cout << buf1[i] << " ";
   }
   std::cout << "\n";
   for (int i = 0; i < 16 ;++i)
   {
      std::cout << buf2[i] << " ";
   }
   std::cout << "\n";
   std::cout << "Loading kernel program into OpenCL device...\n";
   ev = CLCommandQueue::read_buffer(ans);
   ev.wait();
   std::cout << ans[0] << "\n";
}



void Spearman::theMaddening2()
{
   // Compile opencl kernel code and load into memory.
   std::cout << "Loading kernel program into OpenCL device...\n";
   CLProgram::add_source(spearman2_cl);
   if (!CLProgram::compile(""))
   {
      std::cout << CLProgram::log();
      return;
   }
   auto kern = CLProgram::mkernel("calculateSpearmanBlock");
   auto expList = CLContext::buffer<cl_float>(40);
   int i=0;
   expList[i++] = 106;
   expList[i++] = 86;
   expList[i++] = 100;
   expList[i++] = 101;
   expList[i++] = 99;
   expList[i++] = 103;
   expList[i++] = 97;
   expList[i++] = 113;
   expList[i++] = 112;
   expList[i++] = 110;
   expList[i++] = 7;
   expList[i++] = 0;
   expList[i++] = 27;
   expList[i++] = 50;
   expList[i++] = 28;
   expList[i++] = 29;
   expList[i++] = 20;
   expList[i++] = 12;
   expList[i++] = 6;
   expList[i++] = 17;
   expList[i++] = 7;
   expList[i++] = 0;
   expList[i++] = 27;
   expList[i++] = 50;
   expList[i++] = 28;
   expList[i++] = 29;
   expList[i++] = 20;
   expList[i++] = 12;
   expList[i++] = 6;
   expList[i++] = 17;
   std::cout << "Loading kernel program into OpenCL device...\n";
   auto ev = CLCommandQueue::write_buffer(expList);
   ev.wait();
   kern.set_arg(0,10);
   kern.set_arg(1,16);
   kern.set_arg(2,1);
   kern.set_arg(4,&expList);
   auto buf1 = CLContext::buffer<cl_float>(2*16*4);
   kern.set_arg(5,&buf1);
   auto buf2 = CLContext::buffer<cl_int>(16*4);
   kern.set_arg(6,&buf2);
   auto ld = CLContext::buffer<int>(2*4);
   ld[0] = 10;
   ld[1] = 0;
   ld[2] = 20;
   ld[3] = 0;
   ld[4] = 20;
   ld[5] = 10;
   ld[6] = 0;
   ld[7] = 0;
   ev = CLCommandQueue::write_buffer(ld);
   ev.wait();
   auto ans = CLContext::buffer<cl_float>(4);
   kern.set_arg(3,&ld);
   kern.set_arg(7,&ans);
   kern.set_swarm_dims(1);
   kern.set_swarm_size(0,4,4);
   std::cout << "Loading kernel program into OpenCL device...\n";
   ev = CLCommandQueue::add_swarm(kern);
   ev.wait();
   ev = CLCommandQueue::read_buffer(buf1);
   ev.wait();
   ev = CLCommandQueue::read_buffer(buf2);
   ev.wait();
   for (int i = 0+64; i < 16+64 ;++i)
   {
      std::cout << buf1[i] << " ";
   }
   std::cout << "\n";
   for (int i = 16+64; i < 32+64 ;++i)
   {
      std::cout << buf1[i] << " ";
   }
   std::cout << "\n";
   for (int i = 0+32; i < 16+32 ;++i)
   {
      std::cout << buf2[i] << " ";
   }
   std::cout << "\n";
   std::cout << "Loading kernel program into OpenCL device...\n";
   ev = CLCommandQueue::read_buffer(ans);
   ev.wait();
   for (int i = 0; i < 3 ;++i)
   {
      std::cout << ans[i] << " ";
   }
   std::cout << "\n";
}



/**
 * @brief Implements ACE's Analytic::execute_cl.
 *
 * @param ops User arguments ACE object.
 * @param tm ACE terminal object.
 *
 */
void Spearman::execute_cl(Ace::GetOpts& ops, Ace::Terminal& tm)
{
   //theMaddening();
   //theMaddening2();
   // Get all user arguments.
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
      else if (i.is_key("minsize"))
      {
         minSize = i.value<int>();
      }
   }
   // Compile opencl kernel code and load into memory.
   tm << "Loading kernel program into OpenCL device...\n";
   CLProgram::add_source(spearman2_cl);
   if (!CLProgram::compile(""))
   {
      tm << CLProgram::log();
      return;
   }
   auto kern = CLProgram::mkernel("calculateSpearmanBlock");
   Ace::assert<NoDataInput>(_in,f,__LINE__);
   Ace::assert<NoDataOutput>(_out,f,__LINE__);
   int gSize {_in->gene_size()};
   int sSize {_in->sample_size()};
   // Load expression data into memory and then transfer to opencl device global memory.
   tm << "Loading expression data into OpenCL device...\n";
   auto expList = CLContext::buffer<cl_float>(sSize*gSize);
   int inc {0};
   for (auto g = _in->begin();g!=_in->end();++g)
   {
      g.read();
      for (auto i = g.begin();i!=g.end();++i)
      {
         expList[inc++] = *i;
      }
   }
   auto ev = CLCommandQueue::write_buffer(expList);
   ev.wait();
   // Pre-populate the header of the output correlation matrix data object.
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
   // Begin OpenCL processing engine, this is where all the magic happens.
   calculate(tm,kern,expList,sSize,blSize,smSize,minSize);
   // Calculates the time it took for total execution.
   auto t2 = system_clock::now();
   int s = duration_cast<seconds>(t2-t1).count();
   if (s==0)
   {
      tm << "Finished in less than 1 second.\n";
   }
   else
   {
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
}

/**
 * @brief Implements ACE's Analytic::execute_pn.
 *
 * @param ops User arguments ACE object.
 * @param tm ACE terminal object.
 *
 */
void Spearman::execute_pn(Ace::GetOpts& ops, Ace::Terminal& tm)
{
   static const char* f = __PRETTY_FUNCTION__;
   using namespace std::chrono;
   auto t1 = system_clock::now();
   int minSize {30};
   Ace::assert<NoDataInput>(_in,f,__LINE__);
   Ace::assert<NoDataOutput>(_out,f,__LINE__);
   int gSize {_in->gene_size()};
   int sSize {_in->sample_size()};
   // Pre-populate header information of correlation matrix object.
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
   // Begin huge loop to calculate all spearman coefficients.
   _out->initialize(std::move(geneNames),std::move(sampleNames),std::move(correlations),1);
   tm << "Calculating spearman (pn) values and saving to output file[0%]...";
   auto i = _out->begin();
   i.size(1);
   auto bmode = i.modes().begin();
   for (auto m = bmode.begin();m!=bmode.end();++m)
   {
      *m = 1;
   }
   unsigned long total  = CMatrix::diag_size(gSize);
   unsigned long count {0};
   double a[sSize];
   double b[sSize];
   double work[2*sSize];
   int delay {0};
   for (;i!=_out->end();++i)
   {
      int size {0};
      auto gX = _in->find(i.x());
      auto gY = _in->find(i.y());
      gX.read();
      gY.read();
      for (auto x = 0;x<sSize;++x)
      {
         if ( !std::isnan(gX[x]) && !std::isnan(gY[x]) )
         {
            a[size] = gX[x];
            b[size] = gY[x];
            ++size;
         }
      }
      if (size<minSize)
      {
         auto t = i.corrs().find(0);
         t.at(0) = NAN;
      }
      else
      {
         auto t = i.corrs().find(0);
         t.at(0) = gsl_stats_spearman(a,1,b,1,size,work);
      }
      ++count;
      ++delay;
      if (delay==16384)
      {
         tm << "\rCalculating spearman values and saving to output file["
            << 100*count/total << "%]..." << Ace::Terminal::flush;
         delay = 0;
      }
   }
   tm << "\rCalculating spearman values and saving to output file[" << 100*count/total << "%]...\n";
   // Calculate time taken to run all gene pairs.
   auto t2 = system_clock::now();
   int s = duration_cast<seconds>(t2-t1).count();
   if (s==0)
   {
      tm << "Finished in less than 1 second.\n";
   }
   else
   {
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
}


/**
 * @brief Calculates the nearest power of 2  above the value provided.
 *
 * @param i
 *   The value for finding the newarest power of 2.
 */
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


/**
 * @brief Calculates the nearest power of 2  below the value provided.
 *
 * @param i
 *   The value for finding the newarest power of 2.
 */
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
                         int blSize, int smSize, int minSize)
{
   // Initialize the kernel with all arguments that never change.
   enum class State {start,in,exec,out,end};
   unsigned long total = CMatrix::diag_size(_in->gene_size());
   int workSize {pow2_ceil(size)};
   kern.set_arg(0,size);
   kern.set_arg(1,workSize);
   kern.set_arg(2,minSize);
   kern.set_arg(4,&expList);
   auto buf1 = CLContext::buffer<cl_float>(2*workSize*blSize);
   kern.set_arg(5,&buf1);
   auto buf2 = CLContext::buffer<cl_int>(workSize*blSize);
   kern.set_arg(6,&buf2);

   kern.set_swarm_dims(1);
   //kern.set_swarm_size(0,blSize,1024);
   // Define a structure that defines a single kernel execution then make array of them that equal
   // smSize.
   struct
   {
      State st;
      int x;
      int y;
      int size;
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
   // Begin huge loop that continues running until all kernel execution pathways are no longer
   // alive.
   int alive {smSize};
   int si {0};
   auto i = _out->begin();
   unsigned long done {0};
   while (alive>0)
   {
      switch (state[si].st)
      {
      case State::start:
         // Kernel is in START state, load new array targets to the kernel. If there are
         // no gene pairs left to calculate move state to END and decrement alive count.
         if (i!=_out->end())
         {
            state[si].x = i.x();
            state[si].y = i.y();
            int count {0};
            while (i!=_out->end()&&count<blSize)
            {
               state[si].ld[2*count] = i.x()*size;
               state[si].ld[(2*count)+1] = i.y()*size;
               ++count;
               ++i;
            }
            state[si].size = count;
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
         // Kernel is in IN state, once targets have been loaded begin execution of actual
         // kernel.
         if (state[si].ev.is_done())
         {
            int size {pow2_ceil(state[si].size)};
            int wgsize {1};
            if ( size > 1 )
            {
               wgsize = 2;
               while ( wgsize < size && wgsize < 1024 )
               {
                  wgsize *= 2;
               }
            }
            kern.set_swarm_size(0,size,wgsize);
            kern.set_arg(3,&(state[si].ld));
            kern.set_arg(7,&(state[si].ans));
            state[si].ev = CLCommandQueue::add_swarm(kern);
            state[si].st = State::exec;
         }
         break;
      case State::exec:
         // Kernel is in EXEC state, once kernel is done begin loading results back to system
         // memory.
         if (state[si].ev.is_done())
         {
            state[si].ev = CLCommandQueue::read_buffer(state[si].ans);
            state[si].st = State::out;
         }
         break;
      case State::out:
         // Kernel is in OUT state, once kernel is done loading results back system read those
         // results and load them into new correlation matrix object. Lastly, move back to the
         // START state.
         if (state[si].ev.is_done())
         {
            auto wi = _out->find(state[si].x,state[si].y);
            wi.size(1);
            auto bmode = wi.modes().begin();
            for (auto m = bmode.begin();m!=bmode.end();++m)
            {
               *m = 1;
            }
            int count {0};
            while (wi!=_out->end()&&count<blSize)
            {
               auto bcorr = wi.corrs().begin();
               bcorr[0] = state[si].ans[count];
               ++count;
               wi.write();
               ++wi;
               ++done;
            }
            state[si].st = State::start;
         }
         break;
      case State::end:
         // Kernel is in END state, do nothing like zombie.
         break;
      }
      // Increment to next kernel branch, if reached end of array go back to beginning.
      ++si;
      if (si>=smSize)
      {
         si = 0;
      }
      tm << "\rCalculating spearman values and saving to output file[" << 100*done/total
         << "%]..." << Ace::Terminal::flush;
      usleep(100);
   }
   tm << "\n";
}



#ifdef UNIT_TEST



void Spearman::runUnitTests(Ace::Console& console)
{
   console.command("cl set 0:0");
   console.command("open simple:emx --select");
   console.command("load samples.txt");
   console.command("spearman --in=simple --out=simple_spearman:cmx --minsize=0");
   console.command("select simple_spearman");
   console.command("query");
}



#endif
