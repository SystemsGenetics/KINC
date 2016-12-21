
/*****************************************************************************/
#include "pearson.h"
// generate this with mkcl once pearson.cl is done
#include "pearson.cl.h"
/*****************************************************************************/

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
void Pearson::input(Ace::Data* input)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<TooManyInputs>(!_in,f,__LINE__);
   Ace::assert<InvalidDataType>(input->type()==std::string("emx"),f,__LINE__);
   _in = dynamic_cast<EMatrix*>(input);
}


/**
 * @brief Implements ACE's Analytic::output.
 */
void Pearson::output(Ace::Data* output)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<TooManyOutputs>(!_out,f,__LINE__);
   Ace::assert<InvalidDataType>(output->type()==std::string("cmx"),f,__LINE__);
   _out = dynamic_cast<CMatrix*>(output);
}


/**
 * @brief
 *
 * @param ops
 * @param tm
 *
 */
void Pearson::execute_cl(Ace::GetOpts& ops, Ace::Terminal& tm)
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
      else if (i.is_key("minsize"))
      {
         minSize = i.value<int>();
      }
   }
   tm << "Loading kernel program into OpenCL device...\n";
   CLProgram::add_source(pearson_cl);
   if (!CLProgram::compile(""))
   {
      tm << CLProgram::log();
      return;
   }
   auto kern = CLProgram::mkernel("pearson");
   Ace::assert<NoDataInput>(_in,f,__LINE__);
   Ace::assert<NoDataOutput>(_out,f,__LINE__);
   int gSize {_in->gene_size()};
   int sSize {_in->sample_size()};
   tm << "Loading expression data into OpenCL device...\n";

   // A linear list of all the sample scores (sample # size * gene # size)
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
   std::vector<std::string> correlations {"pearson"};
   _out->initialize(std::move(geneNames),std::move(sampleNames),std::move(correlations),1);
   tm << "Calculating pearson values and saving to output file[0%]...";

   // bSize is the sample size rounded to the nearest higer power of 2.
   // TODO: but what is it?
   int bSize {pow2_ceil(sSize)};
   Ace::assert<TooManySamples>(bSize>0,f,__LINE__);
   // wSize is the work group size rounded to the nearest lower power of 2.
   int wSize {pow2_floor(kern.get_wg_size())};
   // TODO: what is the chunk size????
   int chunk {1};
   if ((2*wSize)<bSize)
   {
      chunk = bSize/(2*wSize);
   }
   calculate(tm,kern,expList,sSize,wSize,chunk,blSize,smSize,minSize);

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
 * @brief
 *
 * @param ops
 * @param tm
 *
 */
void Pearson::execute_pn(Ace::GetOpts& ops, Ace::Terminal& tm)
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
   std::vector<std::string> correlations {"pearson"};
   _out->initialize(std::move(geneNames),std::move(sampleNames),std::move(correlations),1);
   tm << "Calculating pearson values and saving to output file[0%]...";
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


         // TODO: find gsl_stats_spearman()
         // Only found in pn ... why?
         t.at(0) = gsl_stats_spearman(a,1,b,1,size,work);
      }
      ++count;
      ++delay;
      if (delay==16384)
      {
         tm << "\rCalculating pearson values and saving to output file["
            << 100*count/total << "%]..." << Ace::Terminal::flush;
         delay = 0;
      }
   }
   tm << "\rCalculating pearson values and saving to output file[" << 100*count/total << "%]...\n";
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
int Pearson::pow2_ceil(int i)
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
int Pearson::pow2_floor(int i)
{
   int ret = pow2_ceil(i);
   if (ret>i)
   {
      ret>>=1;
   }
   return ret;
}


/**
 * @brief Executes the pearson correlation algorithm on the GPU.
 *
 * Divides the data into chunks that can be analyzed one chuck at a time
 * by an OpenCL kernel.
 *
 * @param tm
 *   A pointer to the ACE terminal console.
 * @param kern
 *   A pointer to the CLProgram::mkernel object of "type" pearson.
 * @param expList
 *   A pointer to an CLContext::buffer that has been pre populated with data
 *   from the expression matrix.
 * @param size
 *   The number of samples in the expression matrix.
 * @param wSize
 *   The work group size (i.e. the number of kernels a work group can hold).
 * @param chunk
 *  // TODO: what is the chunk size
 * @param smSize
 *  // TODO: what is the smSize
 * @param minSize
 *  // TODO: what is the minSize
 */

void Pearson::calculate(Ace::Terminal& tm, Ace::CLKernel& kern,
                        /* the giant linear list*/elist& expList, int size,
                        int wSize, int chunk, /*What is blsize?*/int blSize,
                        int smSize, int minSize)
{
   enum class State {start,in,exec,out,end};
   unsigned long total = CMatrix::diag_size(_in->gene_size());
   int bufferSize {wSize*chunk*2};
   kern.set_arg(0,size);
   kern.set_arg(1,chunk);
   kern.set_arg(2,minSize);
   kern.set_arg(4,&expList);
   auto buf1 = CLContext::buffer<cl_float>(blSize*bufferSize);
   kern.set_arg(6,&buf1);
   auto buf2 = CLContext::buffer<cl_float>(blSize*bufferSize);

   //kern.set_arg(7,&buf2);
   //auto buf3 = CLContext::buffer<cl_int>(blSize*bufferSize);
   //kern.set_arg(8,&buf3);
   //auto buf4 = CLContext::buffer<cl_int>(blSize*bufferSize);
   //kern.set_arg(9,&buf4);
   //auto buf5 = CLContext::buffer<cl_long>(blSize*bufferSize);
   kern.set_arg(7,&buf5);
   auto buf6 = CLContext::buffer<cl_float>(blSize*bufferSize);
   kern.set_arg(8,&buf6);
   auto buf7 = CLContext::buffer<cl_float>(blSize*bufferSize);
   kern.set_arg(9,&buf7);
   auto buf8 = CLContext::buffer<cl_int>(blSize*bufferSize);
   kern.set_arg(10,&buf8);
   auto buf9 = CLContext::buffer<cl_int>(blSize*bufferSize);
   kern.set_arg(11,&buf9);
   auto buf10 = CLContext::buffer<cl_int>(blSize*bufferSize);
   kern.set_arg(12,&buf10);
   auto buf11 = CLContext::buffer<cl_int>(blSize*bufferSize);
   kern.set_arg(13,&buf11);
   kern.set_swarm_dims(1);
   kern.set_swarm_size(0,blSize*wSize,wSize);
   struct
   {
      State st;
      // x and y are the two chosen rows that will be compared to get the ans
      // Row '0'
      int x;
      // Row +1
      int y;
      /* Part of event looping yo
       * Part of the open CL listening to see when the kernel is done being
       * crunched
       * TODO: describe more
      */
      Ace::CLEvent ev;
      // an int (is this the offset in the list)?... TODO: what is this?
      AccelCompEng::CLBuffer<int> ld;
      // The OUT data
      AccelCompEng::CLBuffer<cl_float> ans;
   } state[smSize];

   for (int index = 0;index<smSize;++index)
   {
      state[index].st = State::start;
      state[index].x = 0;
      state[index].y = 0;
      state[index].ev = Ace::CLEvent();

      state[index].ld = CLContext::buffer<int>(2*blSize);
      state[index].ans = CLContext::buffer<cl_float>(blSize);
   }
   int alive {smSize};
   int si {0};
   auto i = _out->begin();
   unsigned long done {0};
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
               state[si].ld[2*count] = i.x()*size;
               state[si].ld[(2*count)+1] = i.y()*size;
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
            // An index to a portion of the explist ...
            // fetch list in the cl code actually populates the portion
            // of the 'two' lists
            kern.set_arg(3,&(state[si].ld));
            // The 'answer' the correlation score?
            //TODO: what's the formal term?
            kern.set_arg(5,&(state[si].ans));
            state[si].ev = CLCommandQueue::add_swarm(kern);
            state[si].st = State::exec;
         }
         break;
      case State::exec:
         // Is this polling? spin waiting for the event to be done?
         // Or is there an event firer?
         if (state[si].ev.is_done())
         {
            state[si].ev = CLCommandQueue::read_buffer(state[si].ans);
            state[si].st = State::out;
         }
         break;
         // Placeholder
      case State::out:
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
         break;
      }
      ++si;
      if (si>=smSize)
      {
         si = 0;
      }
      tm << "\rCalculating pearson values and saving to output file[" << 100*done/total
         << "%]..." << Ace::Terminal::flush;
      usleep(100);
   }
   tm << "\n";
}
