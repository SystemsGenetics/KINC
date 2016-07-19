#include "spearman.h"
#include "spearman.cl.h"
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <ctime>



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



void Spearman::execute_cl(GetOpts&, Terminal& tm)
{
   AccelCompEng::assert<NoDataInput>(_in,__LINE__);
   AccelCompEng::assert<NoDataOutput>(_out,__LINE__);
   int gSize {_in->gSize()};
   int sSize {_in->sSize()};
   tm << "Loading kernel program into OpenCL device...\n";
   CLProgram::add_source(spearman_cl);
   if (!CLProgram::compile(""))
   {
      tm << CLProgram::log();
      return;
   }
   auto kern = CLProgram::mkernel("spearman");
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
   int bSize {pow2_wall(sSize)};
   bool cond {bSize>0};
   AccelCompEng::assert<TooManySamples>(cond,__LINE__);
   int wSize {1024};
   int chunk {1};
   if ((2*wSize)<bSize)
   {
      chunk = bSize/(2*wSize);
   }
   tm << chunk << "\n";
}



void Spearman::execute_pn(GetOpts&, Terminal& tm)
{
   tm << "Not yet implemented.\n";
}



int Spearman::pow2_wall(int i)
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

/*

void spearman::execute_cl(GetOpts& ops, Terminal& tm)
{
   assert<NotEnoughInputs>(_in.size()>0,__LINE__);
   assert<NoOutput>(_out,__LINE__);
   tm << "Writing header information..." << Terminal::flush;
   int geneSize {_in[0]->gene_size()};
   int sampleSize {_in[0]->sample_size()};
   if (_in.size()>1)
   {
      for (int i=1;i<_in.size();++i)
      {
         if (sampleSize!=_in[i]->sample_size())
         {
            throw DifferentSampleSizes(__LINE__);
         }
         geneSize += _in[i]->gene_size();
      }
   }
   _out->set_gene_size(geneSize);
   _out->set_sample_size(sampleSize);
   _out->set_correlation_size(1);
   _out->set_correlation_name(0,"spearman");
   for (int y=0;y<_in[0]->sample_size();++y)
   {
      tm << y << Terminal::endl;
      _out->set_sample_name(y,"");
   }
   auto glist = CLContext::buffer<cl_float>(geneSize*sampleSize);
   int x {0};
   int is {0};
   for (auto i:_in)
   {
      i->load_buffer();
      for (int y=0;y<i->gene_size();++y)
      {
         _out->set_gene_name(x++,i->gene_name(y));
         const float* corrs = i->gene(y);
         for (int z=0;z<sampleSize;++z)
         {
            glist[is++] = corrs[z];
         }
      }
   }
   _out->create_data();
   tm << "Done.\n";
   bool alive {true};
   int g1 {1};
   int g2 {0};
   long c {0};
   long t {geneSize*(geneSize-1)/2};
   std::vector<uint8_t> masks(sampleSize,1);
   std::vector<float> corr(1,0);
   CLProgram::add_source(spearman_cl);
   if (!CLProgram::compile(""))
   {
      tm << CLProgram::log();
      return;
   }
   auto ans = CLContext::buffer<cl_float>(1);
   auto kern = CLProgram::mkernel("spearman");
   kern.set_arg(2,sampleSize);
   kern.set_arg(3,1);
   kern.set_arg(4,&glist);
   kern.set_arg(5,&ans);
   kern.set_arg(6,sizeof(cl_float)*1024*2);
   kern.set_arg(7,sizeof(cl_float)*1024*2);
   kern.set_arg(8,sizeof(cl_int)*1024*2);
   kern.set_arg(9,sizeof(cl_long)*1024*2);
   kern.set_swarm_dims(1);
   kern.set_swarm_size(0,1024,1024);
   struct cac
   {
      int g1,g2;
      CLEvent ev;
      CLBuffer<cl_float> buf;
      int which {0};
      // 0 = launch new correlation
      // 1 = read result
      // 2 = write to file
   };
   cac channels[100];
   for (int i=0;i<100;++i)
   {
      channels[i].buf = CLContext::buffer<cl_float>(1);
   }
   CLEvent ev = CLCommandQueue::write_buffer(glist);
   ev.wait();
   kern.set_arg(0,g1*sampleSize);
   kern.set_arg(1,g2*sampleSize);
   ev = CLCommandQueue::add_swarm(kern);
   int chi {0};
   int HAHA {0};
   time_t START = time(NULL);
   while (alive||HAHA>0)
   {
      float perc = (float)c++/(float)t;
      float time_left = (float)(time(NULL)-START)/perc*(1-perc);
      tm << "\r                                                               ";
      tm << "\r[ML=" << (int)time_left/60 << "]Computing... ["
         << perc*100.0
         << "%](" << HAHA << ")" << Terminal::flush;
      switch(channels[chi].which)
      {
      case 0:
         if (alive)
         {
            channels[chi].g1 = g1;
            channels[chi].g2 = g2;
            kern.set_arg(5,&(channels[chi].buf));
            kern.set_arg(0,g1*sampleSize);
            kern.set_arg(1,g2++*sampleSize);
            channels[chi].ev = CLCommandQueue::add_swarm(kern);
            channels[chi].which = 1;
            HAHA++;
         }
         break;
      case 1:
         if (ev.is_done())
         {
            channels[chi].ev = CLCommandQueue::read_buffer(channels[chi].buf);
            channels[chi].which = 2;
         }
         break;
      case 2:
         if (ev.is_done())
         {
            corr[0] = channels[chi].buf[0];
            auto ptr = _out->set_modes(channels[chi].g2,channels[chi].g1,1);
            _out->write_mode(ptr,0,masks,corr);
            channels[chi].which = 0;
            HAHA--;
         }
         break;
      }
      if (g2>=g1)
      {
         g2 = 0;
         g1++;
         if (g1>=geneSize)
         {
            alive = false;
         }
      }
      chi++;
      if (chi>=100) chi = 0;
   }
   tm << "\r                                                               ";
   tm << "\rComputing... Done.\n";
}



*/
