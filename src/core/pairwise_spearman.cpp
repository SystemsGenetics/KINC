#include "pairwise_spearman.h"



using namespace Pairwise;






void Spearman::initialize(ExpressionMatrix* input)
{
   // pre-allocate workspace
   int workSize = nextPower2(input->sampleSize());

   _x.resize(workSize);
   _y.resize(workSize);
   _rank.resize(workSize);
}






float Spearman::computeCluster(
   const QVector<Vector2>& data,
   const QVector<qint8>& labels,
   qint8 cluster,
   int minSamples)
{
   // extract samples in gene pair cluster
   int N_pow2 = nextPower2(labels.size());
   int n = 0;

   for ( int i = 0, j = 0; i < labels.size(); ++i )
   {
      if ( labels[i] >= 0 )
      {
         if ( labels[i] == cluster )
         {
            _x[n] = data[j].s[0];
            _y[n] = data[j].s[1];
            _rank[n] = n + 1;
            ++n;
         }

         ++j;
      }
   }

   for ( int i = n; i < N_pow2; ++i )
   {
      _x[i] = INFINITY;
      _y[i] = INFINITY;
      _rank[i] = 0;
   }

   // compute correlation only if there are enough samples
   float result = NAN;

   if ( n >= minSamples )
   {
      // get new power of 2 floor size
      int n_pow2 = nextPower2(n);

      // execute two bitonic sorts that is beginning of spearman algorithm
      bitonicSort(n_pow2, _x, _y);
      bitonicSort(n_pow2, _y, _rank);

      // go through spearman sorted rank list and calculate difference from 1,2,3,... list
      int diff = 0;

      for ( int i = 0; i < n; ++i )
      {
         int tmp = (i + 1) - _rank[i];
         diff += tmp*tmp;
      }

      // compute spearman coefficient
      result = 1.0 - 6.0 * diff / (n * (n*n - 1));
   }

   return result;
}






int Spearman::nextPower2(int n)
{
   int pow2 = 2;
   while ( pow2 < n )
   {
      pow2 *= 2;
   }

   return pow2;
}






void Spearman::bitonicSort(int size, QVector<float>& sortList, QVector<float>& extraList)
{
   // initialize all variables
   int bsize = size/2;

   // bitonic algorithm, starting with an outer block of 2 and working up to total size of list
   for (int ob = 2; ob <= size; ob *= 2)
   {
      for (int ib = ob; ib >= 2; ib /= 2)
      {
         int t = ib/2;
         for (int i = 0; i < bsize; ++i)
         {
            int dir = -((i/(ob/2))&0x1);
            int a = (i/t)*ib+(i%t);
            int b = a+t;
            if ( ( ( sortList[a] > sortList[b] ) && !dir )
                 || ( ( sortList[a] < sortList[b] ) && dir ) )
            {
               std::swap(sortList[a], sortList[b]);
               std::swap(extraList[a], extraList[b]);
            }
         }
      }
   }
}






void Spearman::bitonicSort(int size, QVector<float>& sortList, QVector<int>& extraList)
{
   // initialize all variables
   int bsize = size/2;

   // bitonic algorithm, starting with an outer block of 2 and working up to total size of list
   for (int ob = 2; ob <= size; ob *= 2)
   {
      for (int ib = ob; ib >= 2; ib /= 2)
      {
         int t = ib/2;
         for (int i = 0; i < bsize; ++i)
         {
            int dir = -((i/(ob/2))&0x1);
            int a = (i/t)*ib+(i%t);
            int b = a+t;
            if ( ( ( sortList[a] > sortList[b] ) && !dir )
                 || ( ( sortList[a] < sortList[b] ) && dir ) )
            {
               std::swap(sortList[a], sortList[b]);
               std::swap(extraList[a], extraList[b]);
            }
         }
      }
   }
}
