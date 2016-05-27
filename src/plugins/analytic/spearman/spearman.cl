


void switch_f(__local float* a, __local float* b)
{
   float c = *a;
   *a = *b;
   *b = c;
}



void switch_i(__local int* a, __local int* b)
{
   int c = *a;
   *a = *b;
   *b = c;
}



void bitonic_sort_ff(__local float* sort, __local float* extra)
{
   int i = get_local_id(0);
   int bsize = get_local_size(0)<<1;
   int ob,ib,dir,a,b,t;
   for (ob=2;ob<=bsize;ob<<=1)
   {
      dir = -((i/(ob>>1))&0x1);
      for (ib=ob;ib>=2;ib>>=1)
      {
         t = ib>>1;
         a = (i/t)*ib+(i%t);
         b = a + t;
         if ((sort[a]>sort[b])^dir)
         {
            switch_f(&sort[a],&sort[b]);
            switch_f(&extra[a],&extra[b]);
         }
         barrier(CLK_LOCAL_MEM_FENCE);
      }
   }
}



void bitonic_sort_fi(__local float* sort, __local int* extra)
{
   int i = get_local_id(0);
   int bsize = get_local_size(0)<<1;
   int ob,ib,dir,a,b,t;
   for (ob=2;ob<=bsize;ob<<=1)
   {
      dir = -((i/(ob>>1))&0x1);
      for (ib=ob;ib>=2;ib>>=1)
      {
         t = ib>>1;
         a = (i/t)*ib+(i%t);
         b = a + t;
         if ((sort[a]>sort[b])^dir)
         {
            switch_f(&sort[a],&sort[b]);
            switch_i(&extra[a],&extra[b]);
         }
         barrier(CLK_LOCAL_MEM_FENCE);
      }
   }
}



float grab_correlation(__global float* correlations, int i, int size)
{
   if (i<size)
   {
      return correlations[i];
   }
   else
   {
      return INFINITY;
   }
}



long calc_rank_diff(long a, long b, int size)
{
   if (b<size)
   {
      return (a-b)*(a-b);
   }
   else
   {
      return 0;
   }
}



void accumulate(__local long* summation)
{
   int i = get_local_id(0);
   int hbsize = get_local_size(0);
   int b;
   for (b=hbsize;b>=1;b>>=1)
   {
      if (i<b)
      {
         summation[i] += summation[i+b];
      }
      barrier(CLK_LOCAL_MEM_FENCE);
   }
}



__kernel void spearman
(
      int aInd,
      int bInd,
      int size,
      __global float* correlations,
      __global float* result,
      __local float* alist,
      __local float* blist,
      __local int* rank,
      __local long* summation
)
{
   int ix = get_local_id(0)<<1;
   int nsize = get_local_size(0)<<1;
   alist[ix] = grab_correlation(correlations,ix+aInd,size);
   alist[ix+1] = grab_correlation(correlations,ix+aInd+1,size);
   blist[ix] = grab_correlation(correlations,ix+bInd,size);
   blist[ix+1] = grab_correlation(correlations,ix+bInd+1,size);
   rank[ix] = ix+1;
   rank[ix+1] = ix+2;
   barrier(CLK_LOCAL_MEM_FENCE);
   bitonic_sort_ff(alist,blist);
   bitonic_sort_fi(blist,rank);
   summation[ix] = calc_rank_diff(rank[ix],ix+1,size);
   summation[ix+1] = calc_rank_diff(rank[ix+1],ix+2,size);
   barrier(CLK_LOCAL_MEM_FENCE);
   accumulate(summation);
   if (get_local_id(0)==0)
   {
      *result = 1.0-(6.0*summation[0]/(float)(nsize*(nsize-1)));
   }
}
