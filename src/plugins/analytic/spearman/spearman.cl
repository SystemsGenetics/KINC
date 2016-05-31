


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



void build_lists(int aInd, int bInd, int size, __local float* alist,
                 __local float* blist, __global float* clist, __local int* rank)
{
   int ix = get_local_id(0)<<1;
   if (ix<size)
   {
      alist[ix] = clist[ix+aInd];
      blist[ix] = clist[ix+bInd];
      rank[ix] = ix+1;
   }
   else
   {
      alist[ix] = INFINITY;
      blist[ix] = INFINITY;
      rank[ix] = 0;
   }
   if ((ix+1)<size)
   {
      alist[ix+1] = clist[ix+aInd+1];
      blist[ix+1] = clist[ix+bInd+1];
      rank[ix+1] = ix+2;
   }
   else
   {
      alist[ix+1] = INFINITY;
      blist[ix+1] = INFINITY;
      rank[ix+1] = 0;
   }
   barrier(CLK_LOCAL_MEM_FENCE);
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
         if (((sort[a]>sort[b])&&!dir)||((sort[a]<sort[b])&&dir))
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
         if (((sort[a]>sort[b])&&!dir)||((sort[a]<sort[b])&&dir))
         {
            switch_f(&sort[a],&sort[b]);
            switch_i(&extra[a],&extra[b]);
         }
         barrier(CLK_LOCAL_MEM_FENCE);
      }
   }
}



long calc_ranks(int size, __local long* sums, __local int* ranks)
{
   int ix = get_local_id(0)<<1;
   long c;
   if (ix<size)
   {
      c = ranks[ix]-(ix+1);
      sums[ix] = c*c;
   }
   else
   {
      sums[ix] = 0.0;
   }
   if ((ix+1)<size)
   {
      c = ranks[ix+1]-(ix+2);
      sums[ix+1] = c*c;
   }
   else
   {
      sums[ix+1] = 0.0;
   }
   barrier(CLK_LOCAL_MEM_FENCE);
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
   build_lists(aInd,bInd,size,alist,blist,correlations,rank);
   bitonic_sort_ff(alist,blist);
   bitonic_sort_fi(blist,rank);
   calc_ranks(size,summation,rank);
   accumulate(summation);
   if (get_local_id(0)==0)
   {
      *result = 1.0-(6.0*summation[0]/((float)(size*(size*size-1))));
   }
}
