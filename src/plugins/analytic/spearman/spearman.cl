


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



void build_lists(int aInd, int bInd, int size, int chunk, __local float* alist,
                 __local float* blist, __global float* clist, __local int* rank)
{
   int i,ix;
   for (i=0;i<chunk;++i)
   {
      ix = (get_local_id(0)*chunk+i)*2;
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
   }
   barrier(CLK_LOCAL_MEM_FENCE);
}



void bitonic_sort_ff(int chunk, __local float* sort, __local float* extra)
{
   int bsize = get_local_size(0)*2*chunk;
   int ob,ib,i,ix,dir,a,b,t;
   for (ob=2;ob<=bsize;ob*=2)
   {
      for (ib=ob;ib>=2;ib/=2)
      {
         for (i=0;i<chunk;++i)
         {
            ix = get_local_id(0)*chunk+i;
            dir = -((ix/(ob/2))&0x1);
            t = ib/2;
            a = (ix/t)*ib+(ix%t);
            b = a + t;
            if (((sort[a]>sort[b])&&!dir)||((sort[a]<sort[b])&&dir))
            {
               switch_f(&sort[a],&sort[b]);
               switch_f(&extra[a],&extra[b]);
            }
         }
         barrier(CLK_LOCAL_MEM_FENCE);
      }
   }
}



void bitonic_sort_fi(int chunk, __local float* sort, __local int* extra)
{
   int bsize = get_local_size(0)*2*chunk;
   int ob,ib,i,ix,dir,a,b,t;
   for (ob=2;ob<=bsize;ob*=2)
   {
      for (ib=ob;ib>=2;ib/=2)
      {
         for (i=0;i<chunk;++i)
         {
            ix = get_local_id(0)*chunk+i;
            dir = -((ix/(ob/2))&0x1);
            t = ib/2;
            a = (ix/t)*ib+(ix%t);
            b = a + t;
            if (((sort[a]>sort[b])&&!dir)||((sort[a]<sort[b])&&dir))
            {
               switch_f(&sort[a],&sort[b]);
               switch_i(&extra[a],&extra[b]);
            }
         }
         barrier(CLK_LOCAL_MEM_FENCE);
      }
   }
}



long calc_ranks(int size, int chunk, __local long* sums, __local int* ranks)
{
   long c;
   int i,ix;
   for (i=0;i<chunk;++i)
   {
      ix = (get_local_id(0)*chunk+i)*2;
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
   }
   barrier(CLK_LOCAL_MEM_FENCE);
}



void accumulate(int chunk, __local long* summation)
{
   int hbsize = get_local_size(0)*chunk;
   int i,ix,b;
   for (b=hbsize;b>=1;b/=2)
   {
      for (i=0;i<chunk;++i)
      {
         ix = get_local_id(0)*chunk+i;
         if (ix<b)
         {
            summation[ix] += summation[ix+b];
         }
      }
      barrier(CLK_LOCAL_MEM_FENCE);
   }
}



__kernel void spearman
(
      int aInd,
      int bInd,
      int size,
      int chunk,
      __global float* correlations,
      __global float* result,
      __local float* alist,
      __local float* blist,
      __local int* rank,
      __local long* summation
)
{
   build_lists(aInd,bInd,size,chunk,alist,blist,correlations,rank);
   bitonic_sort_ff(chunk,alist,blist);
   bitonic_sort_fi(chunk,blist,rank);
   calc_ranks(size,chunk,summation,rank);
   accumulate(chunk,summation);
   if (get_local_id(0)==0)
   {
      *result = 1.0-(6.0*summation[0]/((float)(size*(size*size-1))));
   }
}
