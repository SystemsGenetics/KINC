


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



void build_lists(int aInd, int bInd, int size, int chunk, __local float* aList,
                 __local float* bList, __global float* exprs)
{
   int i,c,ix;
   for (i=0;i<chunk;++i)
   {
      for (c=0;c<2;++c)
      {
         ix = (get_local_id(0)*chunk+i)*2+c;
         if (ix<size)
         {
            if (exprs[aInd+ix]==INFINITY||exprs[bInd+ix]==INFINITY)
            {
               aList[ix] = INFINITY;
               bList[ix] = INFINITY;
            }
            else
            {
               aList[ix] = exprs[aInd+ix];
               bList[ix] = exprs[bInd+ix];
            }
         }
         else
         {
            aList[ix] = INFINITY;
            bList[ix] = INFINITY;
         }
      }
   }
   barrier(CLK_LOCAL_MEM_FENCE);
}



void prune_list(int chunk, __local float* exprs, __local int* work, __local int* point)
{
   int i,c,ix;
   for (i=0;i<chunk;++i)
   {
      for (c=0;c<2;++c)
      {
         ix = (get_local_id(0)*chunk+i)*2+c;
         point[ix] = ix;
         if (exprs[ix]==INFINITY)
         {
            work[ix] = get_local_size(0)*2;
         }
         else
         {
            work[ix] = ix;
         }
      }
   }
}



void double_bitonic_sort_ii(int chunk, __local int* aSort, __local int* aExtra, __local int* bSort,
                            __local int* bExtra)
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
            if (((aSort[a]>aSort[b])&&!dir)||((aSort[a]<aSort[b])&&dir))
            {
               switch_i(&aSort[a],&aSort[b]);
               switch_i(&aExtra[a],&aExtra[b]);
            }
            if (((bSort[a]>bSort[b])&&!dir)||((bSort[a]<bSort[b])&&dir))
            {
               switch_i(&bSort[a],&bSort[b]);
               switch_i(&bExtra[a],&bExtra[b]);
            }
         }
         barrier(CLK_LOCAL_MEM_FENCE);
      }
   }
}



void construct_list(int chunk, __local float* exprs, __local float* list, __local int* point)
{
   int i,c,ix;
   for (i=0;i<chunk;++i)
   {
      for (c=0;c<2;++c)
      {
         ix = (get_local_id(0)*chunk+i)*2+c;
         list[ix] = exprs[point[ix]];
      }
   }
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



void make_debug(int chunk, __local float* alist, __local float* blist, __global float* atest,
                __global float* btest)
{
   int i,ix;
   for (i=0;i<chunk;++i)
   {
      ix = (get_local_id(0)*chunk+i)*2;
      atest[ix] = alist[ix];
      btest[ix] = blist[ix];
      atest[ix+1] = alist[ix+1];
      btest[ix+1] = blist[ix+1];
   }
   barrier(CLK_LOCAL_MEM_FENCE);
}



__kernel void spearman(int size, int chunk, __global int* insts, __global float* exprs,
                       __global float* result, __local float* alist, __local float* blist,
                       __local int* rank, __local long* summation, __global float* atest,
                       __global float* btest, __local float* aTmpList, __local float* bTmpList,
                       __local int* aWork, __local int* bWork, __local int* aPoint,
                       __local int* bPoint)
{
   int i = get_group_id(0)*2;
   build_lists(insts[i],insts[i+1],size,chunk,aTmpList,bTmpList,exprs);
   prune_list(chunk,aTmpList,aWork,aPoint);
   prune_list(chunk,bTmpList,bWork,bPoint);
   barrier(CLK_LOCAL_MEM_FENCE);
   double_bitonic_sort_ii(chunk,aWork,aPoint,bWork,bPoint);
   construct_list(chunk,aTmpList,alist,aPoint);
   construct_list(chunk,bTmpList,blist,bPoint);
   barrier(CLK_LOCAL_MEM_FENCE);
   make_debug(chunk,alist,blist,atest,btest);
   //bitonic_sort_ff(chunk,alist,blist);
   //bitonic_sort_fi(chunk,blist,rank);
   //calc_ranks(size,chunk,summation,rank);
   //accumulate(chunk,summation);
   //if (get_local_id(0)==0)
   //{
   //   *result = 1.0-(6.0*summation[0]/((float)(size*(size*size-1))));
   //}
}
