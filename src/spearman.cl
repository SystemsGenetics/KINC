


void switch_f(__global float* a, __global float* b)
{
   float c = *a;
   *a = *b;
   *b = c;
}



void switch_i(__global int* a, __global int* b)
{
   int c = *a;
   *a = *b;
   *b = c;
}



void fetch_lists(int aInd, int bInd, int size, int chunk, __global float* aList,
                 __global float* bList, __global float* exprs)
{
   int i,c,ix;
   for (i=0;i<chunk;++i)
   {
      for (c=0;c<2;++c)
      {
         ix = (get_local_id(0)*chunk+i)*2+c;
         if (ix<size)
         {
            if (isnan(exprs[aInd+ix])||isnan(exprs[bInd+ix]))
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



void prune_lists(int chunk, __global float* aExprs, __global int* aWork, __global int* aPoint,
                 __global float* bExprs, __global int* bWork, __global int* bPoint)
{
   int i,c,ix;
   for (i=0;i<chunk;++i)
   {
      for (c=0;c<2;++c)
      {
         ix = (get_local_id(0)*chunk+i)*2+c;
         aPoint[ix] = ix;
         bPoint[ix] = ix;
         if (isinf(aExprs[ix]))
         {
            aWork[ix] = get_local_size(0)*4;
         }
         else
         {
            aWork[ix] = ix;
         }
         if (isinf(bExprs[ix]))
         {
            bWork[ix] = get_local_size(0)*4;
         }
         else
         {
            bWork[ix] = ix;
         }
      }
   }
   barrier(CLK_LOCAL_MEM_FENCE);
}



void double_bitonic_sort_ii(int chunk, __global int* aSort, __global int* aExtra, __global int* bSort,
                            __global int* bExtra)
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



void construct_lists(int chunk, __global float* aExprs, __global float* aList, __global int* aPoint,
                     __global float* bExprs, __global float* bList, __global int* bPoint,
                     __global int* rank, __global int* iRank)
{
   int i,c,ix;
   for (i=0;i<chunk;++i)
   {
      for (c=0;c<2;++c)
      {
         ix = (get_local_id(0)*chunk+i)*2+c;
         aList[ix] = aExprs[aPoint[ix]];
         bList[ix] = bExprs[bPoint[ix]];
         if (isinf(aList[ix]))
         {
            rank[ix] = 0;
            iRank[ix] = 0;
         }
         else
         {
            rank[ix] = ix+1;
            iRank[ix] = ix+1;
         }
      }
   }
   barrier(CLK_LOCAL_MEM_FENCE);
}



void bitonic_sort_ff(int chunk, __global float* sort, __global float* extra)
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



void bitonic_sort_fi(int chunk, __global float* sort, __global int* extra)
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



long calc_ranks(int chunk, __global long* sums, __global int* rank, __global int* iRank)
{
   long tmp;
   int i,c,ix;
   for (i=0;i<chunk;++i)
   {
      for (c=0;c<2;++c)
      {
         ix = (get_local_id(0)*chunk+i)*2+c;
         tmp = rank[ix] - iRank[ix];
         sums[ix] = tmp*tmp;
         if (rank[ix]==0)
         {
            iRank[ix] = 0;
         }
         else
         {
            iRank[ix] = 1;
         }
      }
   }
   barrier(CLK_LOCAL_MEM_FENCE);
}



void accumulate(int chunk, __global long* summation, __global int* iRank)
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
            iRank[ix] += iRank[ix+b];
         }
      }
      barrier(CLK_LOCAL_MEM_FENCE);
   }
}



__kernel void spearman(int size, int chunk, int minSize, __global int* insts, __global float* exprs,
                       __global float* result, __global float* alistF, __global float* blistF,
                       __global int* rankF, __global int* iRankF, __global long* summationF,
                       __global float* aTmpListF, __global float* bTmpListF, __global int* aWorkF,
                       __global int* bWorkF, __global int* aPointF, __global int* bPointF)
{
   int i = get_group_id(0)*2;
   int j = get_group_id(0);
   int wsize = get_local_size(0)*2*chunk;
   __global float* alist = &alistF[j*wsize];
   __global float* blist = &blistF[j*wsize];
   __global int* rank = &rankF[j*wsize];
   __global int* iRank = &iRankF[j*wsize];
   __global long* summation = &summationF[j*wsize];
   __global float* aTmpList = &aTmpListF[j*wsize];
   __global float* bTmpList = &bTmpListF[j*wsize];
   __global int* aWork = &aWorkF[j*wsize];
   __global int* bWork = &bWorkF[j*wsize];
   __global int* aPoint = &aPointF[j*wsize];
   __global int* bPoint = &bPointF[j*wsize];
   fetch_lists(insts[i],insts[i+1],size,chunk,aTmpList,bTmpList,exprs);
   prune_lists(chunk,aTmpList,aWork,aPoint,bTmpList,bWork,bPoint);
   double_bitonic_sort_ii(chunk,aWork,aPoint,bWork,bPoint);
   construct_lists(chunk,aTmpList,alist,aPoint,bTmpList,blist,bPoint,rank,iRank);
   bitonic_sort_ff(chunk,alist,blist);
   bitonic_sort_fi(chunk,blist,rank);
   calc_ranks(chunk,summation,rank,iRank);
   accumulate(chunk,summation,iRank);
   if (get_local_id(0)==0)
   {
      size = iRank[0];
      if (size<minSize)
      {
         result[j] = NAN;
      }
      else
      {
         result[j] = 1.0-(6.0*(float)summation[0]/((float)size*(((float)size*(float)size)-1)));
      }
   }
}
