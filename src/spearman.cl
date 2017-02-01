


/*
ALGORITHM DESCRIPTION


STEP 1: POPULATE LISTS
Grab sub-arrays from global super array to populate two lists of samples, making infinity on both
list values if either value from the super array is NaN.

EXAMPLE: (A increment 0, B increment 1, size 8)
global list: 1.1 1.2 NaN 1.3 1.4 1.5 NaN 1.6 1.7
aList:       1.1 Inf Inf 1.3 1.4 Inf Inf 1.6
bList:       1.2 Inf Inf 1.4 1.5 Inf Inf 1.7


STEP 2: PRE-RANK LISTS FOR SORTING OUT INFINITIES
Making ranking and index numbers for A and B lists that will be used to sort out all infinities to
end of list while preserving order of real values.

EXAMPLE:
aList:       1.1 Inf Inf 1.3 1.4 Inf Inf 1.6
bList:       1.2 Inf Inf 1.4 1.5 Inf Inf 1.7
aWork:       1   32  32  4   5   32  32  8
bWork:       1   32  32  4   5   32  32  8
aPoint:      1   2   3   4   5   6   7   8
bPoint:      1   2   3   4   5   6   7   8


STEP 3: SORT WORK AND POINT LISTS TO TAKE OUT INFINITIES
Sort the point(index) and work(ranking) using the work values that will be used to take out all
infinities.

EXAMPLE:
aWork:       1   4   5   8   32  32  32  32
bWork:       1   4   5   8   32  32  32  32
aPoint:      1   4   5   8   2   3   6   7
bPoint:      1   4   5   8   2   3   6   7


STEP 4: MAKE NEW LISTS WITH INFINITY VALUES AT END
Now use point values to generate new sample lists where the infinities are sorted at the end and
order is preserved with real samples. The point value is the index to use for the new value. Also
construct two rank lists that will be used later for the spearman calculation. The rank lists
is simply the same value as the index, If the sample values are infinity
however the rank and irank values are zero.

EXAMPLE:
aList:       1.1 Inf Inf 1.3 1.4 Inf Inf 1.6
bList:       1.2 Inf Inf 1.4 1.5 Inf Inf 1.7
aPoint:      1   4   5   8   2   3   6   7
bPoint:      1   4   5   8   2   3   6   7
aTmpList:    1.1 1.3 1.4 1.6 Inf Inf Inf Inf
bTmpList:    1.2 1.4 1.5 1.7 Inf Inf Inf Inf
rank:        1   2   3   4   0   0   0   0
iRank:       1   2   3   4   0   0   0   0


STEP 5: FIRST BITONIC SORT
Sort both new sample lists by the values of aList.

EXAMPLE:
aTmpList:    1.1 1.3 1.4 1.6 Inf Inf Inf Inf
bTmpList:    1.2 1.4 1.5 1.7 Inf Inf Inf Inf


STEP 6: SECOND BITONIC SORT
Sort both the bTmpList sample list and the rank list using the values of the bTmpList.

EXAMPLE:
bTmpList:    1.2 1.4 1.5 1.7 Inf Inf Inf Inf
rank:        1   2   3   4   0   0   0   0


STEP 7: CALCULATE RANK DIFFERENCE
Now calculate the square of the difference between each rank and iRank value, putting the squared
difference in a new summation array. The iRank is also modified after the fact to have a value of
1 if it is not zero or else zero.

EXAMPLE:
rank:          1   2   3   4   0   0   0   0
iRank(before): 1   2   3   4   0   0   0   0
iRank(after):  1   1   1   1   0   0   0   0
sum:           0   0   0   0   0   0   0   0


STEP 8: SUMMATION
Now calculate the summation of sum and iRank, making the summation of all values be the value of the
first index in teh array.

EXAMPLE:
iRank:         4   1   1   1   0   0   0   0
sum:           0   0   0   0   0   0   0   0


STEP 9: SPEARMAN COEFFICIENT
Now using a single processing unit, calculate the spearman coefficient using the first index of
the summation array and the first index of iRank. The summation value is the summation of rank
differences. The first index of iRank is the total number of samples that have been correlated.
*/



// Switch the values of two variables.
void switch_f(__global float* a, __global float* b)
{
   float c = *a;
   *a = *b;
   *b = c;
}



// Switch the values of two variables.
void switch_i(__global int* a, __global int* b)
{
   int c = *a;
   *a = *b;
   *b = c;
}



// Fetch array of samples for both genes from global memory.
//
// @aInd Index into global mem list for alist
// @bInd Index into global mem list for blist
// @size Size of alist and blist.
// @chunk Number of samples each processing unit is responsible for.
// @aList Pointer to alist.
// @bList Pointer to blist.
// @exprs Pointer to global mem list.
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
              // If either entry is invalid (or both) invalidate that pair
               aList[ix] = INFINITY;
               bList[ix] = INFINITY;
            }
            else
            {
               // Populate aList and bList with the particular row data from the giant
               // 1-D list in global mem
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



// Takes two lists and populates four additional lists. Two lists each are populated using
// one of the input lists. The first list simply points to the index of the same index of the
// input list. The second is ordered the same unless the input list value is infinity; if it
// is infinity it assigned the second list value as a value four times greater than the highest
// value it could otherwise have.
//
// Example;
// Input list:  3.2 1.1 1.2 inf inf 1.1 inf
// First list:  0   1   2   3   4   5   6
// Second list: 0   1   2   28  28  5   28
//
// @chunk Number of items each processing unit is responsible for.
// @aExprs A input list.
// @aWork A second list.
// @aPoint A first list.
// @bExprs B input list.
// @bWork B second list.
// @bPoint B first list.
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
           // TODO: What does this do?
           // Josh: It gives any value with infinity the largest possible rank so it is sorted
           // out of the list.
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



// Bitonically sort two seperate lists at the same time. A second list for each list is also
// reorganized the same as the list being sorted.
//
// @chunk Number of items each processing unit is responsible for.
// @aSort First list to be sorted.
// @bExtra First extra list to be reorganized like the sorted list.
// @bSort Second list to be sorted.
// @bExtra Second extra list to be reorganized like the sorted list.
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



// Build new temporary sample lists that are used for spearman using the point lists. Also build
// rank and iRank lists.
//
// @chunk Number of items each processing unit is responsible for.
// @aExprs Sample list A to take samples from.
// @aList New sample list A to populate with sorted samples.
// @aPoint Point list A to use for building new sample list.
// @bExprs Sample list B to take samples from.
// @bList New sample list B to populate with sorted samples.
// @bPoint Point list B to use for building new sample list.
// @rank Rank list one.
// @iRank Rank list two.
void construct_lists(int chunk, __global float* aExprs, __global float* aList, __global int* aPoint,
                     __global float* bExprs, __global float* bList, __global int* bPoint,
                     __global int* rank, __global int* iRank)
{
   int i,c,ix;
   // The two for loops are the basic way of accessing the two 'rows'
   for (i=0;i<chunk;++i)
   {
      for (c=0;c<2;++c)
      {
         ix = (get_local_id(0)*chunk+i)*2+c;
         aList[ix] = aExprs[aPoint[ix]];
         bList[ix] = bExprs[bPoint[ix]];
         // check if the score is NAN
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



// Bitonically sort two lists at the same time using the values of the first list.
//
// @chunk Number of items each processing unit is responsible for.
// @sort First list to be sorted.
// @extra First extra list to be reorganized like the sorted list.
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



// Bitonically sort two lists at the same time using the values of the first list.
//
// @chunk Number of items each processing unit is responsible for.
// @sort First list to be sorted.
// @extra First extra list to be reorganized like the sorted list.
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



// Calulate differences between rank and iRank, putting the squared difference into a new
// summation array. Reset values of iRank to be one if not zero or else zero.
//
// @chunk Number of items each processing unit is responsible for.
// @sums Summation list where squared differences are stored.
// @rank Rank list.
// @iRank iRank list that will be reset to 1 if not zero else zero.
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



// Calculate summation of summation array and iRank array.
//
// @chunk Number of items each processing unit is responsible for.
// @sums Summation Summation list.
// @iRank iRank list.
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
   // Get data about group id and local size that is used a lot.
   int i = get_group_id(0)*2;
   int j = get_group_id(0);
   int wsize = get_local_size(0)*2*chunk;
   // Increment all pointers based off what workgroup this belongs to.
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
   // Populate sample lists of both genes, making any sample infinify if one or both gene samples
   // do not have a real sample.
   fetch_lists(insts[i],insts[i+1],size,chunk,aTmpList,bTmpList,exprs);
   // This function simple prepares the lists to be sorted, the Work lists are given a rank with
   // any sample that is infinity given the largest possible rank. The point lists are simple given
   // the increment into the list of samples pointing to what sample it represents.
   prune_lists(chunk,aTmpList,aWork,aPoint,bTmpList,bWork,bPoint);
   // This is where the magic happens. This sorts the work and point lists populated with the prune
   // function. As a result, all the infinity samples are sorted at the very end of the list while
   // preserving the order of all samples that are not infinity.
   double_bitonic_sort_ii(chunk,aWork,aPoint,bWork,bPoint);
   // This takes the now sorted work and point lists and uses them to populate the two arrays of
   // samples that will be used to calculate the spearman coefficient. This also builds the ranks
   // of each list that is used to get the spearman coefficient.
   construct_lists(chunk,aTmpList,alist,aPoint,bTmpList,blist,bPoint,rank,iRank);
   // Do the first sorting of the spearman algorithm, sorting by values of alist and changing blist
   // the same way.
   bitonic_sort_ff(chunk,alist,blist);
   // Do the second sorting of the spearman algorithm, sorting by valies of blist and changing rank
   // the same way.
   bitonic_sort_fi(chunk,blist,rank);
   // Now calculate the rank differences along with modifying the irank array.
   calc_ranks(chunk,summation,rank,iRank);
   // Get the summation of the rank differences and total number of iRanks that is one.
   accumulate(chunk,summation,iRank);
   // Lastly calculate the spearman coefficient.
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
