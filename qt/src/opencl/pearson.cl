





int fetchLists(__global float* expressions, int indexA, int indexB, __global float* listA
               , __global float* listB, int size)
{
   // initialize counters and indexes
   int i;
   int j = 0;
   int newSize = 2;
   indexA *= size;
   indexB *= size;

   // go through expression list with given indexes, generating new lists from it
   for (i = 0; i < size ;++i)
   {
      if ( !isnan(expressions[indexA+i]) && !isnan(expressions[indexB+i]) )
      {
         // if both expressions exist add expressions to new lists and increment
         listA[j] = expressions[indexA+i];
         listB[j] = expressions[indexB+i];
         j++;
      }
   }

   // set new size of generated lists and set unused end of lists to infinity
   newSize = j;
   for (i = j; i < size;++i)
   {
      listA[i] = INFINITY;
      listB[i] = INFINITY;
   }

   // return new size for generated lists
   return newSize;
}






__kernel void calculatePearsonBlock(int size, int minimumSize, __global float* expressions
                                    , __global float* workLists , __global int* targetList
                                    , __global float* resultList)
{
   // initialize variables used for pearson
   float aSum;
   float bSum;
   float a;
   float b;
   float ba;
   float bb;
   float ar;
   float br;

   // initialize counters and other variables
   int x;
   int newSize;
   int i = get_global_id(0);

   // initialize pointers to both work lists
   __global float* listA = &workLists[2*i*size];
   __global float* listB = &workLists[(2*i+1)*size];

   // populate work lists from expression data and make sure minimum size is reached
   newSize = fetchLists(expressions,targetList[2*i],targetList[2*i + 1],listA,listB,size);
   if ( newSize >= minimumSize )
   {
      // calculate summation of both lists
      aSum = 0.0;
      bSum = 0.0;
      for (x = 0; x < newSize ;++x)
      {
         aSum += listA[x];
         bSum += listB[x];
      }

      // divide summations by list size and initialize temporary variables
      aSum /= newSize;
      bSum /= newSize;
      a = 0;
      ba = 0;
      bb = 0;

      // iterate through both lists
      for (x = 0; x < newSize ;++x)
      {
         // conduct all math required for pearson
         ar = listA[x] - aSum;
         br = listB[x] - bSum;
         a += ar + br;
         ba += ar*ar;
         bb += br*br;
      }

      // calculate denominator of pearson and check to see if it is zero
      b = sqrt(ba*bb);
      if ( b == 0.0 )
      {
         // set result to zero
         resultList[i] = 0.0;
      }

      // else denominator is non-zero so calculate pearson with a/b
      else
      {
         resultList[i] = a/b;
      }
   }

   // else minimum is not reached and NAN is returned
   else
   {
      resultList[i] = NAN;
   }
}
