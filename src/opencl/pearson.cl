





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





/*
 *
 *
 *
 */
__kernel void calculatePearsonBlock(int size, int minimumSize, __global float* expressions
                                    , __global float* workLists , __global int* targetList
                                    , __global float* resultList)
{
   // Initialize variables used for pearson
   float aSum;
   float bSum;
   float aaSum;
   float bbSum;
   float abSum;
   float a;
   float b;


   // initialize counters and other variables
   int x;
   int newSize;
   int i = get_global_id(0);

   // initialize pointers to both work lists
   __global float* listA = &workLists[2*i*size];
   __global float* listB = &workLists[(2*i+1)*size];

   // Populate work lists from expression data and make sure minimum size is reached
   newSize = fetchLists(expressions,targetList[2*i],targetList[2*i + 1],listA,listB,size);
   if ( newSize >= minimumSize )
   {
      // Calculate summation of both lists
      aSum = aaSum = 0.0;
      bSum = bbSum = 0.0;
      abSum = 0.0;
      for (x = 0; x < newSize ;++x)
      {
         aSum += listA[x];
         bSum += listB[x];

         aaSum += (listA[x]*listA[x]);
         bbSum += (listB[x]*listB[x]);

         abSum += (listA[x]*listB[x]);
      }

      // Calculate numerator of pearson
      a = ((newSize)*abSum) - (aSum*bSum);
      // Calculate denominator of pearson and check to see if it is zero
      b = sqrt(((newSize*aaSum)-(aSum*aSum))*((newSize*bbSum)-(bSum*bSum)));
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
