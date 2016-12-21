// http://www.stat.wmich.edu/s216/book/node122.html

// Precondition: File data will have to be read in
// done timestamp: 4:00 pm 11-22-1016



/********************************DATA PARTITIONING*****************************/
// Initially we'll need dummy data (2 - 80 double entries)
/*****************************END DATA PARTITIONING****************************/



/*****************************PROCESS PARTITIONING*****************************/
// The processing block listed below will be 'placed' in here
// this part bridges with the data partitioning labeled above
/*****************************END PROCESS PARTITIONING*************************/



/*********************************PROCESSING**********************************/

/******************************PSEUDO CODE************************************/

/* Human readable Steps to get the Coefficient */
/*
   1. Multiply the current row to subsequent rows beneath it
    a. Add the the products from multiplying two rows into one sum
       This will be variable A (Sum of multiplying 2 rows)
   2. Sum the the current row and the row beneath the current seperately
    a. This will be variable B (sum of entries in current row) and C (sum of
       entries in lower row)
    b. Multiply B and C to get variable D
    c. Divide D by E (the integral number of entries (this should be the same
       for every row))
    d. Call the result of this F
   3. Subtract: 1. - 2. (or A - F) This gives the covariance of the two
      particular samples. We will call this end result G.
    a. If there is a number there that is not 'too close' to 0, then the two
       genes/transcripts may have some sort of relationship in their respective
       synthesis/construction/expression
   4. Square and add each entry of the respective row and the subsequent row as
      well
    a. Respectively these variables will be called H and I
  5. Sum each row and square the sum of the respective row and the subsequent
     row as well
    a. Respectively these variables will be called J and K
  6. Divide 5. (J and K) by E (the integral number of entries of the specific
     row (this number should be the same for EVERY row))
    a. Call the results of J/E and K/E: L and M respectively
  7. Subtract H by L (H - L) and I by M (I - M)
    a. Call the results N and O respectively
  8. Multiply N and O and take the square root of the product
    a. Call the result P
  9. Divide F by P (F / P)
    a. Call this the result the Pearson Coefficient
    b. OH SNAP! We're done
*/

// Possible # of Intermediate Variables: 16 intermediate variables

/*****************************CL FUNCTION CALLS*******************************/

/* Toes into OpenCL */
/* Precondition: Have two different rows in scope of this kernel */
/*
  Update: This pseudo code is using vectors, but in the actual implementation
          the code is using Josh Burns' indexed based access of the data
*/

/*
  Functions that can used in acquiring the above variables
  Variable A: Use the dot product
    The largest vector using the 1.2 standard is the 16 wide vector
    16*5 = 80. Split the row into 5 sets of 'dotting' double16 vectors add the
    5 sepearate scalars together.
    ***(SEPEARATE FUNCTION)
  Variable B and C: Sum the rows (5 sets of double16 vectors) into scalars
    ***(SEPEARATE FUNCTION)
  Variable D: Scalar multiplication B * C store in D
    ***(IN PEARSON MAIN FUNCTION)
  Variable E: This should be the number 80 (the sSize!)
    ***(IN PEARSON MAIN FUNCTION)
  Variable F: Scalar division D / E store in F
    ***(IN PEARSON MAIN FUNCTION)
  Variable G: Scalar subtraction A - F store in G
    ***(IN PEARSON MAIN FUNCTION)
  Variable H and I: Square atomically (each of the elements in the vectors) the
    5 double16 vectors add their respective squares together.
    ***(SEPEARATE FUNCTION)
  Variable J and K: Square scalars B and C and store in variables (J & K)
    ***(IN PEARSON MAIN FUNCTION)
  Variable L and M: Scalar Division (L = J / E) and (M = K / E)
    ***(IN PEARSON MAIN FUNCTION)
  Variable N and O: Scalar subtraction (N = H - L) & (O = I - M)
    ***(IN PEARSON MAIN FUNCTION)
  Variable P: Scalar multiplication and root (sqrt(N*O))
    ***(IN PEARSON MAIN FUNCTION)

  Variable FIN: Scalar Division (F / P)

*/

/*******************************END PROCESSING********************************/


/*
  The data prep code simply extracts two of the rows from the 1-D list that was
  passed to a given kernel
  TODO: Ask Josh for clarification on how the offsetting works EXACTLY
*/
/*******************************BEGIN DATA PREP*******************************/
/*
  Plagerized code from Josh Burns' spearman cl code
  Fetch the expression scores from two rows (these are the two rows to compare)
*/

void fetch_lists(int aInd, int bInd, int size, int chunk, __global float* aList,
                 __global float* bList, __global float* exprs)
{
   int i,c,ix;
   for (i=0;i<chunk;++i)
   {
     // Why 2? for the max? for the 'row'?
     // row 0 and row 1?
     //
      for (c=0;c<2;++c)
      {
         ix = (get_local_id(0)*chunk+i)*2+c;
         if (ix<size)
         {
            if (isnan(exprs[aInd+ix]) || isnan(exprs[bInd+ix]))
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


/*
  Plagerized code from Joshi Burns' Spearman ACE plugin implementation
  Another piece of prepping the data (from NAN's)
*/
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
/*******************************END DATA PREP*********************************/


/*****************************BEGIN DOT ROWS**********************************/
/*
  To be used in the construction of variable A
  var dot_product will be passed in initialized to 0
*/
void dot_rows(int chunk,  __global float* aList,  __global float* bList,
             __global float* dot_product)
{
  int i,c,ix;
   // Plagerized from Josh Burns' codebase
   // These two for loops are the basic way of accessing the two 'rows'
   for (i=0;i<chunk;++i)
   {
      for (c=0;c<2;++c)
      {
         ix = (get_local_id(0)*chunk+i)*2+c;
         // Multiply element from a with element from b and add to scalar total
         dot_product += aList[ix] * bList[ix];
      }
   }
   barrier(CLK_LOCAL_MEM_FENCE);
}
/*******************************END DOT ROWS**********************************/


/*****************************BEGIN SUM ROWS**********************************/
/*
  To be used in the construction of variables B & C
  var sum will be passed in initialized to 0
*/
void sum_rows_individually(int chunk,  __global float* aList,  __global float* bList,
                           __global float* asum, __global float* bsum)
{
  int i,c,ix;
   // Plagerized from Josh Burns' codebase
   // These two for loops are the basic way of accessing the two 'rows'
   for (i=0;i<chunk;++i)
   {
      for (c=0;c<2;++c)
      {
         ix = (get_local_id(0)*chunk+i)*2+c;
         // Sum the two rows individually and store two seperate sums
         asum += aList[ix];
         bsum += bList[ix];
      }
   }
   barrier(CLK_LOCAL_MEM_FENCE);
}
/*******************************END SUM ROWS**********************************/


/*
  To be used in the construction of variables H & I
*/
void sqaure_and_sum_row_individually(int chunk,  __global float* aList,
                                     __global float* bList, __global float* asq_sum,
                                     __global float* bsq_sum)
{
  int i,c,ix;
   // Plagerized from Josh Burns' codebase
   // These two for loops are the basic way of accessing the two 'rows'
   for (i=0;i<chunk;++i)
   {
      for (c=0;c<2;++c)
      {
         ix = (get_local_id(0)*chunk+i)*2+c;
         // Sum the square of each element in the two rows individually and
         // store two seperate sums
         asq_sum += aList[ix] * aList[ix];
         bsq_sum += bList[ix] * bList[ix];
      }
   }
   barrier(CLK_LOCAL_MEM_FENCE);
}



/*
  Blindly relying that the partitioning will properly slice my
  data properly
  ***praise be to the Yoshi of Burné
*/

__kernel void spearman(/*kern_arg0*/int size,/*kern_arg1*/ int chunk,
                       /*kern_arg2*/int minSize, /*kern_arg3 TODO: what is this?*/__global int* insts,
                       /*kern_arg4 (INPUT 1-D list)*/__global float* exprs,
                       /*kern_arg5 (OUTPUT, relationship of the two rows)*/__global float* result,
                       /*kern_arg6*/ __global float* alistF, /*kern_arg7*/__global float* blistF,
                      //  __global int* rankF, __global int* iRankF, __global long* summationF,
                       /*kern_arg8*/__global float* aTmpListF, /*kern_arg9*/__global float* bTmpListF,
                       /*kern_arg10*/__global int* aWorkF, /*kern_arg11*/__global int* bWorkF,
                       /*kern_arg12*/__global int* aPointF, /*kern_arg13*/__global int* bPointF)
{
   int i = get_group_id(0)*2;
   int j = get_group_id(0);
   int wsize = get_local_size(0)*2*chunk;
   // '1st' row initialize size
   __global float* alist = &alistF[j*wsize];
   // '2nd' row initialize size
   __global float* blist = &blistF[j*wsize];
   // not needed, used in ranking of spearman
   //__global int* rank = &rankF[j*wsize];
   //__global int* iRank = &iRankF[j*wsize];
   //__global long* summation = &summationF[j*wsize];

   // These are to be the data prepped list of the two rows of the particular
   // kernel before bing pruned
   __global float* aTmpList = &aTmpListF[j*wsize];
   __global float* bTmpList = &bTmpListF[j*wsize];
   // The pruned  of each row (after dealing with the NANs?)
   // TODO: Verify this
   __global int* aWork = &aWorkF[j*wsize];
   __global int* bWork = &bWorkF[j*wsize];
   // TODO: What are these? (Looks like another buffer for scratch space)
   // I need them for the prune_lists ...
   __global int* aPoint = &aPointF[j*wsize];
   __global int* bPoint = &bPointF[j*wsize];

   // what is insts? the 'row?' it looks like it
   fetch_lists(insts[i],insts[i+1],size,chunk,aTmpList,bTmpList,exprs);
   // I don't think I need prune_lists, TODO: what does it do then? 
   prune_lists(chunk,aTmpList,aWork,aPoint,bTmpList,bWork,bPoint);




   // Not needed
   /*
   double_bitonic_sort_ii(chunk,aWork,aPoint,bWork,bPoint);
   construct_lists(chunk,aTmpList,alist,aPoint,bTmpList,blist,bPoint,rank,iRank);
   bitonic_sort_ff(chunk,alist,blist);
   bitonic_sort_fi(chunk,blist,rank);
   calc_ranks(chunk,summation,rank,iRank);
   accumulate(chunk,summation,iRank);
  */

   // At this point I assume that alist and blist hold my two rows to compared
   // extracted from fetch_lists and screened NANs from prune_lists
   // here I run my sub functions


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
