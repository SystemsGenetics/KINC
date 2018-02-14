





// Fetch and build array of expressions for both genes, skipping any expressions that are missing
// for either gene.
//
// @param expressions Array of all expressions for all genes to generate gene lists from.
// @param size Size of lists for both genes.
// @param vector Index into expression list for gene A.
// @param sampleMask Array of values denoting membership in a cluster.
// @param x New array of expressions for gene A that this function builds.
// @param y New array of expressions for gene B that this function builds.
// @return Returns size of newly generated arrays which excludes any missing expression values.
int fetchData(
   __global const float *expressions,
   int size,
   int2 vector,
   __global const char *sampleMask,
   __global float *x,
   __global float *y)
{
   // initialize counters and indexes
   __global const float *gene1 = &expressions[vector.x * size];
   __global const float *gene2 = &expressions[vector.y * size];

   // populate a and b with shared expressions of gene pair
   int numSamples = 0;

   if ( sampleMask[0] != -1 )
   {
      // add samples that are in the cluster
      for ( int i = 0; i < size; ++i )
      {
         if ( sampleMask[i] == 1 )
         {
            x[numSamples] = gene1[i];
            y[numSamples] = gene2[i];
            ++numSamples;
         }
      }
   }
   else
   {
      // add samples that are valid
      for ( int i = 0; i < size; ++i )
      {
         if ( !isnan(gene1[i]) && !isnan(gene2[i]) )
         {
            x[numSamples] = gene1[i];
            y[numSamples] = gene2[i];
            ++numSamples;
         }
      }
   }

   // return new size for generated lists
   return numSamples;
}






// Calculate a block of pearson coefficients given a block of gene pairs.
//
// @param expressions Gene expression matrix
// @param size The size of the expressions/samples per gene.
// @param pairs Array of gene pairs.
// @param sampleMasks Array of sample masks for each gene pair.
// @param work Work space to be used for pearson calculations.
// @param results Array of output spearman coefficients for each gene pair.
__kernel void calculatePearsonBlock(
   __global const float *expressions,
   int size,
   __global const int2 *pairs,
   __global const char *sampleMasks,
   int minSamples,
   __global float *work,
   __global float *results)
{
   int i = get_global_id(0);

   if ( pairs[i].x == 0 && pairs[i].y == 0 )
   {
      return;
   }

   // initialize workspace variables
   __global const char *sampleMask = &sampleMasks[i * size];
   __global float *x = &work[(2*i+0) * size];
   __global float *y = &work[(2*i+1) * size];

   // fetch x and y arrays from expression matrix
   int n = fetchData(
      expressions, size,
      pairs[i],
      sampleMask,
      x, y
   );

   // compute correlation only if there are enough samples
   results[i] = NAN;

   if ( n >= minSamples )
   {
      // compute intermediate sums
      float sumx = 0;
      float sumy = 0;
      float sumx2 = 0;
      float sumy2 = 0;
      float sumxy = 0;

      for ( int i = 0; i < n; ++i )
      {
         sumx += x[i];
         sumy += y[i];
         sumx2 += x[i] * x[i];
         sumy2 += y[i] * y[i];
         sumxy += x[i] * y[i];
      }

      // compute Pearson correlation coefficient
      results[i] = (n*sumxy - sumx*sumy) / sqrt((n*sumx2 - sumx*sumx) * (n*sumy2 - sumy*sumy));
   }
}
