#include "exportparametermatrix.h"
#include "exportparametermatrix_input.h"
#include "cpmatrix_pair.h"
#include "expressionmatrix_gene.h"
#include "pairwise_linalg.h"



using namespace Pairwise;






/*!
 * Return the total number of blocks this analytic must process as steps
 * or blocks of work. This implementation uses a work block for writing
 * each pair to the output file.
 */
int ExportParameterMatrix::size() const
{
   EDEBUG_FUNC(this);

   return _ccm->size();
}






/*!
 * Process the given index with a possible block of results if this analytic
 * produces work blocks. This implementation uses only the index of the result
 * block to determine which piece of work to do.
 *
 * @param result
 */
void ExportParameterMatrix::process(const EAbstractAnalyticBlock*)
{
   EDEBUG_FUNC(this);

   // read next pair
   _ccmPair.readNext();

   // read in gene expressions
   ExpressionMatrix::Gene gene1(_emx);
   ExpressionMatrix::Gene gene2(_emx);

   gene1.read(_ccmPair.index().getX());
   gene2.read(_ccmPair.index().getY());

   // initialize pairwise iterator for parameter matrix
   CPMatrix::Pair cpmPair;

   cpmPair.addCluster(_ccmPair.clusterSize());

   // compute parameters for each pairwise cluster
   for ( int k = 0; k < cpmPair.clusterSize(); k++ )
   {
      auto& component {cpmPair.at(k)};

      // compute mixture proportion
      int N = 0;
      int n_k = 0;

      for ( int i = 0; i < _ccm->sampleSize(); i++ )
      {
         N = (_ccmPair.at(k, i) == 0 || _ccmPair.at(k, i) == 1);
         n_k += (_ccmPair.at(k, i) == 1);
      }

      component.pi = n_k / N;

      // compute mean
      vectorInitZero(component.mu);

      for ( int i = 0; i < _ccm->sampleSize(); i++ )
      {
         // include only samples in cluster k
         if ( _ccmPair.at(k, i) == 1 )
         {
            Pairwise::Vector2 x_i = { gene1.at(i), gene2.at(i) };

            vectorAdd(component.mu, x_i);
         }
      }

      vectorScale(component.mu, 1.0f / n_k);

      // compute covariance matrix
      matrixInitZero(component.sigma);

      for ( int i = 0; i < _ccm->sampleSize(); i++ )
      {
         // include only samples in cluster k
         if ( _ccmPair.at(k, i) == 1 )
         {
            // compute xm = (x_i - mu_k)
            Vector2 xm = { gene1.at(i), gene2.at(i) };
            vectorSubtract(xm, component.mu);

            // compute Sigma_ki = (x_i - mu_k) (x_i - mu_k)^T
            matrixAddOuterProduct(component.sigma, 1, xm);
         }
      }

      matrixScale(component.sigma, 1.0f / n_k);
   }

   cpmPair.write(_ccmPair.index());
}






/*!
 * Make a new input object and return its pointer.
 */
EAbstractAnalyticInput* ExportParameterMatrix::makeInput()
{
   EDEBUG_FUNC(this);

   return new Input(this);
}






/*!
 * Initialize this analytic. This implementation checks to make sure the input
 * data objects and output file have been set.
 */
void ExportParameterMatrix::initialize()
{
   EDEBUG_FUNC(this);

   // make sure input/output arguments are valid
   if ( !_emx || !_ccm || !_cpm )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Did not get valid input and/or output arguments."));
      throw e;
   }

   // initialize pairwise iterators
   _ccmPair = CCMatrix::Pair(_ccm);
}
