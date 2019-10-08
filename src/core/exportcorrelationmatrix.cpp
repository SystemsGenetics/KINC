#include "exportcorrelationmatrix.h"
#include "exportcorrelationmatrix_input.h"
#include "datafactory.h"
#include "expressionmatrix_gene.h"



using namespace std;






/*!
 * Return the total number of blocks this analytic must process as steps
 * or blocks of work. This implementation uses a work block for writing
 * each pair to the output file.
 */
int ExportCorrelationMatrix::size() const
{
   EDEBUG_FUNC(this);

   return _cmx->size();
}






/*!
 * Process the given index with a possible block of results if this analytic
 * produces work blocks. This implementation uses only the index of the result
 * block to determine which piece of work to do.
 *
 * @param result
 */
void ExportCorrelationMatrix::process(const EAbstractAnalyticBlock*)
{
   EDEBUG_FUNC(this);

   // initialize workspace
   QString sampleMask(_ccm->sampleSize(), '0');

   // read next pair
   _cmxPair.readNext();
   _ccmPair.read(_cmxPair.index());

   // write pairwise data to output file
   for ( int k = 0; k < _cmxPair.clusterSize(); k++ )
   {
      float correlation = _cmxPair.at(k);
      int numSamples = 0;
      int numMissing = 0;
      int numPostOutliers = 0;
      int numPreOutliers = 0;
      int numThreshold = 0;

      // if cluster data exists then use it
      if ( _ccmPair.clusterSize() > 0 )
      {
         // compute summary statistics
         for ( int i = 0; i < _ccm->sampleSize(); i++ )
         {
            switch ( _ccmPair.at(k, i) )
            {
            case 1:
               numSamples++;
               break;
            case 6:
               numThreshold++;
               break;
            case 7:
               numPreOutliers++;
               break;
            case 8:
               numPostOutliers++;
               break;
            case 9:
               numMissing++;
               break;
            }
         }

         // write sample mask to string
         for ( int i = 0; i < _ccm->sampleSize(); i++ )
         {
            sampleMask[i] = '0' + _ccmPair.at(k, i);
         }
      }

      // otherwise use expression data if provided
      else if ( _emx )
      {
         // read in gene expressions
         ExpressionMatrix::Gene gene1(_emx);
         ExpressionMatrix::Gene gene2(_emx);

         gene1.read(_cmxPair.index().getX());
         gene2.read(_cmxPair.index().getY());

         // determine sample mask, summary statistics from expression data
         for ( int i = 0; i < _emx->sampleSize(); ++i )
         {
            if ( isnan(gene1.at(i)) || isnan(gene2.at(i)) )
            {
               sampleMask[i] = '9';
               numMissing++;
            }
            else
            {
               sampleMask[i] = '1';
               numSamples++;
            }
         }
      }

      // otherwise throw an error
      else
      {
         E_MAKE_EXCEPTION(e);
         e.setTitle(tr("Invalid Input"));
         e.setDetails(tr("Expression Matrix was not provided but Cluster Matrix is missing sample data."));
         throw e;
      }

      // write cluster to output file
      _stream
         << _cmxPair.index().getX()
         << "\t" << _cmxPair.index().getY()
         << "\t" << k
         << "\t" << _cmxPair.clusterSize()
         << "\t" << numSamples
         << "\t" << numMissing
         << "\t" << numPostOutliers
         << "\t" << numPreOutliers
         << "\t" << numThreshold
         << "\t" << correlation
         << "\t" << sampleMask
         << "\n";
   }

   // make sure writing output file worked
   if ( _stream.status() != QTextStream::Ok )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Qt Text Stream encountered an unknown error."));
      throw e;
   }
}






/*!
 * Make a new input object and return its pointer.
 */
EAbstractAnalyticInput* ExportCorrelationMatrix::makeInput()
{
   EDEBUG_FUNC(this);

   return new Input(this);
}






/*!
 * Initialize this analytic. This implementation checks to make sure the input
 * data objects and output file have been set.
 */
void ExportCorrelationMatrix::initialize()
{
   EDEBUG_FUNC(this);

   // make sure input/output arguments are valid
   if ( !_ccm || !_cmx || !_output )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Did not get valid input and/or output arguments."));
      throw e;
   }

   // initialize pairwise iterators
   _ccmPair = CCMatrix::Pair(_ccm);
   _cmxPair = CorrelationMatrix::Pair(_cmx);

   // initialize output file stream
   _stream.setDevice(_output);
   _stream.setRealNumberPrecision(8);
}
