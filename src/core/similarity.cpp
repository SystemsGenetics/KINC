#include "similarity.h"
#include "similarity_input.h"
#include "similarity_resultblock.h"
#include "similarity_serial.h"
#include "similarity_workblock.h"
#include "similarity_opencl.h"
#include "similarity_cuda.h"
#include "ccmatrix_pair.h"
#include "correlationmatrix_pair.h"
#include <ace/core/ace_qmpi.h>
#include <ace/core/elog.h>



using namespace std;






/*!
 * Compute the next power of 2 which occurs after a number.
 *
 * @param n
 */
int Similarity::nextPower2(int n)
{
   int pow2 = 2;
   while ( pow2 < n )
   {
      pow2 *= 2;
   }

   return pow2;
}






/*!
 * Return the total number of pairs that must be processed for a given
 * expression matrix.
 *
 * @param emx
 */
qint64 Similarity::totalPairs(const ExpressionMatrix* emx)
{
   EDEBUG_FUNC(this,emx);

   return static_cast<qint64>(emx->geneSize()) * (emx->geneSize() - 1) / 2;
}






/*!
 * Return the total number of work blocks this analytic must process.
 */
int Similarity::size() const
{
   EDEBUG_FUNC(this);

   return (totalPairs(_input) + _workBlockSize - 1) / _workBlockSize;
}






/*!
 * Create and return a work block for this analytic with the given index. This
 * implementation creates a work block with a start index and size denoting the
 * number of pairs to process.
 *
 * @param index
 */
std::unique_ptr<EAbstractAnalyticBlock> Similarity::makeWork(int index) const
{
   EDEBUG_FUNC(this,index);

   if ( ELog::isActive() )
   {
      ELog() << tr("Making work index %1 of %2.\n").arg(index).arg(size());
   }

   qint64 start {index * static_cast<qint64>(_workBlockSize)};
   qint64 size {min(totalPairs(_input) - start, static_cast<qint64>(_workBlockSize))};

   return unique_ptr<EAbstractAnalyticBlock>(new WorkBlock(index, start, size));
}






/*!
 * Create an empty and uninitialized work block.
 */
std::unique_ptr<EAbstractAnalyticBlock> Similarity::makeWork() const
{
   EDEBUG_FUNC(this);

   return unique_ptr<EAbstractAnalyticBlock>(new WorkBlock);
}






/*!
 * Create an empty and uninitialized result block.
 */
std::unique_ptr<EAbstractAnalyticBlock> Similarity::makeResult() const
{
   EDEBUG_FUNC(this);

   return unique_ptr<EAbstractAnalyticBlock>(new ResultBlock);
}






/*!
 * Read in a block of results made from a block of work with the corresponding
 * index. This implementation takes the Pair objects in the result block and
 * saves them to the output correlation matrix and cluster matrix.
 *
 * This function converts the labels for each pair into the sample string format
 * that is used by the cluster matrix. In the pairwise label format, a label k
 * means that the sample belongs to cluster k, and a negative value means that
 * the sample was excluded for some reason. In the sample string format, there
 * is a separate sample string for each pairwise cluster, 0 and 1 denote whether
 * a sample belongs to that cluster, and all other values (6, 7, 8, 9) denote
 * that a sample was excluded for some other reason and correspond to the negative
 * values in the label format. For example, the following label array:
 *
 *   0 0 1 2 1 -9 2 -6
 *
 * is converted to the following sample string:
 *
 *   1 1 0 0 0 9 0 6 ,
 *   0 0 1 0 1 9 0 6 ,
 *   0 0 0 1 0 9 1 6
 *
 * @param result
 */
void Similarity::process(const EAbstractAnalyticBlock* result)
{
   EDEBUG_FUNC(this,result);

   if ( ELog::isActive() )
   {
      ELog() << tr("Processing result %1 of %2.\n").arg(result->index()).arg(size());
   }

   const ResultBlock* resultBlock {result->cast<ResultBlock>()};

   // iterate through all pairs in result block
   Pairwise::Index index {resultBlock->start()};

   for ( auto& pair : resultBlock->pairs() )
   {
      // save correlations that are within thresholds
      CCMatrix::Pair ccmPair(_ccm);
      CorrelationMatrix::Pair cmxPair(_cmx);

      for ( qint8 k = 0; k < pair.K; ++k )
      {
         // determine whether correlation is within thresholds
         float corr = pair.correlations[k];

         if ( !isnan(corr) && _minCorrelation <= abs(corr) && abs(corr) <= _maxCorrelation )
         {
            // save sample string
            ccmPair.addCluster();

            for ( int i = 0; i < _input->sampleSize(); ++i )
            {
               // convert label format to sample string format
               ccmPair.at(ccmPair.clusterSize() - 1, i) = (pair.labels[i] >= 0)
                  ? (k == pair.labels[i])
                  : -pair.labels[i];
            }

            // save correlation
            cmxPair.addCluster();
            cmxPair.at(cmxPair.clusterSize() - 1) = corr;
         }
      }

      if ( ccmPair.clusterSize() > 0 )
      {
         ccmPair.write(index);
      }

      if ( cmxPair.clusterSize() > 0 )
      {
         cmxPair.write(index);
      }

      ++index;
   }
}






/*!
 * Make a new input object and return its pointer.
 */
EAbstractAnalyticInput* Similarity::makeInput()
{
   EDEBUG_FUNC(this);

   return new Input(this);
}






/*!
 * Make a new serial object and return its pointer.
 */
EAbstractAnalyticSerial* Similarity::makeSerial()
{
   EDEBUG_FUNC(this);

   return new Serial(this);
}






/*!
 * Make a new OpenCL object and return its pointer.
 */
EAbstractAnalyticOpenCL* Similarity::makeOpenCL()
{
   EDEBUG_FUNC(this);

   return new OpenCL(this);
}






/*!
 * Make a new CUDA object and return its pointer.
 */
EAbstractAnalyticCUDA* Similarity::makeCUDA()
{
   EDEBUG_FUNC(this);

   return new CUDA(this);
}






/*!
 * Initialize this analytic. This implementation checks to make sure that valid
 * arguments were provided.
 */
void Similarity::initialize()
{
   EDEBUG_FUNC(this);

   // get MPI instance
   auto& mpi {Ace::QMPI::instance()};

   // only the master process needs to validate arguments
   if ( !mpi.isMaster() )
   {
      return;
   }

   // make sure input data is valid
   if ( !_input )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Did not get a valid input data object."));
      throw e;
   }

   // make sure cluster range is valid
   if ( _maxClusters < _minClusters )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Minimum clusters must be less than or equal to maximum clusters."));
      throw e;
   }

   // make sure kernel work sizes are valid
   if ( _globalWorkSize % _localWorkSize != 0 )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Global work size must be divisible by local work size."));
      throw e;
   }

   // initialize work block size
   if ( _workBlockSize == 0 )
   {
      int numWorkers = max(1, mpi.size() - 1);

      _workBlockSize = min(32768LL, totalPairs(_input) / numWorkers);
   }
}






/*!
 * Initialize the output data objects of this analytic.
 */
void Similarity::initializeOutputs()
{
   EDEBUG_FUNC(this);

   // make sure output data is valid
   if ( !_ccm || !_cmx )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Did not get valid output data objects."));
      throw e;
   }

   // initialize cluster matrix
   _ccm->initialize(_input->geneNames(), _maxClusters, _input->sampleNames());

   // initialize correlation matrix
   _cmx->initialize(_input->geneNames(), _maxClusters, _corrName);
}
