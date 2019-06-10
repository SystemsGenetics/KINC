#include "cluster_filter.h"
#include "cluster_filter_input.h"
#include "cluster_filter_serial.h"
#include "cluster_filter_resultblock.h"
#include "cluster_filter_workblock.h"
#include "datafactory.h"
#include "expressionmatrix_gene.h"
#include "ccmatrix_pair.h"
#include "correlationmatrix_pair.h"
#include <ace/core/ace_qmpi.h>
#include <ace/core/elog.h>



using namespace std;






/*!
 * Return the total number of blocks this analytic must process as steps
 * or blocks of work.
 */
int ClusterFilter::size() const
{
   EDEBUG_FUNC(this);

   return (_ccm->size() + _workBlockSize - 1) / _workBlockSize;

}


/*!
 * Create and return a work block for this analytic with the given index. This
 * implementation creates a work block with a start index and size denoting the
 * number of pairs to process.
 *
 * @param index
 */
std::unique_ptr<EAbstractAnalyticBlock> ClusterFilter::makeWork(int index) const
{
   EDEBUG_FUNC(this,index);

   if ( ELog::isActive() )
   {
      ELog() << tr("Making work index %1 of %2.\n").arg(index).arg(size());
   }

   qint64 start {index * static_cast<qint64>(_workBlockSize)};
   qint64 size {min(_ccm->size() - start, static_cast<qint64>(_workBlockSize))};

   return unique_ptr<EAbstractAnalyticBlock>(new WorkBlock(index, start, size));
}


/*!
 * Create an empty and uninitialized result block.
 */
std::unique_ptr<EAbstractAnalyticBlock> ClusterFilter::makeResult() const
{
   EDEBUG_FUNC(this);

   return unique_ptr<EAbstractAnalyticBlock>(new ResultBlock);
}



/*!
 * Create an empty and uninitialized work block.
 */
std::unique_ptr<EAbstractAnalyticBlock> ClusterFilter::makeWork() const
{
   EDEBUG_FUNC(this);

   return unique_ptr<EAbstractAnalyticBlock>(new WorkBlock);
}



/*!
 * Process the given index with a possible block of results if this analytic
 * produces work blocks. This implementation uses only the index of the result
 * block to determine which piece of work to do.
 *
 * @param result
 */
void ClusterFilter::process(const EAbstractAnalyticBlock*)
{
   EDEBUG_FUNC(this);


}

/*!
 * Make a new serial object and return its pointer.
 */
EAbstractAnalyticSerial* ClusterFilter::makeSerial()
{
   EDEBUG_FUNC(this);

   return new Serial(this);
}




/*!
 * Make a new input object and return its pointer.
 */
EAbstractAnalyticInput* ClusterFilter::makeInput()
{
   EDEBUG_FUNC(this);

   return new Input(this);
}






/*!
 * Initialize this analytic. This implementation checks to make sure the input
 * data objects and output file have been set.
 */
void ClusterFilter::initialize()
{
   EDEBUG_FUNC(this);

   // get MPI instance
   auto& mpi {Ace::QMPI::instance()};

   // only the master process needs to validate arguments
   if ( !mpi.isMaster() )
   {
      return;
   }

   // make sure input/output arguments are valid
   if ( !_emx || !_ccm || !_cmx )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Did not get valid input and/or output arguments."));
      throw e;
   }

   // initialize pairwise iterators

   // initialize output file stream
   //_stream.setDevice(_output);
   //_stream.setRealNumberPrecision(8);

   // initialize work block size
   if ( _workBlockSize == 0 )
   {
      int numWorkers = max(1, mpi.size() - 1);

      _workBlockSize = min(32768LL, _ccm->size() / numWorkers);
   }
}
