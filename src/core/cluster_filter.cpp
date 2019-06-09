#include "cluster_filter.h"
#include "cluster_filter_input.h"
#include "cluster_filter_serial.h"
#include "datafactory.h"
#include "expressionmatrix_gene.h"
#include <ace/core/ace_qmpi.h>
#include <ace/core/elog.h>



using namespace std;






/*!
 * Return the total number of blocks this analytic must process as steps
 * or blocks of work. This implementation uses a work block for writing
 * each pair to the output file.
 */
int ClusterFilter::size() const
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

   // make sure input/output arguments are valid
   if ( !_emx || !_ccm || !_cmx )
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
   //_stream.setDevice(_output);
   //_stream.setRealNumberPrecision(8);
}
