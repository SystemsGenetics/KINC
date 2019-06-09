#include "cluster_filter_serial.h"
#include "cluster_filter_resultblock.h"
#include "cluster_filter_workblock.h"
#include "expressionmatrix_gene.h"
#include "pairwise_gmm.h"
#include "pairwise_pearson.h"
#include "pairwise_spearman.h"
#include <ace/core/elog.h>



using namespace std;






/*!
 * Construct a new serial object with the given analytic as its parent.
 *
 * @param parent
 */
ClusterFilter::Serial::Serial(ClusterFilter* parent):
   EAbstractAnalyticSerial(parent),
   _base(parent)
{
   EDEBUG_FUNC(this,parent);

}






/*!
 * Read in the given work block and save the results in a new result block. This
 * implementation takes the starting pairwise index and pair size from the work
 * block and processes those pairs.
 *
 * @param block
 */
std::unique_ptr<EAbstractAnalyticBlock> ClusterFilter::Serial::execute(const EAbstractAnalyticBlock* block)
{
   EDEBUG_FUNC(this,block);


}




