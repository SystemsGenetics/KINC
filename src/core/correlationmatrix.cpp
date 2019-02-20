#include "correlationmatrix.h"
#include "correlationmatrix_model.h"
#include "correlationmatrix_pair.h"



/*!
 * Return a qt table model that represents this data object as a table.
 */
QAbstractTableModel* CorrelationMatrix::model()
{
   EDEBUG_FUNC(this);

   if ( !_model )
   {
      _model = new Model(this);
   }
   return _model;
}






/*!
 * Initialize this correlation matrix with a list of gene names, the max cluster
 * size, and a correlation name.
 *
 * @param geneNames
 * @param maxClusterSize
 * @param correlationName
 */
void CorrelationMatrix::initialize(const EMetaArray& geneNames, int maxClusterSize, const QString& correlationName)
{
   EDEBUG_FUNC(this,&geneNames,maxClusterSize,&correlationName);

   // save correlation names to metadata
   EMetaObject metaObject {meta().toObject()};
   metaObject.insert("correlation", correlationName);
   setMeta(metaObject);

   // initialize base class
   Matrix::initialize(geneNames, maxClusterSize, sizeof(float), SUBHEADER_SIZE);
}






/*!
 * Return the correlation name for this correlation matrix.
 */
QString CorrelationMatrix::correlationName() const
{
   EDEBUG_FUNC(this);

   return meta().toObject().at("correlation").toString();
}






/*!
 * Return a list of correlation pairs in raw form.
 */
std::vector<CorrelationMatrix::RawPair> CorrelationMatrix::dumpRawData() const
{
   EDEBUG_FUNC(this);

   // create list of raw pairs
   std::vector<RawPair> pairs;
   pairs.reserve(size());

   // iterate through all pairs
   Pair pair(this);

   while ( pair.hasNext() )
   {
      // read in next pair
      pair.readNext();

      // copy pair to raw list
      RawPair rawPair;
      rawPair.index = pair.index();
      rawPair.correlations.resize(pair.clusterSize());

      for ( int k = 0; k < pair.clusterSize(); ++k )
      {
         rawPair.correlations[k] = pair.at(k);
      }
      
      pairs.push_back(rawPair);
   }

   return pairs;
}
