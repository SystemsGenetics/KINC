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
 * size, and a list of correlation names.
 *
 * @param geneNames
 * @param maxClusterSize
 * @param correlationNames
 */
void CorrelationMatrix::initialize(const EMetaArray& geneNames, int maxClusterSize, const EMetaArray& correlationNames)
{
   EDEBUG_FUNC(this,&geneNames,maxClusterSize,&correlationNames);

   // make sure correlation names is not empty
   if ( correlationNames.isEmpty() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Correlation names metadata is empty."));
      throw e;
   }

   // save correlation names to metadata
   EMetaObject metaObject {meta().toObject()};
   metaObject.insert("correlations", correlationNames);
   setMeta(metaObject);

   // save correlation size and initialize base class
   _correlationSize = correlationNames.size();
   Matrix::initialize(geneNames, maxClusterSize, _correlationSize * sizeof(float), SUBHEADER_SIZE);
}






/*!
 * Return the list of correlation names in this correlation matrix.
 */
EMetaArray CorrelationMatrix::correlationNames() const
{
   EDEBUG_FUNC(this);

   return meta().toObject().at("correlations").toArray();
}






/*!
 * Return a list of correlation pairs in raw form.
 */
QVector<CorrelationMatrix::RawPair> CorrelationMatrix::dumpRawData() const
{
   EDEBUG_FUNC(this);

   // create list of raw pairs
   QVector<RawPair> pairs;
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
         rawPair.correlations[k] = pair.at(k, 0);
      }
      
      pairs.append(rawPair);
   }

   return pairs;
}
