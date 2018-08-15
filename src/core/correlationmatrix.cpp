#include "correlationmatrix.h"
#include "correlationmatrix_model.h"
#include "correlationmatrix_pair.h"



QAbstractTableModel* CorrelationMatrix::model()
{
   if ( !_model )
   {
      _model = new Model(this);
   }
   return _model;
}






void CorrelationMatrix::initialize(const EMetadata &geneNames, int maxClusterSize, const EMetadata &correlationNames)
{
   // make sure correlation names is an array and is not empty
   if ( !correlationNames.isArray() || correlationNames.toArray().isEmpty() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Correlation names metadata is not an array or is empty."));
      throw e;
   }

   // save correlation names to metadata
   EMetaObject metaObject {meta().toObject()};
   metaObject.insert("correlations", correlationNames);
   setMeta(metaObject);

   // save correlation size and initialize base class
   _correlationSize = correlationNames.toArray().size();
   Matrix::initialize(geneNames, maxClusterSize, _correlationSize * sizeof(float), DATA_OFFSET);
}






EMetadata CorrelationMatrix::correlationNames() const
{
   return meta().toObject().at("correlations");
}






QVector<float> CorrelationMatrix::dumpRawData() const
{
   // if there are no genes do nothing
   if ( geneSize() == 0 )
   {
      return QVector<float>();
   }

   // create new correlation matrix
   QVector<float> data(geneSize() * geneSize() * maxClusterSize());

   // iterate through all pairs
   Pair pair(this);

   while ( pair.hasNext() )
   {
      // read in next pair
      pair.readNext();

      // load cluster data
      int i = pair.index().getX();
      int j = pair.index().getY();

      for ( int k = 0; k < pair.clusterSize(); ++k )
      {
         float correlation = pair.at(k, 0);

         data[i * geneSize() * maxClusterSize() + j * maxClusterSize() + k] = correlation;
         data[j * geneSize() * maxClusterSize() + i * maxClusterSize() + k] = correlation;
      }
   }

   return data;
}
