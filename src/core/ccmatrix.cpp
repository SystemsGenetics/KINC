#include "ccmatrix.h"
#include "ccmatrix_model.h"



/*!
 * Return a qt table model that represents this data object as a table.
 */
QAbstractTableModel* CCMatrix::model()
{
   if ( !_model )
   {
      _model = new Model(this);
   }
   return _model;
}






/*!
 * Initialize this cluster matrix with a list of gene names, the max cluster
 * size, and a list of sample names.
 *
 * @param geneNames
 * @param maxClusterSize
 * @param sampleNames
 */
void CCMatrix::initialize(const EMetadata& geneNames, int maxClusterSize, const EMetadata& sampleNames)
{
   // make sure sample names is an array and is not empty
   if ( !sampleNames.isArray() || sampleNames.toArray().isEmpty() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Sample names metadata is not an array or is empty."));
      throw e;
   }

   // save sample names to metadata
   EMetaObject metaObject {meta().toObject()};
   metaObject.insert("samples", sampleNames);
   setMeta(metaObject);

   // save sample size and initialize base class
   _sampleSize = sampleNames.toArray().size();
   Matrix::initialize(geneNames, maxClusterSize, (_sampleSize + 1) / 2 * sizeof(qint8), SUBHEADER_SIZE);
}






/*!
 * Return the list of correlation names in this correlation matrix.
 */
EMetadata CCMatrix::sampleNames() const
{
   return meta().toObject().at("samples");
}
