#include "ccmatrix.h"
#include "ccmatrix_model.h"



QAbstractTableModel* CCMatrix::model()
{
   if ( !_model )
   {
      _model = new Model(this);
   }
   return _model;
}






void CCMatrix::initialize(const EMetadata &geneNames, int maxClusterSize, const EMetadata &sampleNames)
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
   Matrix::initialize(geneNames, maxClusterSize, (_sampleSize + 1) / 2 * sizeof(qint8), DATA_OFFSET);
}






EMetadata CCMatrix::sampleNames() const
{
   return meta().toObject().at("samples");
}
