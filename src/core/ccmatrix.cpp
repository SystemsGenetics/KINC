#include "ccmatrix.h"
#include "ccmatrix_pair.h"



using namespace std;






QAbstractTableModel* CCMatrix::model()
{
   return nullptr;
}






QVariant CCMatrix::headerData(int section, Qt::Orientation orientation, int role) const
{
   // orientation is not used
   Q_UNUSED(orientation);

   // if role is not display return nothing
   if ( role != Qt::DisplayRole )
   {
      return QVariant();
   }

   // get genes metadata and make sure it is an array
   const EMetadata& genes {geneNames()};
   if ( genes.isArray() )
   {
      // make sure section is within limits of gene name array
      if ( section >= 0 && section < genes.toArray().size() )
      {
         // return gene name
         return genes.toArray().at(section).toString();
      }
   }

   // no gene found return nothing
   return QVariant();
}






int CCMatrix::rowCount(const QModelIndex&) const
{
   return geneSize();
}






int CCMatrix::columnCount(const QModelIndex&) const
{
   return geneSize();
}






QVariant CCMatrix::data(const QModelIndex &index, int role) const
{
   // if role is not display return nothing
   if ( role != Qt::DisplayRole )
   {
      return QVariant();
   }

   // if row and column are equal return empty string
   if ( index.row() == index.column() )
   {
      return "";
   }

   // get constant pair and read in values
   const Pair pair(this);
   int x {index.row()};
   int y {index.column()};
   if ( y > x )
   {
      swap(x,y);
   }
   pair.read({x,y});

   // Return value of pair as a string
   return pair.toString();
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
