#include "cpmatrix.h"
#include "cpmatrix_model.h"
#include "cpmatrix_pair.h"



/*!
 * Return a qt table model that represents this data object as a table.
 */
QAbstractTableModel* CPMatrix::model()
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
 */
void CPMatrix::initialize(const EMetaArray& geneNames, int maxClusterSize)
{
   EDEBUG_FUNC(this,&geneNames,maxClusterSize);

   // initialize base class
   Matrix::initialize(geneNames, maxClusterSize, sizeof(float), SUBHEADER_SIZE);
}






/*!
 * Return a list of correlation pairs in raw form.
 */
std::vector<CPMatrix::RawPair> CPMatrix::dumpRawData() const
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
      rawPair.components.resize(pair.clusterSize());

      for ( int k = 0; k < pair.clusterSize(); ++k )
      {
         rawPair.components[k] = pair.at(k);
      }

      pairs.push_back(rawPair);
   }

   return pairs;
}
