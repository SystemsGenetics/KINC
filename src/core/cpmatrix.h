#ifndef CPMATRIX_H
#define CPMATRIX_H
#include "pairwise_linalg.h"
#include "pairwise_matrix.h"



/*!
 * This class implements the correlation matrix data object. A correlation matrix
 * is a pairwise matrix where each pair-cluster element is a correlation value. The
 * matrix data can be accessed using the pairwise iterator for this class.
 */
class CPMatrix : public Pairwise::Matrix
{
   Q_OBJECT
public:
   class Pair;
public:
   struct Component
   {
      float pi;
      Pairwise::Vector2 mu;
      Pairwise::Matrix2x2 sigma;
   };
   struct RawPair
   {
      Pairwise::Index index;
      std::vector<Component> components;
   };
public:
   virtual QAbstractTableModel* model() override final;
public:
   void initialize(const EMetaArray& geneNames, int maxClusterSize);
   std::vector<RawPair> dumpRawData() const;
private:
   class Model;
private:
   /*!
    * Write the sub-header to the data object file.
    */
   virtual void writeHeader() override final {}
   /*!
    * Read the sub-header from the data object file.
    */
   virtual void readHeader() override final {}
   /*!
    * The size (in bytes) of the sub-header.
    */
   constexpr static int SUBHEADER_SIZE {0};
   /*!
    * Pointer to a qt table model for this class.
    */
   Model* _model {nullptr};
};



#endif
