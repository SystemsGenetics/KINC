#ifndef CCMATRIX_H
#define CCMATRIX_H
#include "pairwise_matrix.h"



/*!
 * This class implements the cluster matrix data object. A cluster matrix is a
 * pairwise matrix where each pair-cluster element is a sample mask denoting
 * whether a sample belongs in the cluster. The matrix data can be accessed
 * using the pairwise iterator for this class.
 */
class CCMatrix : public Pairwise::Matrix
{
   Q_OBJECT
public:
   class Pair;
public:
   virtual QAbstractTableModel* model() override final;
public:
   void initialize(const EMetaArray& geneNames, int maxClusterSize, const EMetaArray& sampleNames);
   EMetaArray sampleNames() const;
   /*!
    * Return the number of samples in the cluster matrix.
    */
   int sampleSize() const { return _sampleSize; }
private:
   class Model;
private:
   /*!
    * Write the sub-header to the data object file.
    */
   virtual void writeHeader() { stream() << _sampleSize; }
   /*!
    * Read the sub-header from the data object file.
    */
   virtual void readHeader() { stream() >> _sampleSize; }
   /*!
    * The size (in bytes) of the sub-header. The sub-header consists of the
    * sample size.
    */
   constexpr static int SUBHEADER_SIZE {4};
   /*!
    * The number of samples in each sample mask.
    */
   qint32 _sampleSize {0};
   /*!
    * Pointer to a qt table model for this class.
    */
  Model* _model {nullptr};
};



#endif
