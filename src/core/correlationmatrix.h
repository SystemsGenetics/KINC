#ifndef CORRELATIONMATRIX_H
#define CORRELATIONMATRIX_H
#include "pairwise_matrix.h"



/*!
 * This class implements the correlation matrix data object. A correlation matrix
 * is a pairwise matrix where each pair-cluster element is a correlation value. The
 * matrix data can be accessed using the pairwise iterator for this class.
 */
class CorrelationMatrix : public Pairwise::Matrix
{
   Q_OBJECT
public:
   class Pair;
public:
   virtual QAbstractTableModel* model() override final;
public:
   void initialize(const EMetaArray& geneNames, int maxClusterSize, const EMetaArray& correlationNames);
   EMetaArray correlationNames() const;
   QVector<float> dumpRawData() const;
private:
   class Model;
private:
   /*!
    * Write the sub-header to the data object file.
    */
   virtual void writeHeader() { stream() << _correlationSize; }
   /*!
    * Read the sub-header from the data object file.
    */
   virtual void readHeader() { stream() >> _correlationSize; }
   /*!
    * The size (in bytes) of the sub-header. The sub-header consists of the
    * correlation size.
    */
   constexpr static int SUBHEADER_SIZE {1};
   /*!
    * The number of correlations in each pair-cluster.
    */
   qint8 _correlationSize {0};
   /*!
    * Pointer to a qt table model for this class.
    */
  Model* _model {nullptr};
};



#endif
