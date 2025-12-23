#ifndef CSM_H
#define CSM_H
#include "pairwise_matrix.h"



/*!
 * This class implements the condition-specific results matrix data object. A correlation matrix
 * is a pairwise matrix where each pair-cluster element is a combination of two values: p-value
 * and an r2 value from either from a statitical test denoting if the cluster is signicant for
 * any number of annotation categories/features. The matrix data can be accessed
 * using the pairwise iterator for this class.
 */
class CSMatrix : public Pairwise::Matrix
{
    Q_OBJECT
public:
    class Pair;
public:
    virtual QAbstractTableModel* model() override final;
public:
    void initialize(
        const EMetaArray& features,
        const QVector<EMetaArray>& featureInfo,
        const QVector<EMetaArray>& data,
        int numTests,
        QString testNames,
        const EMetaArray& geneNames,
        int maxClusterSize,
        int subheader);

    int sampleSize() const;
    QString getTestName(int index) const;

    /*!
     * Retrieves the type of test given it's index. Can
     * be 'Quantitative', 'Ordinal', 'Categorical', 'None'
     * or 'Unknown'.
     */
    QString getTestType(int index) const;

    qint32 getTestCount();

private:
    class Model;
private:
    /*!
     * Write the sub-header to the data object file.
     */
    virtual void writeHeader() override final;
    /*!
     * Read the sub-header from the data object file.
     */
    virtual void readHeader() override final;
    /*!
     * The number of statistical tests performed.
     */
    qint32 _testcount {0};
    /*!
     * The number of samples.
     */
    qint64 _sampleSize {0};
    /*!
     * The size (in bytes) of the sub-header. The sub-header consists of the
     * sample size.
     */
    constexpr static qint16 SUBHEADER_SIZE {12};
    /*!
     * Pointer to a qt table model for this class.
     */
    Model* _model {nullptr};
};



#endif
