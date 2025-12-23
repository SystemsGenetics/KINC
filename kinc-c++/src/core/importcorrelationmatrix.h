#ifndef IMPORTCORRELATIONMATRIX_H
#define IMPORTCORRELATIONMATRIX_H
#include <ace/core/core.h>

#include "ccmatrix_pair.h"
#include "ccmatrix.h"
#include "correlationmatrix_pair.h"
#include "correlationmatrix.h"



/*!
 * This class implements the import correlation matrix analytic. This analytic
 * reads in a text file of correlations, where each line is a correlation that
 * includes the pairwise index, correlation value, sample mask, and several
 * other summary statistics. This analytic produces two data objects: a
 * correlation matrix containing the pairwise correlations, and a cluster matrix
 * containing the sample masks for each pairwise cluster. There are several
 * fields which are not represented in the text file and therefore must be
 * specified manually, including the gene size, sample size, max cluster size,
 * and correlation name.
 */
class ImportCorrelationMatrix : public EAbstractAnalytic
{
    Q_OBJECT
public:
    class Input;
    virtual int size() const override final;
    virtual void process(const EAbstractAnalyticBlock* result) override final;
    virtual EAbstractAnalyticInput* makeInput() override final;
    virtual void initialize();
private:
    /**
     * Workspace variables to read from the input file.
     */
    QTextStream _stream;
    int _numLines {0};
    Pairwise::Index _index {0};
    CCMatrix::Pair _ccmPair;
    CorrelationMatrix::Pair _cmxPair;
    /*!
     * Pointer to the input text file.
     */
    QFile* _input {nullptr};
    /*!
     * Pointer to the output cluster matrix.
     */
    CCMatrix* _ccm {nullptr};
    /*!
     * Pointer to the output correlation matrix.
     */
    CorrelationMatrix* _cmx {nullptr};
    /*!
     * The number of genes in the correlation matrix.
     */
    qint32 _geneSize {0};
    /*!
     * The maximum number of clusters allowed in a single pair of the
     * correlation matrix.
     */
    qint32 _maxClusterSize {1};
    /*!
     * The number of samples in the sample masks of the cluster matrix.
     */
    qint32 _sampleSize {0};
    /*!
     * The name of the correlation used in the correlation matrix.
     */
    QString _correlationName;
};



#endif
