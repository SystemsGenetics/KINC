#ifndef NET_WRITER_H
#define NET_WRITER_H

#include <QObject>
#include <QTextStream>
#include "ccmatrix_pair.h"
#include "ccmatrix.h"
#include "correlationmatrix_pair.h"
#include "correlationmatrix.h"
#include "expressionmatrix.h"
#include "conditionspecificclustersmatrix.h"
#include "conditionspecificclustersmatrix_pair.h"

/*!
 * The parent class for writing a network to an output file.
 *
 * Child classes need only implement the initialize(), writeEdge() and
 * finish() functions.
 */
class NetWriter : public QObject
{
    Q_OBJECT
protected:
    QTextStream * _stream {nullptr};


    /*!
     * Pointer to the input expression matrix.
     */
    ExpressionMatrix* _emx {nullptr};
    /*!
     * Pointer to the input cluster matrix and corresponding pair object.
     */
    CCMatrix* _ccm {nullptr};
    CCMatrix::Pair _ccmPair;
    /*!
     * Pointer to the input correlation matrix and corresponding pair object.
     */
    CorrelationMatrix* _cmx {nullptr};
    CorrelationMatrix::Pair _cmxPair;
    /*!
     * Pointer to the input condition specific cluster matrix and
     * corresponding pair object.
     */
    CSMatrix* _csm {nullptr};
    CSMatrix::Pair _csmPair;

    /*!
     * Ther index of the current pair.
     */
    Pairwise::Index _index;

    /**
     * The variable that houses the test names for easy lookup and
     * the function that should be used to set the test names.
     */
    void setTestNames();
    QVector<QString> _testNames;

public:

    void setPair(Pairwise::Index cmx_index);
    QString getEdgeGene1();
    QString getEdgeGene2();
    QString getEdgeSampleString(int cluster_index);
    int getEdgeNumSamples(QString sample_string);
    float getEdgeSimilarity(int cluster_index);
    float getEdgeTestValue(int cluster_index, int test_index);
    QVector<QString> getTestNames() { return _testNames; }

    NetWriter(QTextStream * stream, ExpressionMatrix * emx,
                      CorrelationMatrix * cmx, CCMatrix * ccm,
                      CSMatrix * csm);
    virtual ~NetWriter() {}
    virtual void initialize() = 0;
    virtual void writeEdgeCluster(int cluster_index, QVector<QString> passed) = 0;
    virtual void finish() = 0;
};

class GMLNetWriter: public NetWriter
{
private:
    /*!
     * A hash lookup for storing nodes in the network.
     */
    QHash<QString, bool> _nodes;
public:
    GMLNetWriter(QTextStream * stream, ExpressionMatrix * emx, CorrelationMatrix * cmx,
                 CCMatrix * ccm, CSMatrix * csm) : NetWriter(stream, emx, cmx, ccm, csm) {}
    void initialize();
    void writeEdgeCluster(int cluster_index, QVector<QString> passed);
    void finish();
};

class FullNetWriter: public NetWriter
{
public:
    FullNetWriter(QTextStream * stream, ExpressionMatrix * emx, CorrelationMatrix * cmx,
                  CCMatrix * ccm, CSMatrix * csm) : NetWriter(stream, emx, cmx, ccm, csm) {}
    void initialize();
    void writeEdgeCluster(int cluster_index, QVector<QString> passed);
    void finish() {}
};

class TidyNetWriter: public NetWriter
{
public:
    TidyNetWriter(QTextStream * stream, ExpressionMatrix * emx, CorrelationMatrix * cmx,
                  CCMatrix * ccm, CSMatrix * csm) : NetWriter(stream, emx, cmx, ccm, csm) {}
    void initialize();
    void writeEdgeCluster(int cluster_index, QVector<QString> passed);
    void finish() {}
};

class MinimalNetWriter: public NetWriter
{
public:
    MinimalNetWriter(QTextStream * stream, ExpressionMatrix * emx, CorrelationMatrix * cmx,
                     CCMatrix * ccm, CSMatrix * csm) : NetWriter(stream, emx, cmx, ccm, csm) {}
    void initialize();
    void writeEdgeCluster(int cluster_index, QVector<QString> passed);
    void finish() {}
};

#endif // NET_WRITER_H
