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

    Pairwise::Index _index;

    QVector<QString> _sampleStrings;
    QVector<int> _numSamples;
    void setEdgeSampleStrings();
    void setTestNames();

    QVector<QString> _testNames;

public:

    void setEdge(Pairwise::Index cmx_index);
    QString getEdgeGene1();
    QString getEdgeGene2();
    QString getEdgeSampleString(int cluster_index);
    int getEdgeNumSamples(int cluster_index);
    float getEdgeSimilarity(int cluster_index);
    float getEdgeTestValue(int cluster_index, int test_index);
    QVector<QString> getTestNames() { return _testNames; }

    NetWriter(QTextStream * stream, ExpressionMatrix * emx,
                      CorrelationMatrix * cmx, CCMatrix * ccm,
                      CSMatrix * csm);
    virtual ~NetWriter() {}
    virtual void initialize();
    virtual void writeEdgeCluster(int cluster_index);
    virtual void finish();
};

class GMLNetWriter: public NetWriter
{
public:
    GMLNetWriter(QTextStream * stream, ExpressionMatrix * emx, CorrelationMatrix * cmx,
                 CCMatrix * ccm, CSMatrix * csm) : NetWriter(stream, emx, cmx, ccm, csm) {}
    void initialize();
    void writeEdgeCluster(int cluster_index);
    void finish();
};

class FullNetWriter: public NetWriter
{
public:
    FullNetWriter(QTextStream * stream, ExpressionMatrix * emx, CorrelationMatrix * cmx,
                  CCMatrix * ccm, CSMatrix * csm) : NetWriter(stream, emx, cmx, ccm, csm) {}
    void initialize();
    void writeEdgeCluster(int cluster_index);
    void finish() {}
};

class MinimalNetWriter: public NetWriter
{
public:
    MinimalNetWriter(QTextStream * stream, ExpressionMatrix * emx, CorrelationMatrix * cmx,
                     CCMatrix * ccm, CSMatrix * csm) : NetWriter(stream, emx, cmx, ccm, csm) {}
    void initialize();
    void writeEdgeCluster(int cluster_index);
    void finish() {}
};

#endif // NET_WRITER_H
