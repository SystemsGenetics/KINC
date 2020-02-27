#ifndef NETWORKWRITER_H
#define NETWORKWRITER_H
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
class NetworkWriter
{
public:
    NetworkWriter(
        ExpressionMatrix* emx,
        CorrelationMatrix* cmx,
        CCMatrix* ccm,
        CSMatrix* csm,
        QFile* output);
    virtual ~NetworkWriter() {}
public:
    int readNext();
    QString getEdgeGene1() const;
    QString getEdgeGene2() const;
    QString getEdgeSampleString(int cluster_index) const;
    int getEdgeNumSamples(QString sample_string) const;
    float getEdgeSimilarity(int cluster_index) const;
    float getEdgeTestValue(int cluster_index, int test_index) const;
    QVector<QString> getTestNames() const { return _testNames; }
    void checkStatus() const;
public:
    virtual void initialize() = 0;
    virtual void writeEdgeCluster(int cluster_index, QVector<QString> passed) = 0;
    virtual void finish() = 0;
protected:
    /**
     * Text stream for the output text file.
     */
    QTextStream _stream;
    /*!
     * Pointer to the input expression matrix.
     */
    ExpressionMatrix* _emx {nullptr};
    /*!
     * Pointer to the input correlation matrix and corresponding pair object.
     */
    CorrelationMatrix* _cmx {nullptr};
    CorrelationMatrix::Pair _cmxPair;
    /*!
     * Pointer to the input cluster matrix and corresponding pair object.
     */
    CCMatrix* _ccm {nullptr};
    CCMatrix::Pair _ccmPair;
    /*!
     * Pointer to the input condition specific cluster matrix and
     * corresponding pair object.
     */
    CSMatrix* _csm {nullptr};
    CSMatrix::Pair _csmPair;
    /*!
     * Pointer to the output text file.
     */
    QFile* _output {nullptr};
    /**
     * The variable that houses the test names for easy lookup and
     * the function that should be used to set the test names.
     */
    void setTestNames();
    QVector<QString> _testNames;
};



class GMLNetworkWriter: public NetworkWriter
{
private:
    /*!
     * A hash lookup for storing nodes in the network.
     */
    QHash<QString, bool> _nodes;
public:
    GMLNetworkWriter(
        ExpressionMatrix* emx,
        CorrelationMatrix* cmx,
        CCMatrix* ccm,
        CSMatrix* csm,
        QFile* output): NetworkWriter(emx, cmx, ccm, csm, output) {}
public:
    void initialize();
    void writeEdgeCluster(int cluster_index, QVector<QString> passed);
    void finish();
};



class FullNetworkWriter: public NetworkWriter
{
public:
    FullNetworkWriter(
        ExpressionMatrix* emx,
        CorrelationMatrix* cmx,
        CCMatrix* ccm,
        CSMatrix* csm,
        QFile* output): NetworkWriter(emx, cmx, ccm, csm, output) {}
public:
    void initialize();
    void writeEdgeCluster(int cluster_index, QVector<QString> passed);
    void finish() {}
};



class TidyNetworkWriter: public NetworkWriter
{
public:
    TidyNetworkWriter(
        ExpressionMatrix* emx,
        CorrelationMatrix* cmx,
        CCMatrix* ccm,
        CSMatrix* csm,
        QFile* output): NetworkWriter(emx, cmx, ccm, csm, output) {}
public:
    void initialize();
    void writeEdgeCluster(int cluster_index, QVector<QString> passed);
    void finish() {}
};



class MinimalNetworkWriter: public NetworkWriter
{
public:
    MinimalNetworkWriter(
        ExpressionMatrix* emx,
        CorrelationMatrix* cmx,
        CCMatrix* ccm,
        CSMatrix* csm,
        QFile* output): NetworkWriter(emx, cmx, ccm, csm, output) {}
public:
    void initialize();
    void writeEdgeCluster(int cluster_index, QVector<QString> passed);
    void finish() {}
};



#endif // NETWORKWRITER_H
