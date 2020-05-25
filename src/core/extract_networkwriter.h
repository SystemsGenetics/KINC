#ifndef EXTRACT_NETWORKWRITER_H
#define EXTRACT_NETWORKWRITER_H
#include <QObject>
#include <QTextStream>
#include "ccmatrix_pair.h"
#include "ccmatrix.h"
#include "conditionspecificclustersmatrix.h"
#include "conditionspecificclustersmatrix_pair.h"
#include "correlationmatrix_pair.h"
#include "correlationmatrix.h"
#include "expressionmatrix.h"
#include "extract.h"



/*!
 * The parent class for writing a network to an output file.
 *
 * Child classes need only implement the initialize(), writeEdge() and
 * finish() functions.
 */
class Extract::NetworkWriter
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
    QString getEdgeSampleString(int clusterIndex) const;
    int getEdgeNumSamples(QString sample_string) const;
    float getEdgeSimilarity(int clusterIndex) const;
    float getEdgeTestValue(int clusterIndex, const QString& fullTestName) const;
    QVector<QString> getTestNames() const { return _testNames; }
    void checkStatus() const;
public:
    virtual void initialize() = 0;
    virtual void writeEdgeCluster(int clusterIndex, QVector<QString> passed) = 0;
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



/**
 * The full text network writer includes data from the CCM and CMX, and it
 * writes each conditional test from the CSM in separate columns.
 */
class Extract::FullNetworkWriter : public Extract::NetworkWriter
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
    void writeEdgeCluster(int clusterIndex, QVector<QString> passed);
    void finish() {}
};



/**
 * The tidy text network writer includes data from the CCM and CMX, and it
 * writes a single column for the conditional tests from the CSM. Edges that
 * are significant for multiple tests are duplicated for each test.
 */
class Extract::TidyNetworkWriter : public Extract::NetworkWriter
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
    void writeEdgeCluster(int clusterIndex, QVector<QString> passed);
    void finish() {}
};



/**
 * The minimal text network writer includes data from the CMX only, it does
 * not include sample strings from the CCM or conditional tests from the CSM.
 */
class Extract::MinimalNetworkWriter : public Extract::NetworkWriter
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
    void writeEdgeCluster(int clusterIndex, QVector<QString> passed);
    void finish() {}
};



/**
 * The GraphML network writer includes all information from the CMX, CCM, and
 * CSM in GraphML format. The conditional tests are written in the same manner
 * as the full text network writer.
 */
class Extract::GMLNetworkWriter : public Extract::NetworkWriter
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
    void writeEdgeCluster(int clusterIndex, QVector<QString> passed);
    void finish();
};



#endif // EXTRACT_NETWORKWRITER_H
