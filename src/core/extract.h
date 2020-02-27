#ifndef EXTRACT_H
#define EXTRACT_H
#include <ace/core/core.h>

#include "ccmatrix_pair.h"
#include "ccmatrix.h"
#include "correlationmatrix_pair.h"
#include "correlationmatrix.h"
#include "expressionmatrix.h"
#include "conditionspecificclustersmatrix.h"
#include "conditionspecificclustersmatrix_pair.h"
#include "networkwriter.h"


/*!
 * This class implements the extract analytic. This analytic is very similar to
 * the export correlation matrix analytic, except for a few differences: (1) this
 * analytic uses a slightly different format for the text file, (2) this analytic
 * can apply a correlation threshold, and (3) this analytic can optionally write
 * a GraphML file. The key difference is that this analytic "extracts" a network
 * from the correlation matrix and writes an edge list rather than a correlation
 * list.
 */
class Extract : public EAbstractAnalytic
{
    Q_OBJECT
public:
    class Input;
public:
    virtual int size() const override final;
    virtual void process(const EAbstractAnalyticBlock* result) override final;
    virtual EAbstractAnalyticInput* makeInput() override final;
    virtual void initialize();
private:
    /*!
     * Defines the output formats this analytic supports.
     */
    enum class OutputFormat
    {
        /*!
         * Text format
         */
        Text
        /*!
         * Minimal format (does not use CCM)
         */
        ,Minimal
        /*!
         * GraphML format
         */
        ,GraphML
        /*!
         * GraphML format
         */
        ,Tidy
    };
private:
    /**
     * Workspace variables to write to the output file
     */
    QTextStream _stream;
    CCMatrix::Pair _ccmPair;
    CorrelationMatrix::Pair _cmxPair;
    CSMatrix::Pair _csmPair;
    /*!
     * Pointer to the input expression matrix.
     */
    ExpressionMatrix* _emx {nullptr};
    /*!
     * Pointer to the input cluster matrix.
     */
    CCMatrix* _ccm {nullptr};
    /*!
     * Pointer to the input correlation matrix.
     */
    CorrelationMatrix* _cmx {nullptr};
    /*!
     * Pointer to the input condition specific cluster matrix.
     */
    CSMatrix* _csm {nullptr};
    /*!
     * Pointer to the annotation file.
     */
    QFile* _amx {nullptr};
    /*!
     * The output format to use.
     */
    OutputFormat _outputFormat {OutputFormat::Text};
    /*!
     * Pointer to the output text file.
     */
    QFile* _output {nullptr};
    /*!
     * The minimum (absolute) correlation threshold.
     */
    float _minCorrelation {0.85f};
    /*!
     * The maximum (absolute) correlation threshold.
     */
    float _maxCorrelation {1.00f};

    /*!
     * Condition-Specific Cluster Matrix p-value and r-squared filter input.
     */
    QString _csmPValueFilter {""};
    QString _csmRSquareFilter {""};


    /*!
     * An instance of a NetworkWriter class that ensures that
     * the same edges are always written in any file format.
     */
    NetworkWriter * _networkWriter {nullptr};

    /*!
     * Stores the names of the condition-specific testing
     * that was performed.
     */
    QVector<QString> _testNames;

    /*!
     * An associative array where the key is the test name
     * and the value is a pair where the first element of
     * the pair is the logical comparison value (e.g. gt, lt)
     * and the second element is the value for comparing.
     */
    QHash<QString, QPair<QString, float>> _filters;

    /*!
     * Sets the _filters element.
     */
    void setFilters(QString input_filters, QString type);

    /*!
     * Performs filtering of a cluster of the current edge.
     */
    QVector<QString> filterEdge(int cluster_index);
};



#endif
