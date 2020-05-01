#ifndef ConditionalTest_H
#define ConditionalTest_H
#include <ace/core/core.h>
#include "conditionspecificclustersmatrix.h"
#include "ccmatrix.h"
#include "correlationmatrix.h"
#include "expressionmatrix.h"



class ConditionalTest : public EAbstractAnalytic
{
    Q_OBJECT
public:
    class Input;
    class Serial;
    class WorkBlock;
    class ResultBlock;

    enum TESTTYPE
    {
        CATEGORICAL = 0,
        ORDINAL,
        QUANTITATIVE,
        UNKNOWN,
        NONE
    };

    struct Pair
    {
        /*!
         * The pairwise index of the pair.
         */
        Pairwise::Index index;
        /*!
         * The p values for each cluster in a pair.
         */
        QVector<QVector<double>> pValues;
        /*!
         * The r^2 values for each cluster in a pair.
         */
        QVector<QVector<double>> r2;
    };

    virtual int size() const override final;
    virtual void process(const EAbstractAnalyticBlock* result) override final;
    virtual std::unique_ptr<EAbstractAnalyticBlock> makeWork(int index) const override final;
    virtual std::unique_ptr<EAbstractAnalyticBlock> makeWork() const override final;
    virtual std::unique_ptr<EAbstractAnalyticBlock> makeResult() const override final;
    virtual EAbstractAnalyticSerial* makeSerial() override final;
    virtual EAbstractAnalyticInput* makeInput() override final;
    virtual void initialize() override final;
    virtual void initializeOutputs() override final;

    /*!
     * Reads in the annotation matrix populating the metadata when its done.
     */
    void readInAMX(QVector<QVector<QString>>& amxdata,
                   QVector<QVector<QVariant>>& data,
                   QVector<TESTTYPE>& dataTestType);
    void configureTests(QVector<TESTTYPE>& dataTestType);
    int max(QVector<qint32> &counts) const;
    /*!
     * Test overrides.
     */
    void Test();
    void override();
    /*!
     * Crates a dtring of test names, delimited by a colon.
     */
    QString testNames();

    void initialize(qint32 &maxClusterSize, qint32 &subHeaderSize,QVector<QVector<QString>> &amxData, QVector<TESTTYPE> &testType, QVector<QVector<QVariant>> &data);

    void rearrangeSamples();

private:
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
     * Pointer to the input annotation matrix file.
     */
    QFile* _amx {nullptr};
    /*!
     * Pointer to the output cluster annotation matrix.
     */
    CSMatrix* _out {nullptr};
    /*!
     * User provided features not to test.
     */
    QString _Testing{""};
    QVector<QString> _Test;
    /*!
     * User provided test type overrides.
     */
    QString _testOverride{""};
    QVector<QVector<QString>> _override;
    /*!
     * Assosiated stream for the annotation matrix input file.
     */
    QTextStream _stream;
    /*!
     * Number of lines in the input annotation matrix file.
     */
    qint32 _amxNumLines {0};
    /*!
     * The number of pairs to process in each work block.
     */
    int _workBlockSize {0};
    /*!
     * Current start index used to create work blocks.
     */
    mutable qint64 _workBlockStart {0};
    /*!
     * Annotation matrix data.
     */
    QVector<QVector<QString>> _features;
    QVector<QVector<QVariant>> _data;
    QVector<TESTTYPE> _testType;
    int _numTests {0};
    qint32 _geneSize {0};
    qint32 _sampleSize {0};
    QString _delimiter = "tab";
    QString _missing = "NA";
    /*!
     * Cluster information
     */
    QVector<QVector<Pairwise::Index>> _clusters;
};



#endif
