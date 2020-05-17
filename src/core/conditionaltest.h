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

    enum TestType
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

public:
    virtual int size() const override final;
    virtual void process(const EAbstractAnalyticBlock* result) override final;
    virtual std::unique_ptr<EAbstractAnalyticBlock> makeWork(int index) const override final;
    virtual std::unique_ptr<EAbstractAnalyticBlock> makeWork() const override final;
    virtual std::unique_ptr<EAbstractAnalyticBlock> makeResult() const override final;
    virtual EAbstractAnalyticSerial* makeSerial() override final;
    virtual EAbstractAnalyticInput* makeInput() override final;
    virtual void initialize() override final;
    virtual void initializeOutputs() override final;

public:
    int max(QVector<qint32> &counts) const;
    /*!
     * Crates a string of test names, delimited by a colon.
     */
    QString testNames();

private:

    /*!
     * Test overrides.
     */
    void setUserTests();
    void setUserTestTypes();

    /*!
     * Reads in the annotation matrix populating the metadata when its done.
     */
    void setFeatures();
    void setTestTypes();
    void setTestCount();
    void setData();
    void orderLabelsBySample();
    void setNumTests();

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
    QString _userTestsStr{""};
    QVector<QString> _userTests;
    /*!
     * User provided test type overrides.
     */
    QString _userTestTypesStr{""};
    QVector<QVector<QString>> _userTestTypes;
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
    QVector<TestType> _testType;
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
