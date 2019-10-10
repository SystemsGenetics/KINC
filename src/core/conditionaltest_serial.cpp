#include "conditionaltest_serial.h"
#include "conditionaltest_workblock.h"
#include "conditionaltest_resultblock.h"
#include "ccmatrix_pair.h"
#include "correlationmatrix_pair.h"
#include "correlationmatrix.h"
#include "correlationmatrix_pair.h"
#include "expressionmatrix_gene.h"
#include <ace/core/elog.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
//

/*!
*  Implements an interface to create a serial object.
*/
ConditionalTest::Serial::Serial(ConditionalTest* parent) : EAbstractAnalyticSerial(parent), _base(parent)
{
    EDEBUG_FUNC(this,parent);
}





/*!
*  Implements an interface to perform the desired work on a result block.
*
* @param block The work block you are working on.
*
* @return The populated result block.
*/
std::unique_ptr<EAbstractAnalyticBlock> ConditionalTest::Serial::execute(const EAbstractAnalyticBlock* block)
{
    EDEBUG_FUNC(this, block);

    //crate the work and result blocks
    const WorkBlock* workBlock {block->cast<WorkBlock>()};
    ResultBlock* resultBlock {new ResultBlock(workBlock->index(), _base->_numTests, workBlock->startpair())};

    // Create iterators for the CCM data object.
    CCMatrix::Pair ccmPair = CCMatrix::Pair(_base->_ccm);
    CorrelationMatrix::Pair cmxPair = CorrelationMatrix::Pair(_base->_cmx);

    // Iterate through the pairs the workblock is assigned to.
    qint64 start = workBlock->startpair();
    qint64 size = workBlock->size();
    Pairwise::Index index(workBlock->start());

    //iterate through each pair in the matrix
    for ( qint64 ccmIndex = start; ccmIndex < start + size; ccmIndex++ )
    {
        //reads the first value in the ccm

        if(ccmIndex == start)
        {
            ccmPair.read(index);
            cmxPair.read(ccmPair.index());
        }
        if((ccmPair.index().getX() != 1 && ccmPair.index().getY()!= 0) || ccmIndex != start)
        {
            ccmPair.readNext();
            cmxPair.read(ccmPair.index());
        }
        //if the first one isnt in the cluster we should not count it.
        if(ccmPair.clusterSize() == 0)
        {
            size++;
            continue;
        }

        // Initialize new pvalues, one set of pvalues for each cluster.
        QVector<QVector<double>> pValues;
        pValues.resize(ccmPair.clusterSize());

        //for each cluster in the pair, run the binomial and linear regresion tests.
        for(qint32 clusterIndex = 0; clusterIndex < ccmPair.clusterSize(); clusterIndex++)
        {
            //resize for room for each test.
            pValues[clusterIndex].resize(_base->_numTests);

            for(qint32 featureIndex = 0, testIndex = 0; featureIndex < _base->_features.size(); featureIndex++)
            {
                if(_base->_testType.at(featureIndex) == NONE || _base->_testType.at(featureIndex) == UNKNOWN)
                {
                    continue;
                }

                for(qint32 labelIndex = 0; labelIndex < _base->_features.at(featureIndex).size(); labelIndex++)
                {

                    prepAnxData(_base->_features.at(featureIndex).at(labelIndex), featureIndex);
                    //if there are sub labels to test for the feature
                    if(_base->_features.at(featureIndex).size() > 1)
                    {
                        if(labelIndex == 0)
                        {
                            labelIndex = 1;
                        }
                        test(ccmPair, clusterIndex, testIndex, featureIndex, labelIndex, pValues);
                    }
                    //if only the feature needs testing (no sub labels)
                    else
                    {
                        test(ccmPair, clusterIndex, testIndex, featureIndex, 0, pValues);
                    }
                }
            }
        }
        if(!isEmpty(pValues))
        {
            CSMPair pair;
            pair.pValues = pValues;
            pair.x_index = ccmPair.index().getX();
            pair.y_index = ccmPair.index().getY();
            resultBlock->append(pair);
        }
    }
    return std::unique_ptr<EAbstractAnalyticBlock>(resultBlock);
}






/*!
*  Implements an interface to prerpare the annotation matrix data for testing.
*
* @param testLabel The label you are testing on.
*
* @param dataIndex The feature the label is part of.
*
* @return The number of samples in total of the test label.
*/
int ConditionalTest::Serial::prepAnxData(QString testLabel, int dataIndex)
{
    EDEBUG_FUNC(this, testLabel, dataIndex);

    //get the needed data fro the comparison
    _catCount = 0;
    _anxData.resize(_base->_data.at(dataIndex).size());
    //populate array with annotation data
    for(int j = 0; j < _base->_data.at(dataIndex).size(); j++)
    {
        _anxData[j] = _base->_data.at(dataIndex).at(j).toString();
        //if data is the same as the test label add one to the catagory counter
        if(_anxData[j] == testLabel)
        {
            _catCount++;
        }
    }
    return _catCount;
}




/*!
*  Implements an interface to check to see if a matrix is empty.
*
* @param vector The matrix you want to check.
*
* @return True if the matrix is empty, false otherwise.
*/
bool ConditionalTest::Serial::isEmpty(QVector<QVector<double>>& matrix)
{
    EDEBUG_FUNC(this, &matrix);

    int index = 0;
    while(index < matrix.size())
    {
        if(!matrix.at(index++).isEmpty())
        {
            return false;
        }
    }
    return true;
}





/*!
*  Implements an interface to prepare the cluster catagory count information.
*
* @param ccmPair The gene pair that we are counting the labels for.
*
* @param clusterIndex The number cluster we are in in the pair
*
* @return The number of labels in the given cluster.
*/
int ConditionalTest::Serial::clusterInfo(CCMatrix::Pair& ccmPair, int clusterIndex, QString label)
{
    _catCount = _clusterInMask = _catInCount = 0;

    //look through all the samples in the mask
    for(qint32 i = 0; i < _base->_emx->sampleSize(); i++)
    {
        //if the sample label matches with the given label'
        if(_anxData.at(i) == label)
        {
            _catCount++;
        }
        if(ccmPair.at(clusterIndex, i) == 1)
        {
            _clusterInMask++;
            if(_anxData.at(i) == label)
            {
                   _catInCount++;
            }
        }
    }
    return _catInCount;
}




/*!
*  An interface to choose and run the correct tests on the pair.
*
* @param ccmPair The pair thats going to be tested.
*
* @param clusterIndex The cluster to test inside the pair, each cluster is
*        tested.
*
* @param testIndex Which test we are currently performing.
*
* @param featureIndex The current feature we are testing.
*
* @param labelIndex The label in the feature we are running a test on.
*
* @param pValues The two dimensional array holding all of the results from the
*        tests.
*
* @return The test that was just conducted.
*/
int ConditionalTest::Serial::test(CCMatrix::Pair& ccmPair,
                             qint32 clusterIndex,
                             qint32& testIndex,
                             qint32 featureIndex,
                             qint32 labelIndex,
                             QVector<QVector<double>>& pValues)
{
    EDEBUG_FUNC(this,&ccmPair, clusterIndex, &testIndex, featureIndex, labelIndex, &pValues);

    //get informatiopn on the mask
    clusterInfo(ccmPair, clusterIndex, _base->_features.at(featureIndex).at(labelIndex));

    //conduct the correct test based on the type of data
    switch(_base->_testType.at(featureIndex))
    {
        case CATEGORICAL :
            pValues[clusterIndex][testIndex] = binomial();
            testIndex++;
        break;
        case ORDINAL :
            pValues[clusterIndex][testIndex] = regresion(_anxData, ccmPair, clusterIndex);
            testIndex++;
        break;
        case QUANTATATIVE :
            pValues[clusterIndex][testIndex] = regresion(_anxData, ccmPair, clusterIndex);
            testIndex++;
        break;
        default:; //quash a compiler warning
    }
    return _base->_testType.at(featureIndex);
}





/*!
*  Implements an interface to run the binomial test for given data.
*
* @param alpha The threshold for keeping the data.
*
* @return Pvalue corrosponding to the test, negative if not accepted.
*/
double ConditionalTest::Serial::binomial()
{
    EDEBUG_FUNC(this);

    //calculate the pvalue, the max out of the two tests
    double test1 = testOne();
    double test2 = testTwo();
    double pvalue = std::max(test1, test2);

    return pvalue;
}






/*!
*  Implements an interface to run the first binomial test for given data.
*
* @return Pvalue corrosponding to the test.
*/
double ConditionalTest::Serial::testOne()
{

    EDEBUG_FUNC(this);

    // Test #1
    // successes = number of non-category samples in the cluster
    // failures = number of non-category samples not in the cluster
    // Ho: successes >= 0.15
    // Ha: successes < 0.15

    return gsl_cdf_binomial_P(_clusterInMask - _catInCount, 1 - _base->_probabilitySuccess, _base->_emx->sampleSize() - _catCount);
}





/*!
*  Implements an interface to run the second binomial test for given data.
*
* @return Pvalue corrosponding to the test.
*/
double ConditionalTest::Serial::testTwo()
{
    EDEBUG_FUNC(this);

    // Test #2
    // successes = number of category samples in the cluster
    // failures = number of category samples out of the cluster
    // Ho: successes = 0.85
    // Ha: successes > 0.85

    return gsl_cdf_binomial_Q(_catInCount, _base->_probabilitySuccess, _catCount - _catInCount);
}





/*!
*  Implements an interface to run the regresion test for given data, the
*  regresion line is genex vs geney vs label data.
*
* @param anxInfo Annotation matrix information.
*
* @param ccmPair Cluster matrix pair.
*
* @param clusterIndex The index of the cluster, used to get the right info
*        from the cluster pair
*
* @return Pvalue corrosponding to the test.
*/
double ConditionalTest::Serial::regresion(QVector<QString> &anxInfo, CCMatrix::Pair& ccmPair, int clusterIndex)
{
    EDEBUG_FUNC(this, &anxInfo, &ccmPair, clusterIndex);

    Q_UNUSED(anxInfo);
    Q_UNUSED(ccmPair);
    Q_UNUSED(clusterIndex);

    return -1;
}
