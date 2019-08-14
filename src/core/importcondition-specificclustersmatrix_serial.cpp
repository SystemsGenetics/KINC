#include "importcondition-specificclustersmatrix_serial.h"
#include "importcondition-specificclustersmatrix_workblock.h"
#include "importcondition-specificclustersmatrix_resultblock.h"
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
#include <iostream>

#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>

#include <iomanip>
//

/*!
*  Implements an interface to create an serial object.
*/
importCSCM::Serial::Serial(importCSCM* parent) : EAbstractAnalyticSerial(parent), _base(parent)
{
    EDEBUG_FUNC(this,parent);
}





/*!
*  Implements an interface to perfomr the desired work on a result block.
*
* @param block The work block you are working on.
*/
std::unique_ptr<EAbstractAnalyticBlock> importCSCM::Serial::execute(const EAbstractAnalyticBlock* block)
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
        if((ccmPair.index().getX() != 1 && ccmPair.index().getY() != 0) || ccmIndex != start)
        {
            ccmPair.readNext();
            cmxPair.read(ccmPair.index());
        }
        //if the first one isnt in the cluster we should not count it.
        if(ccmPair.clusterSize() == 0)
        {
            size++;
        }

        // Initialize new pvalues, one set of pvalues for each cluster.
        QVector<QVector<double>> pValues;
        pValues.resize(ccmPair.clusterSize());

        //for each cluster in the pair, run the binomial and linear regresion tests.
        for(qint32 clusterIndex = 0; clusterIndex < ccmPair.clusterSize(); clusterIndex++)
        {
            //check to see if the corrolation is high enough to test.
            if(cmxPair.at(clusterIndex) < _base->_corrthresh)
            {
                //if not skip this cluster
                continue;
            }

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
                        test(cmxPair, ccmPair, clusterIndex, testIndex, featureIndex, labelIndex, pValues);
                    }
                    //if only the feature needs testing (no sub labels)
                    else
                    {
                        test(cmxPair,ccmPair, clusterIndex, testIndex, featureIndex, 0, pValues);
                    }
                }
            }
        }

        //if the matrix is empty we dont need to keep its data.
        if(!isEmpty(pValues))
        {
            CSCMPair pair;
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
* @param testLabel The label you are testing.
*
* @param dataIndex The feature the label is in.
*
* @return True if the matrix is empty, false otherwise.
*/
int importCSCM::Serial::prepAnxData(QString testLabel, int dataIndex)
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
* @param vector The vector you want to check.
*
* @return True if the matrix is empty, false otherwise.
*/
bool importCSCM::Serial::isEmpty(QVector<QVector<double>>& vector)
{
    EDEBUG_FUNC(this, &vector);

    int index = 0;
    while(index < vector.size())
    {
        if(!vector.at(index++).isEmpty())
        {
            return false;
        }
    }
    return true;
}





/*!
*  Implements an interface to prepare the cluster catagopry count information.
*
* @param ccmPair The gene pair that we are counting the labels for.
*
* @param clusterIndex The number cluster we are in in the pair
*
* @return The number of labels in the given cluster.
*/
int importCSCM::Serial::clusterInfo(CCMatrix::Pair& ccmPair, int clusterIndex, QString label)
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
*  An interface to run the correct tests on the pair..
*
* @param ccmPair The pair thats going to be tested.
*
* @param clusterIndex The cluster to test inside the pair, each cluster is tested.
*
* @param testIndex Which test we are currently performing.
*
* @param featureIndex The current feature we are testing.
*
* @param labelIndex The sub label in the feature we are running a test on.
*
* @param pValues The two dimensional array holding all of the results from the tests.
*
* @return The test that was just conducted.
*/
int importCSCM::Serial::test(CorrelationMatrix::Pair cmxPair,
                             CCMatrix::Pair& ccmPair,
                             qint32 clusterIndex,
                             qint32& testIndex,
                             qint32 featureIndex,
                             qint32 labelIndex,
                             QVector<QVector<double>>& pValues)
{
    EDEBUG_FUNC(this, &cmxPair,&ccmPair, clusterIndex, &testIndex, featureIndex, labelIndex, &pValues);

    //get informatiopn on the mask
    clusterInfo(ccmPair, clusterIndex, _base->_features.at(featureIndex).at(labelIndex));

    //conduct the correct test based on the type of data
    switch(_base->_testType.at(featureIndex))
    {
        case CATEGORICAL :
            pValues[clusterIndex][testIndex] = binomial(_base->_alpha);
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
* @param clusterInfo The cluster to run the test on.
*
* @return Pvalue corrosponding to the test, negative if not accepted.
*/
double importCSCM::Serial::binomial(double alpha)
{
    EDEBUG_FUNC(this, alpha);

    //calculate the pvalue, the max out of the two tests
    double test1 = testOne();
    double test2 = testTwo();
    double pvalue = std::max(test1, test2);

    if(pvalue < alpha)
    {
        pvalue *= -1;
    }
    return pvalue;
}






/*!
*  Implements an interface to run the first binomial test for given data.
*
* @param clusterInfo The cluster to run the test on.
*
* @param anxInfo the column information for the feature, used to compare the label
*        with.
*
* @param label The label to compare the anxData with, its the current test label.
*
* @param numLabels The number of samples with the same label as the test label.
*
* @return Pvalue corrosponding to the test, negative if not accepted.
*/
double importCSCM::Serial::testOne()
{

    EDEBUG_FUNC(this);

    // Test #1
    // successes = number of non-category samples in the cluster
    // failures = number of non-category samples not in the cluster
    // Ho: successes >= 0.15
    // Ha: successes < 0.15
    double pvalue = 0.0;

    pvalue = gsl_cdf_binomial_P(_clusterInMask - _catInCount, .25, _clusterInMask - _catInCount);

    return pvalue;
}





/*!
*  Implements an interface to run the second binomial test for given data.
*
* @param clusterInfo The cluster to run the test on.
*
* @param anxInfo the column information for the feature, used to compare the label
*        with.
*
* @param label The label to compare the anxData with, its the current test label.
*
* @param numLabels The number of samples with the same label as the test label.
*
* @return Pvalue corrosponding to the test, negative if not accepted.
*/
double importCSCM::Serial::testTwo()
{
    EDEBUG_FUNC(this);

    // Test #2
    // successes = number of category samples in the cluster
    // failures = number of category samples out of the cluster
    // Ho: successes = 0.85
    // Ha: successes > 0.85
    //gsl_ran_binomial_pdf(k, p, n)
    double pvalue = 0.0;

    pvalue = gsl_cdf_binomial_Q(_catInCount, .75, _catCount);

    return pvalue;
}





/*!
*  Implements an interface to run the regresion test for given data, the regresion line is
*  genex vs geney vs label data.
*
* @param anxInfo Annotation matrix information.
*
* @param ccmPair Cluster matrix pair.
*
* @param clusterIndex The index of the cluster, used to get the right info from the cluster pair
*
* @return Pvalue corrosponding to the test, negative if not accepted.
*/
double importCSCM::Serial::regresion(QVector<QString> &anxInfo, CCMatrix::Pair& ccmPair, int clusterIndex)
{
    EDEBUG_FUNC(this, &anxInfo, &ccmPair, clusterIndex);

    //temp containers
    QVector<double> labelInfo;

    //regression model containers
    double chisq;
    gsl_matrix *X, *cov;
    gsl_vector *Y, *C;
    double pValue = 0.0;

    //allocate a matrix to hold the predictior variables, in this cas the gene expression data
    X = gsl_matrix_alloc (_clusterInMask, 3);

    //allocate a vector to hold observation data, in this case the data corrosponding ot the features
    Y = gsl_vector_alloc (_clusterInMask);

    //allocate a vector and matrix for the slop info
    C = gsl_vector_alloc (3);
    cov = gsl_matrix_alloc (_clusterInMask, 3);

    //Read in the gene pairs expression information
    ExpressionMatrix::Gene geneX(_base->_emx);
    ExpressionMatrix::Gene geneY(_base->_emx);

    geneX.read(ccmPair.index().getX());
    geneY.read(ccmPair.index().getY());

    //look through all the samples in the mask
    for(qint32 i = 0, j = 0; i < _base->_emx->sampleSize(); i++)
    {
        //if the sample label matches with the given label
        if(ccmPair.at(clusterIndex, i) == 1)
        {
            //add emx data to the predictions if the sample is in the cluster
            gsl_matrix_set(X, j, 0, static_cast<double>(1));
            gsl_matrix_set(X, j, 1, static_cast<double>(geneX.at(i)));
            gsl_matrix_set(X, j, 2, static_cast<double>(geneY.at(i)));

            //add the feature data into the observations vector
            if(!labelInfo.contains(anxInfo.at(i).toInt()))
            {
                labelInfo.append(anxInfo.at(i).toInt());
            }
            for(int k = 0; k < labelInfo.size(); k++)
            {
                if(labelInfo.at(k) == anxInfo.at(i).toInt())
                {
                    gsl_vector_set(Y, j, k + 1);
                }
            }
            j++;
        }
    }

    //create the workspace for the gnu scientific library to work in
    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (_clusterInMask, 3);

    //regrassion calculation
    gsl_multifit_linear(X, Y, C, cov, &chisq, work);


    /* vvv For Debugging vvv */
    std::cout << std::endl << "coefficients:" << std::endl;
    for (int j = 0; j < 3; j++)
    {
        std::cout << "c" << j << " = " << std::setprecision(9);
        std::cout << gsl_vector_get(C, j) << std::endl;
    }
    std::cout << "R squared: " << chisq << std::endl << std::endl;

    std::cout << std::endl;
    std::cout << "expected <=> predicted" << std::endl;
    for (int i = 0 ; i < _clusterInMask; i++)
    {
        double r = gsl_vector_get(C, 0);
        r += gsl_matrix_get(X, i, 1) * gsl_vector_get(C, 1);
        r += gsl_matrix_get(X, i, 2) * gsl_vector_get(C, 2);
        std::cout << gsl_vector_get(Y, i) << " <=> " << std::setprecision(9) << r << std::endl << std::endl;
    }

    for(int i = 0; i < _clusterInMask; i++)
    {
        //print out the predictor variables matrix
        std::cout << "Matrix of predictor values:" <<
        std::cout << gsl_matrix_get(X, i, 0) << "\t";
        std::cout << gsl_matrix_get(X, i, 1) << "\t";
        std::cout << gsl_matrix_get(X, i, 2) << std::endl;

        //print out the observation variables matrix
        std::cout << gsl_vector_get(Y, i) << std::endl;
    }
    /* ^^^ For Debugging ^^^ */

    //grab the variable
    pValue = gsl_vector_get(C, 1);

    //free all of the data
    gsl_matrix_free(X);
    gsl_vector_free(Y);
    gsl_matrix_free(cov);
    gsl_vector_free(C);
    gsl_multifit_linear_free(work);

    //return the slope of the line
    return pValue;
}
