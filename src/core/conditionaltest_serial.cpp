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

        if ( ccmIndex == start )
        {
            ccmPair.read(index);
            cmxPair.read(ccmPair.index());
        }
        if ( (ccmPair.index().getX() != 1 && ccmPair.index().getY()!= 0) || ccmIndex != start )
        {
            ccmPair.readNext();
            cmxPair.read(ccmPair.index());
        }
        //if the first one isnt in the cluster we should not count it.
        if ( ccmPair.clusterSize() == 0 )
        {
            size++;
            continue;
        }

        // Initialize new pvalues, one set of pvalues for each cluster.
        QVector<QVector<double>> pValues;
        pValues.resize(ccmPair.clusterSize());

        //for each cluster in the pair, run the binomial and linear regresion tests.
        for ( qint32 clusterIndex = 0; clusterIndex < ccmPair.clusterSize(); clusterIndex++ )
        {
            //resize for room for each test.
            pValues[clusterIndex].resize(_base->_numTests);

            for ( qint32 featureIndex = 0, testIndex = 0; featureIndex < _base->_features.size(); featureIndex++ )
            {
                if ( _base->_testType.at(featureIndex) == NONE || _base->_testType.at(featureIndex) == UNKNOWN )
                {
                    continue;
                }
                if(_base->_testType.at(featureIndex) == QUANTATATIVE || _base->_testType.at(featureIndex) == ORDINAL)
                {
                    prepAnxData(_base->_features.at(featureIndex).at(0), featureIndex, _base->_testType.at(featureIndex));
                    test(ccmPair, clusterIndex, testIndex, featureIndex, 0, pValues);
                }
                else if (_base->_testType.at(featureIndex) == CATEGORICAL)
                {
                    for ( qint32 labelIndex = 0; labelIndex < _base->_features.at(featureIndex).size(); labelIndex++ )
                    {

                        prepAnxData(_base->_features.at(featureIndex).at(labelIndex), featureIndex, _base->_testType.at(featureIndex));
                        //if there are sub labels to test for the feature
                        if ( _base->_features.at(featureIndex).size() > 1 )
                        {
                            if ( labelIndex == 0 )
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
        }
        if ( !isEmpty(pValues) )
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
int ConditionalTest::Serial::prepAnxData(QString testLabel, int dataIndex, TESTTYPE testType)
{
    EDEBUG_FUNC(this, testLabel, dataIndex);

    //get the needed data fro the comparison
    _catCount = 0;
    _anxData.resize(_base->_data.at(dataIndex).size());
    //populate array with annotation data
    for ( int j = 0; j < _base->_data.at(dataIndex).size(); j++ )
    {
        _anxData[j] = _base->_data.at(dataIndex).at(j).toString();
        //if data is the same as the test label add one to the catagory counter
        if ( testType == _base->CATEGORICAL && _anxData[j] == testLabel )
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
        if ( !matrix.at(index++).isEmpty() )
        {
            return false;
        }
    }
    return true;
}





/*!
*  Implements an interface to prepare the cluster category count information.
*
* @param ccmPair The gene pair that we are counting the labels for.
*
* @param clusterIndex The number cluster we are in in the pair
*
* @return The number of labels in the given cluster.
*/
int ConditionalTest::Serial::clusterInfo(CCMatrix::Pair& ccmPair, int clusterIndex, QString label, TESTTYPE testType)
{
    _catCount = _clusterSize = _catInCluster = 0;

    // Look through all the samples in the mask.
    for ( qint32 i = 0; i < _base->_emx->sampleSize(); i++ )
    {
        // If the sample label matches with the given label.
        if ( testType == _base->CATEGORICAL && _anxData.at(i) == label )
        {
            _catCount++;
        }
        if ( ccmPair.at(clusterIndex, i) == 1 )
        {
            _clusterSize++;
            if ( testType == _base->CATEGORICAL && _anxData.at(i) == label )
            {
                   _catInCluster++;
            }
        }
    }
    return _catInCluster;
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
    clusterInfo(ccmPair, clusterIndex, _base->_features.at(featureIndex).at(labelIndex), _base->_testType.at(featureIndex));

    //conduct the correct test based on the type of data
    switch(_base->_testType.at(featureIndex))
    {
        case CATEGORICAL :
            pValues[clusterIndex][testIndex] = hypergeom(ccmPair, clusterIndex);
            testIndex++;
        break;
        case ORDINAL :
            pValues[clusterIndex][testIndex] = regresion(_anxData, ccmPair, clusterIndex, ORDINAL);
            testIndex++;
        break;
        case QUANTATATIVE :
            pValues[clusterIndex][testIndex] = regresion(_anxData, ccmPair, clusterIndex, QUANTATATIVE);
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

    // Calculate the pvalue, the max out of the two tests.
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

    // Number of trials, n: is the number of non-category samples.
    int n =  _base->_emx->sampleSize() - _catCount;

    // Number of successes, k: is the number of non-category samples in the cluster.
    int k = _clusterSize - _catInCluster;

    double p = 1 - _base->_probabilitySuccess;

    // The gsl_cdf_binomial_P function uses the lower-tail of the CDF.
    // gives the probablity of the variate taking a value less than k.
    double pvalue = gsl_cdf_binomial_P(k, p, n);
    return pvalue;
}





/*!
*  Implements an interface to run the second binomial test for given data.
*
* @return Pvalue corrosponding to the test.
*/
double ConditionalTest::Serial::testTwo()
{
    EDEBUG_FUNC(this);

    // Test #2:
    // successes = number of category samples in the cluster
    // failures = number of category samples out of the cluster
    // Ho: successes = 0.85
    // Ha: successes > 0.85

    // Number of trials, n: is the number of category samples.
    int n = _catCount;

    // Number of successes, k: is the number of category samples in the cluster.
    int k = _catInCluster;

    double p = _base->_probabilitySuccess;

    // The gsl_cdf_binomial_Q function uses the upper-tail of the CDF.
    // gives the probablity of the variate taking a value greater than k.
    double pvalue = gsl_cdf_binomial_Q(k, p, n);
    return pvalue;
}






/*!
*  Implements an interface to run the first binomial test for given data.
*
* @return Pvalue corrosponding to the test.
*/
double ConditionalTest::Serial::hypergeom(CCMatrix::Pair& ccmPair, int clusterIndex)
{
    EDEBUG_FUNC(this);
    // QString samples = ccmPair.toString();

    // We use the hypergeometric distribution because the samples are
    // selected from the population for membership in the cluster without
    // replacement.

    // If a population contains n_1 elements of “type 1” and n_2 elements of
    // “type 2” then the hypergeometric distribution gives the probability
    // of obtaining k elements of “type 1” in t samples from the population
    // without replacement.

    // Population contains n1 elements of Type 1.
    int n1 = _catCount;
    // Population contains n2 elements of Type 2.
    int n2 = _base->_emx->sampleSize() - _catCount;
    // k elements of Type 1 were selected.
    int k = _catInCluster;
    // t total elements were selected.
    int t = _clusterSize;

    // If n1 == k we will always get a zero because we've
    // reached the end of the distribution, so the Ho that
    // X > k is always 0.  This happens if the cluster is 100%
    // comprised of the category we're looking for.
    if (k == n1) {
      return 1;
    }

    // If our dataset is large, the power to detect the effect
    // size increases, resulting in potentially insignificant
    // proportions having signficant p-values. Using Cohen's H
    // to set a large-effect size (e.g. 0.8) with a sig.level of
    // 0.001 and a power of 0.95 we need at least 31 samples.
    // So, we'll perform a jacknife resampling of our data
    // to calculate an average proportion of 31 samples
    if (t > 31)
    {
        // Initialize the uniform random number generator.
        const gsl_rng_type * T;
        gsl_rng * r;
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);

        // Holds the jacknife average proportion.
        int jkap = 0;

        // To perform the Jacknife resampling we will
        // perform 30 iterations (central limit thereom)
        int in = 30;
        for (int i = 0; i < in; i++) {

            // Keeps track of the number of successes, and the
            // new proportion average.
            int ns = 0;

            // Generate 31 random numbers between 0 and the size of
            // the sample string.  We will use these numbers to
            // randomly select a sample and if it is a 1 we consider
            // it a success.
            for (int j = 0; j < 31; j++)
            {
                int u = static_cast<int>(gsl_rng_uniform(r) * _base->_emx->sampleSize());
                if (ccmPair.at(clusterIndex, u) == 1)
                {
                    ns = ns + 1;
                }
            }
            jkap += ns;
        }
        jkap = jkap/in;
        gsl_rng_free(r);

        // Now reset the sample size and the proporiton of success.
        k = jkap;
        t = 31;
    }

    // The gsl_cdf_hypergeometric_Q function uses the upper-tail of the CDF.
    double pvalue = gsl_cdf_hypergeometric_Q(k, n1, n2, t);
    return pvalue;
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
double ConditionalTest::Serial::regresion(QVector<QString> &anxInfo, CCMatrix::Pair& ccmPair, int clusterIndex, TESTTYPE testType)
{
    EDEBUG_FUNC(this, &anxInfo, &ccmPair, clusterIndex);

    //temp containers
    QVector<double> labelInfo;

    //regression model containers
    double chisq;
    gsl_matrix *X, *cov;
    gsl_vector *Y, *C;
    double pValue = 0.0;

    //allocate a matrix to hold the predictior variables, in this cas the gene
    //expression data
    X = gsl_matrix_alloc (_clusterSize, 3);

    //allocate a vector to hold observation data, in this case the data
    //corrosponding ot the features
    Y = gsl_vector_alloc (_clusterSize);

    //allocate a vector and matrix for the slop info
    C = gsl_vector_alloc (3);
    cov = gsl_matrix_alloc (3, 3);

    //Read in the gene pairs expression information
    ExpressionMatrix::Gene geneX(_base->_emx);
    ExpressionMatrix::Gene geneY(_base->_emx);

    geneX.read(ccmPair.index().getX());
    geneY.read(ccmPair.index().getY());

    //look through all the samples in the mask
    for(int i = 0, j = 0; i < _base->_emx->sampleSize(); i++)
    {
        //if the sample label matches with the given label
        if(ccmPair.at(clusterIndex, i) == 1)
        {
            //add emx data to the predictions if the sample is in the cluster
            gsl_matrix_set(X, j, 0, static_cast<double>(1)); //for the intercept
            gsl_matrix_set(X, j, 1, static_cast<double>(geneX.at(i)));
            gsl_matrix_set(X, j, 2, static_cast<double>(geneY.at(i)));

            if(testType == ORDINAL)
            {
                //convert the observation data into a "design vector"
                //each unique number being a ssigned a unique integer.
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
            }
            else
            {
                gsl_vector_set(Y, j, anxInfo.at(i).toFloat());
            }
            j++;
        }
    }

    //create the workspace for the gnu scientific library to work in
    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (_clusterSize, 3);

    //regrassion calculation
    gsl_multifit_linear(X, Y, C, cov, &chisq, work);

    //From here we use the F-tests to calculate the p-Value for the entire model
    pValue = fTest(chisq, X, cov, C);

    //free all of the data
    gsl_matrix_free(X);
    gsl_vector_free(Y);
    gsl_matrix_free(cov);
    gsl_vector_free(C);
    gsl_multifit_linear_free(work);

    //return the slope of the line
    return pValue;
}





/*!
*  Implements an interface to run the F test for a linear regresion model.
*
* @param chisq Sum of the squares of residuals from gsl_multifit_linear.
*
* @param X Matrix of predictor variables used in gsl_multifit_linear.
*
* @param dov Dovariance matrix poulated by gsl_multifit_linear.
*
* @param C Observation vector used in gsl_multifit_linear.
*
* @return Pvalue corrosponding to the test.
*/
double ConditionalTest::Serial::fTest(double chisq, gsl_matrix* X, gsl_matrix* cov, gsl_vector* C)
{
    // gsl_cdf_fdist_P(double x, double nu1, double nu2)
    // x  : F statistic = Mean of squares (model) / Mean of squares (error)
    //                  = chisq / sum(Y - Ypred) for all samples in cluster
    // nu1: Degrees of freedom (Model) = 2 - 1 = 1
    // nu2: Degrees of freedom (Error) = _clusterInMask - 2

    //Degrees of freedom (Model)
    int DFM = 2 - 1;

    //Degrees of freedom (Error)
    int DFE = _clusterSize - 2;

    //Mean of square (Model)
    double MSM = chisq / DFM;

    //Mean of squares (Error)
    double SSE = 0.0;
    gsl_vector *testPoint;
    testPoint = gsl_vector_alloc (3);
    for(int i = 0; i < _clusterSize; i++)
    {
        double sumSq = 0.0, stDev = 0.0;
        gsl_vector_set(testPoint, 0, gsl_matrix_get(X, i, 0));
        gsl_vector_set(testPoint, 1, gsl_matrix_get(X, i, 1));
        gsl_vector_set(testPoint, 2, gsl_matrix_get(X, i, 2));
        gsl_multifit_linear_est(testPoint, C, cov, &sumSq, &stDev);
        SSE += sumSq;
    }
    double MSE = SSE / DFE;

    //F statistic
    double fStat = MSM / MSE;

    //Calc F test
    double pValue = gsl_cdf_fdist_P(fStat, DFM, DFE);

    return pValue;
}
