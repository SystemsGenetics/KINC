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



/*!
 * Create a serial object.
 */
ConditionalTest::Serial::Serial(ConditionalTest* parent) : EAbstractAnalyticSerial(parent), _base(parent)
{
    EDEBUG_FUNC(this,parent);
}



/*!
 * Perform the desired work on a result block.
 *
 * @param block The work block you are working on.
 *
 * @return The populated result block.
 */
std::unique_ptr<EAbstractAnalyticBlock> ConditionalTest::Serial::execute(const EAbstractAnalyticBlock* block)
{
    EDEBUG_FUNC(this, block);

    // crate the work and result blocks
    const WorkBlock* workBlock {block->cast<WorkBlock>()};
    ResultBlock* resultBlock {new ResultBlock(workBlock->index())};

    // Create iterators for the CCM data object.
    CCMatrix::Pair ccmPair = CCMatrix::Pair(_base->_ccm);
    CorrelationMatrix::Pair cmxPair = CorrelationMatrix::Pair(_base->_cmx);

    // Iterate through the pairs the workblock is assigned to.
    Pairwise::Index index(workBlock->start());

    // iterate through each pair in the matrix
    for ( qint64 rawIndex = 0; rawIndex < workBlock->size(); ++rawIndex, ++index )
    {
        // reads the first value in the ccm
        if ( rawIndex == 0 )
        {
            ccmPair.read(index);
            cmxPair.read(ccmPair.index());
        }
        else
        {
            ccmPair.readNext();
            cmxPair.readNext();
        }

        // print warning if pairwise indices do not match
        if ( ccmPair.index() != cmxPair.index() )
        {
            qInfo() << "warning: ccm and cmx files are out of sync";
        }

        // if the first one isnt in the cluster we should not count it.
        if ( ccmPair.clusterSize() == 0 )
        {
            continue;
        }

        // Initialize new pvalues, one set of pvalues for each cluster.
        QVector<QVector<double>> pValues(ccmPair.clusterSize());

        // Initialize new r2, one set of pvalues for each cluster.
        QVector<QVector<double>> r2(ccmPair.clusterSize());

        // for each cluster in the pair, run the binomial and linear regression tests.
        for ( qint32 clusterIndex = 0; clusterIndex < ccmPair.clusterSize(); clusterIndex++ )
        {
            // resize for room for each test.
            pValues[clusterIndex].resize(_base->_numTests);
            r2[clusterIndex].resize(_base->_numTests);

            for ( qint32 featureIndex = 0, testIndex = 0; featureIndex < _base->_features.size(); featureIndex++ )
            {
                if ( _base->_testType.at(featureIndex) == NONE || _base->_testType.at(featureIndex) == UNKNOWN )
                {
                    continue;
                }

                if ( _base->_testType.at(featureIndex) == QUANTITATIVE || _base->_testType.at(featureIndex) == ORDINAL )
                {
                    prepAnxData(_base->_features.at(featureIndex).at(0), featureIndex, _base->_testType.at(featureIndex));
                    test(ccmPair, clusterIndex, testIndex, featureIndex, 0, pValues, r2);
                }
                else if ( _base->_testType.at(featureIndex) == CATEGORICAL )
                {
                    for ( qint32 labelIndex = 0; labelIndex < _base->_features.at(featureIndex).size(); labelIndex++ )
                    {
                        prepAnxData(_base->_features.at(featureIndex).at(labelIndex), featureIndex, _base->_testType.at(featureIndex));

                        // if there are sub labels to test for the feature
                        if ( _base->_features.at(featureIndex).size() > 1 )
                        {
                            if ( labelIndex == 0 )
                            {
                                labelIndex = 1;
                            }
                            test(ccmPair, clusterIndex, testIndex, featureIndex, labelIndex, pValues, r2);
                        }

                        // if only the feature needs testing (no sub labels)
                        else
                        {
                            test(ccmPair, clusterIndex, testIndex, featureIndex, 0, pValues, r2);
                        }
                    }
                }
            }
        }

        if ( !isEmpty(pValues) )
        {
            resultBlock->append(Pair {
                index,
                pValues,
                r2
            });
        }
    }

    return std::unique_ptr<EAbstractAnalyticBlock>(resultBlock);
}



/*!
 * Prepare the annotation matrix data for testing.
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

    // get the needed data fro the comparison
    _catCount = 0;
    _amxData.resize(_base->_data.at(dataIndex).size());

    // populate array with annotation data
    for ( int j = 0; j < _base->_data.at(dataIndex).size(); j++ )
    {
        _amxData[j] = _base->_data.at(dataIndex).at(j).toString();

        // if data is the same as the test label add one to the catagory counter
        if ( testType == _base->CATEGORICAL && _amxData[j] == testLabel )
        {
            _catCount++;
        }
    }

    return _catCount;
}



/*!
 * Check to see if a matrix is empty.
 *
 * @param vector The matrix you want to check.
 *
 * @return True if the matrix is empty, false otherwise.
 */
bool ConditionalTest::Serial::isEmpty(QVector<QVector<double>>& matrix)
{
    EDEBUG_FUNC(this, &matrix);

    int index = 0;

    while ( index < matrix.size() )
    {
        if ( !matrix.at(index++).isEmpty() )
        {
            return false;
        }
    }

    return true;
}



/*!
 * Prepare the cluster category count information.
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
        if ( testType == _base->CATEGORICAL && _amxData.at(i) == label )
        {
            _catCount++;
        }

        if ( ccmPair.at(clusterIndex, i) == 1 )
        {
            _clusterSize++;

            if ( testType == _base->CATEGORICAL && _amxData.at(i) == label )
            {
                _catInCluster++;
            }
        }
    }

    return _catInCluster;
}



/*!
 * An interface to choose and run the correct tests on the pair.
 *
 * @param ccmPair The pair thats going to be tested.
 *
 * @param clusterIndex The cluster to test inside the pair, each cluster is
 *       tested.
 *
 * @param testIndex Which test we are currently performing.
 *
 * @param featureIndex The current feature we are testing.
 *
 * @param labelIndex The label in the feature we are running a test on.
 *
 * @param pValues The two dimensional array holding all of the results from the
 *       tests.
 *
 * @return The test that was just conducted.
 */
int ConditionalTest::Serial::test(
    CCMatrix::Pair& ccmPair,
    qint32 clusterIndex,
    qint32& testIndex,
    qint32 featureIndex,
    qint32 labelIndex,
    QVector<QVector<double>>& pValues,
    QVector<QVector<double>>& r2)
{
    EDEBUG_FUNC(this,&ccmPair, clusterIndex, &testIndex, featureIndex, labelIndex, &pValues);

    // get informatiopn on the mask
    clusterInfo(ccmPair, clusterIndex, _base->_features.at(featureIndex).at(labelIndex), _base->_testType.at(featureIndex));

    // For linear regresssion we need a variable that will hold the
    // pvalue and the r2 value.
    QVector<double> results(2);

    // conduct the correct test based on the type of data
    switch(_base->_testType.at(featureIndex))
    {
        case CATEGORICAL:
            pValues[clusterIndex][testIndex] = hypergeom(ccmPair, clusterIndex,
                                                         _base->_features.at(featureIndex).at(labelIndex));
            r2[clusterIndex][testIndex] = qQNaN();
            testIndex++;
            break;
        case ORDINAL:
            regression(_amxData, ccmPair, clusterIndex, ORDINAL, results);
            pValues[clusterIndex][testIndex] = results.at(0);
            r2[clusterIndex][testIndex] = results.at(1);
            testIndex++;
            break;
        case QUANTITATIVE:
            regression(_amxData, ccmPair, clusterIndex, QUANTITATIVE, results);
            pValues[clusterIndex][testIndex] = results.at(0);
            r2[clusterIndex][testIndex] = results.at(1);
            testIndex++;
            break;
        default:
            break;
    }

    return _base->_testType.at(featureIndex);
}



/*!
 * Run the first binomial test for given data.
 *
 * @return Pvalue corrosponding to the test.
 */
double ConditionalTest::Serial::hypergeom(CCMatrix::Pair& ccmPair, int clusterIndex, QString test_label)
{
    EDEBUG_FUNC(this);

    // We use the hypergeometric distribution because the samples are
    // selected from the population for membership in the cluster without
    // replacement.

    // If a population contains n_1 elements of “type 1” and n_2 elements of
    // “type 2” then the hypergeometric distribution gives the probability
    // of obtaining k elements of “type 1” in t samples from the population.

    int sampleSize =  _base->_emx->sampleSize();

    // Population contains n1 elements of Type 1.
    int n1 = _catCount;
    // Population contains n2 elements of Type 2.
    int n2 = sampleSize - _catCount;
    // k elements of Type 1 were selected.
    int k = _catInCluster;
    // t total elements were selected.
    int t = _clusterSize;
    // Holds the pvalue
    double pvalue = 1;

    // If n1 == k we will always get a zero because we've
    // reached the end of the distribution, so the Ho that
    // X > k is always 0.  This happens if the cluster is 100%
    // comprised of the category we're looking for.
    if ( k == n1 )
    {
        return 1;
    }

    // If our dataset is large, the power to detect the effect
    // size increases, resulting in potentially insignificant
    // proportions having signficant p-values. Using Cohen's H
    // to set a large-effect size (e.g. 0.8) with a sig.level of
    // 0.001 and a power of 0.95 we need at least 31 samples.
    // So, we'll perform a bootstrap resampling of our data
    // to calculate an average proportion of 31 samples
    int bs_t = 31;
    if ( t > bs_t )
    {
        // Initialize the uniform random number generator.
        const gsl_rng_type * T;
        gsl_rng * r;
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);

        // Holds the bootstrap average proportion.
        int jkap = 0;

        // To perform the bootstrap resampling we will
        // perform 30 iterations (central limit thereom)
        int in = 30;
        for ( int i = 0; i < in; i++ )
        {
            // Keeps track of the number of successes for each iteration.
            int ns = 0;

            // Generate 31 random indexes between 0 and the size of
            // the sample string.  We will use these numbers to
            // randomly select a sample and if it is a 1 and of the
            // testing category then we consider it a success.
            int indexes[t];
            int chosen[bs_t];
            int jx = 0;
            for ( int j = 0; j < sampleSize; j++ )
            {
                if ( ccmPair.at(clusterIndex, j) == 1 )
                {
                    indexes[jx++] = j;
                }
            }

            // The gsl_ran_sample function randomly chooses samples with replacement from a list.
            gsl_ran_sample(r, chosen, bs_t, indexes, t, sizeof(int));

            // Now count the number of randomly selected items from the cluster.
            for ( int j = 0; j < bs_t; j++ )
            {
                if ( ccmPair.at(clusterIndex, chosen[j]) == 1 && _amxData.at(chosen[j]) == test_label )
                {
                    ns = ns + 1;
                }
            }

            jkap += ns;
        }

        // Calculate the average proportion from all iterations
        // and free the random number struct.
        jkap = jkap/in;
        gsl_rng_free(r);

        // Since we want to test P(X >= jkap) not P(X > jkap) we should
        // subtract 1 from jkap but we shouldn't go below zero.
        jkap = jkap - 1;
        if (jkap < 0)
        {
            jkap = 0;
        }

        // Now reset the sample size and the proporiton of success.
        pvalue = gsl_cdf_hypergeometric_Q(jkap, n1, n2, bs_t);
        return pvalue;
    }

    // Since we want to test P(X >= k) not P(X > k) we should
    // subtract 1 from k but we shouldn't go below zero.
    k = k - 1;
    if (k < 0)
    {
        k = 0;
    }

    // The gsl_cdf_hypergeometric_Q function uses the upper-tail of the CDF.
    pvalue = gsl_cdf_hypergeometric_Q(k, n1, n2, t);
    return pvalue;
}



/*!
 * Run the regression test for given data, the regression line is genex vs geney vs label data.
 *
 * @param amxInfo Annotation matrix information.
 *
 * @param ccmPair Cluster matrix pair.
 *
 * @param clusterIndex The index of the cluster, used to get the right info
 *       from the cluster pair
 *
 * @return Pvalue corrosponding to the test.
 */
void ConditionalTest::Serial::regression(QVector<QString> &amxInfo, CCMatrix::Pair& ccmPair, int clusterIndex, TESTTYPE testType, QVector<double>& results)
{
    EDEBUG_FUNC(this, &amxInfo, &ccmPair, clusterIndex);

    // Temp containers.
    QVector<double> labelInfo;

    // Regression model containers.
    double chisq;
    gsl_matrix *X, *cov;
    gsl_vector *Y, *C;

    // Before allocating memory for the linear regression let's remove any
    // values in the expression matrix that have missing values. This will
    // tell us the full size we need to reserve.
    int test_cluster_size = 0;
    for ( int i = 0; i < _base->_emx->sampleSize(); i++ )
    {
        // If the sample label matches with the given label.
        if ( ccmPair.at(clusterIndex, i) == 1 )
        {
            if (QString::compare(amxInfo.at(i), _base->_missing) == 0) {
                continue;
            }
            test_cluster_size++;
        }
    }

    // Allocate a matrix to hold the predictior variables, in this case the gene
    // Expression data.
    X = gsl_matrix_alloc(test_cluster_size, 4);

    // Allocate a vector to hold observation data, in this case the data
    // corrosponding to the features.
    Y = gsl_vector_alloc(test_cluster_size);

    // Allocate a vector and matrix for the slope info.
    C = gsl_vector_alloc(4);
    cov = gsl_matrix_alloc(4, 4);

    // Read in the gene pairs expression information.
    ExpressionMatrix::Gene geneX(_base->_emx);
    ExpressionMatrix::Gene geneY(_base->_emx);
    geneX.read(ccmPair.index().getX());
    geneY.read(ccmPair.index().getY());

    // QString g1Name = geneX.toString();
    // QString g2Name = geneY.toString();

    // Look through all the samples in the mask.
    for ( int i = 0, j = 0; i < _base->_emx->sampleSize(); i++ )
    {
        // If the sample label matches with the given label.
        if ( ccmPair.at(clusterIndex, i) == 1 )
        {
            // Skip samples with a missing value in the annotation matrix
            if (QString::compare(amxInfo.at(i), _base->_missing) == 0) {
                continue;
            }

            // Add emx data as the predictors but only for samples in the cluster. We
            // add a 1 for the intercept, gene1, gene2 and the interaction: gene1*gene2.
            // The interaction term handles the case where the relationship is dependent
            // on the values of both genes.
            double g1 = static_cast<double>(geneX.at(i));
            double g2 = static_cast<double>(geneY.at(i));

            gsl_matrix_set(X, j, 0, static_cast<double>(1)); // for the intercept
            gsl_matrix_set(X, j, 1, g1);
            gsl_matrix_set(X, j, 2, g2);
            gsl_matrix_set(X, j, 3, g1*g2);
            gsl_vector_set(Y, j, amxInfo.at(i).toFloat());

            j++;
        }
    }

    // Create the workspace for the gnu scientific library to work in.
    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (X->size1, X->size2);

    // Regression calculation.
    gsl_multifit_linear(X, Y, C, cov, &chisq, work);

    // Calculate some quantities we'll need for R^2 and p-value
    // calculations.
    // Sum of Squares Total
    double SST = gsl_stats_tss(Y->data, Y->stride, Y->size);
    // Sum of Squares of the Error (residuals)
    double SSE = chisq;

    // Sum of Squares for the model.
    double SSM = SST - SSE;

    // Degrees of Freedom for the Model. This is p-1 where
    // p is the number of parameters in the model.
    double DFM = 3;
    // Degrees of Freedom for the Error. This is n-p where
    // n is the number of samples.
    double DFE = test_cluster_size - 4;

    // Mean of Squares for the Model.
    double MSM = SSM/DFM;
    // Mean of Squares for Error.
    double MSE = SSE/DFE;

    // Calculate R^2 and R^2-adjusted
    double R2 = 1 - (SSE / SST);
    double R2adj = 1.0 - ((double) test_cluster_size - 1) / DFE * (1 - R2);

    // Calculate the p-value. We will do this using the F-test.
    double Fstat = MSM/MSE;
    double pValue = 1 - gsl_cdf_fdist_P(Fstat, DFM, DFE);


    // TODO: we should check the assumptions of the linear regression and
    // not return if the assumptions are not met.

    // TODO: it would be nice to return a rate of change of the conditioal mean.

    // TODO: it would be nice to return the r2 value along with the p-value.

    // Four scenarios:
    // 1) low R-square and low p-value (p-value <= 0.05).  Model doesn't explain
    //    the variation but it does follow the trend or regression line well.
    // 2) low R-square and high p-value (p-value > 0.05).  Model doesn't explai
    //    the variation and it doesn't follow the trend line.
    // 3) high R-square and low p-value.  Model explains the variance it
    //    follows the trend line.
    // 4) high R-square and high p-value. Model explains the variance well but
    //    it does not follow the trend line very well.
    // In summary, low p-values still indicate a real relationship between the
    // predictors and the observed values.

    // Free all of the data.
    gsl_matrix_free(X);
    gsl_vector_free(Y);
    gsl_matrix_free(cov);
    gsl_vector_free(C);
    gsl_multifit_linear_free(work);

    // Set the results array
    if ( qIsNaN(pValue) )
    {
        results[0] = 1;
        results[1] = 0;
    }
    else {
        results[0] = pValue;
        results[1] = R2;
    }
}
