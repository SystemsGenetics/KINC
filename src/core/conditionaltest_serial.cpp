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
#include <boost/math/distributions/non_central_f.hpp>
#include <math.h>
#include <iostream>


/*!
 * Create a serial object.
 */
ConditionalTest::Serial::Serial(ConditionalTest* parent) :
    EAbstractAnalyticSerial(parent),
    _base(parent)
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

    // create iterators for the CCM and CMX data objects
    CCMatrix::Pair ccmPair(_base->_ccm);
    CorrelationMatrix::Pair cmxPair(_base->_cmx);

    // read the first pair of the work block
    Pairwise::Index index(workBlock->start());

    ccmPair.read(index);
    cmxPair.read(index);

    // iterate through each pair in the work block
    for ( qint64 wbIndex = 0; wbIndex < workBlock->size(); ++wbIndex )
    {

        // print warning if iterator indices do not match
        if ( ccmPair.index() != cmxPair.index() )
        {
            QString source = _base->_cmx->geneNames().at(cmxPair.index().getX()).toString();
            QString target = _base->_cmx->geneNames().at(cmxPair.index().getY()).toString();
            qInfo() << "warning: ccm and cmx files are out of sync at cmx coordinate:"
                    << source << "," << target << " ("
                    << cmxPair.index().getX() << "," << cmxPair.index().getY() <<").";
        }

        // initialize set of p values and r2 values for each cluster
        QVector<QVector<double>> pValues(ccmPair.clusterSize());
        QVector<QVector<double>> r2(ccmPair.clusterSize());

        // for each cluster in the pair, run the binomial or linear regression tests
        for ( qint32 clusterIndex = 0; clusterIndex < ccmPair.clusterSize(); clusterIndex++ )
        {

            // resize for room for each test.
            pValues[clusterIndex].resize(_base->_numTests);
            r2[clusterIndex].resize(_base->_numTests);


            // Iterate through all of the tests.
            int test_index = 0;
            for ( qint32 featureIndex = 0; featureIndex < _base->_features.size(); featureIndex++ )
            {
                // Get the column data from the annotation matrix for this feature.
                int num_samples = _base->_data.at(featureIndex).size();
                QVector<QString> amx_column(num_samples);
                for ( int j = 0; j < num_samples; j++ )
                {
                    amx_column[j] = _base->_data.at(featureIndex).at(j).toString();
                }

                if ( _base->_testType.at(featureIndex) == QUANTITATIVE ||
                     _base->_testType.at(featureIndex) == ORDINAL )
                {
                    // For linear regresssion we need a variable that will hold the
                    // pvalue and the r2 value.
                    QVector<double> results(2);
                    regression(amx_column, ccmPair, clusterIndex, results);
                    pValues[clusterIndex][test_index] = results.at(0);
                    r2[clusterIndex][test_index] = results.at(1);                    
                    test_index++;
                }
                else if ( _base->_testType.at(featureIndex) == CATEGORICAL )
                {
                    // Loop through each label (category) of the feature.
                    for ( qint32 labelIndex = 1; labelIndex < _base->_features.at(featureIndex).size(); labelIndex++ )
                    {
                        double pval;
                        //hypergeom(amx_column, ccmPair, clusterIndex, featureIndex, labelIndex, pval);
                        test_proportions(amx_column, ccmPair, clusterIndex, featureIndex, labelIndex, pval);
                        pValues[clusterIndex][test_index] = pval;
                        r2[clusterIndex][test_index] = qQNaN();
                        test_index++;
                    }
                }
            }
        }

        // append pair to result block
        resultBlock->append(Pair {
            cmxPair.index(),
            pValues,
            r2
        });

        // read the next pair
        cmxPair.readNext();
        ccmPair.read(cmxPair.index());
    }

    return std::unique_ptr<EAbstractAnalyticBlock>(resultBlock);
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
 * The hypergeometric test
 *
 */
void ConditionalTest::Serial::hypergeom(
    const QVector<QString>& amx_column,
    const CCMatrix::Pair& ccmPair,
    int clusterIndex,
    int featureIndex,
    int labelIndex,
    double &result)
{
    EDEBUG_FUNC(this);

    // Get the test name and the test type.
    QString test_label = _base->_features.at(featureIndex).at(labelIndex);
    TestType test_type = _base->_testType.at(featureIndex);

    // We need to get the data for the feature being tested and
    // count how many entries there are for the label in the
    // annotation data.
    int num_samples = _base->_data.at(featureIndex).size();
    int cluster_size = 0;
    int total_label_count = 0;
    int labels_in_cluster = 0;
    for ( int j = 0; j < num_samples; j++ )
    {

        // if data is the same as the test label add one to the catagory counter
        if ( test_type == _base->CATEGORICAL && amx_column[j] == test_label )
        {
            total_label_count++;
        }
        if ( ccmPair.at(clusterIndex, j) == 1 )
        {
            cluster_size++;
            if ( test_type == _base->CATEGORICAL && amx_column.at(j) == test_label )
            {
                labels_in_cluster++;
            }
        }
    }

    // If there are no matching labels in this cluster then return. This
    // could happen if the annotation matrix has all NAs for the cluster.
    if ( labels_in_cluster == 0 )
    {
        result = qQNaN();
        return;
    }


    // We use the hypergeometric distribution because the samples are
    // selected from the population for membership in the cluster without
    // replacement.

    // If a population contains n_1 elements of “type 1” and n_2 elements of
    // “type 2” then the hypergeometric distribution gives the probability
    // of obtaining k elements of “type 1” in t samples from the population.

    int sampleSize =  _base->_emx->sampleSize();

    // Population contains n1 elements of Type 1.
    int n1 = total_label_count;
    // Population contains n2 elements of Type 2.
    int n2 = sampleSize - total_label_count;
    // k elements of Type 1 were selected.
    int k = labels_in_cluster;
    // t total elements were selected.
    int t = cluster_size;

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
                if ( ccmPair.at(clusterIndex, chosen[j]) == 1 && amx_column.at(chosen[j]) == test_label )
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

        // The gsl_cdf_hypergeometric_Q function uses the upper-tail of the CDF.
        result = gsl_cdf_hypergeometric_Q(jkap, n1, n2, bs_t);
        return;
    }

    // Since we want to test P(X >= k) not P(X > k) we should
    // subtract 1 from k but we shouldn't go below zero.
    k = k - 1;
    if (k < 0)
    {
        k = 0;
    }

    // The gsl_cdf_hypergeometric_Q function uses the upper-tail of the CDF.
    result = gsl_cdf_hypergeometric_Q(k, n1, n2, t);
    return;
}

/*!
 * The z test of proprotions
 *
 */
void ConditionalTest::Serial::test_proportions (
    const QVector<QString>& amx_column,
    const CCMatrix::Pair& ccmPair,
    int clusterIndex,
    int featureIndex,
    int labelIndex,
    double &result)
{
    EDEBUG_FUNC(this);

    // Get the test name and the test type.
    QString test_label = _base->_features.at(featureIndex).at(labelIndex);
    TestType test_type = _base->_testType.at(featureIndex);

    /* For debugging purposes
    QString source_name = _base->_cmx->geneNames().at(ccmPair.index().getX()).toString();
    QString target_name = _base->_cmx->geneNames().at(ccmPair.index().getY()).toString();

    if (source_name == "pycom14g13710" && target_name == "pycom02g20220") {
      std::cout << "Found it";
    }*/

    // We need to get the data for the feature being tested and
    // count how many entries there are for the label in the
    // annotation data.
    int num_samples = _base->_data.at(featureIndex).size();
    int cluster_size = 0;
    int total_label_count = 0;
    int total_non_missing = 0;
    int labels_in_cluster = 0;
    int labels_out_cluster = 0;
    int jx = 0;
    for ( int j = 0; j < num_samples; j++ )
    {

        // if data is the same as the test label add one to the catagory counter
        if ( test_type == _base->CATEGORICAL && amx_column[j] == test_label )
        {
            total_label_count++;
        }
        if ( ccmPair.at(clusterIndex, j) == 1 )
        {
            cluster_size++;
            if ( test_type == _base->CATEGORICAL && amx_column.at(j) == test_label )
            {
                labels_in_cluster++;
            }
        }
        if ( ccmPair.at(clusterIndex, j) != 1 && ccmPair.at(clusterIndex, j) != 9)
        {
            if ( test_type == _base->CATEGORICAL && amx_column.at(j) == test_label )
            {
                labels_out_cluster++;
            }
        }
        if ( ccmPair.at(clusterIndex, j) != 9 )
        {
          total_non_missing++;
        }
    }

    // If there are no matching labels in this cluster then return. This
    // could happen if the annotation matrix has all NAs for the cluster.
    if ( labels_in_cluster == 0 )
    {
        result = qQNaN();
        return;
    }

    int sampleSize =  _base->_emx->sampleSize();


    // Before bootstraping get the list of indexes that belong to this cluster.
    // We will randomly select from these indexes.
    int cluster_indexes[cluster_size];
    for ( int j = 0; j < sampleSize; j++ )
    {
        if ( ccmPair.at(clusterIndex, j) == 1 )
        {
            cluster_indexes[jx++] = j;
        }
    }

    // Initialize the uniform random number generator.
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // Perform 30 iterations of the boostrap sampling.
    int num_i = 30;
    // Select 30 items for each iteration
    int bs_t = 30;
    // Holds the bootstrap average proportion.
    double avg_prop = 0.0;
    for ( int i = 0; i < num_i; i++ )
    {
        // Keeps track of the number of successes for each iteration.
        int ns = 0;

        // The gsl_ran_sample function randomly chooses samples with replacement from a list.
        int chosen[bs_t];
        gsl_ran_sample(r, chosen, bs_t, cluster_indexes, cluster_size, sizeof(int));

        // Now count the number of randomly selected items from the cluster.
        for ( int j = 0; j < bs_t; j++ )
        {
            QString slabel = amx_column.at(chosen[j]);
            if ( ccmPair.at(clusterIndex, chosen[j]) == 1 && amx_column.at(chosen[j]) == test_label )
            {
                ns = ns + 1;
            }
        }

        avg_prop += static_cast<double>(ns) / static_cast<double>(bs_t);
    }

    // Calculate the average proportion from all iterations
    // and free the random number struct.
    avg_prop = avg_prop / num_i;
    gsl_rng_free(r);

    // First test.
    // The null hypothesis for this test is that the cluster is
    // comprised of equal or less of the desired label. This would
    // happen if all labels have equal chance, or there is a bias
    // against the label.
    // Ho: avg_prop <= 0.5. The value 0.5 is chosen with the asumption
    //     that there is equal or less proportions of the desired label.
    // Ha: avg_prop > 0.5
    double exp_prop = 0.5;
    double z = (avg_prop - exp_prop) / sqrt((exp_prop * (1 - exp_prop)) / bs_t);
    double pval1 = 1 - gsl_cdf_gaussian_P(z, 1);

    // Second, test.
    // The null hypothesis is that the proportion of labels in the cluster
    // is equal or less than the proportion out of the cluster.  This would
    // happen if all the labels have equal chance of being in the cluster or
    // if there is bias against the label.
    double obs_prop = (avg_prop * cluster_size) / static_cast<double>(total_label_count);
    obs_prop = std::min(1.0, obs_prop);
    exp_prop = static_cast<double>(total_label_count) / static_cast<double>(num_samples);
    z = (obs_prop - exp_prop) / sqrt((exp_prop * (1 - exp_prop)) / bs_t);
    double pval2 = 1 - gsl_cdf_gaussian_P(z, 1);

    // Return the maximum p-value
    result = std::max(pval1, pval2);
}

/*!
 * Performs the regression test for quantitative data.
 */
void ConditionalTest::Serial::regression(
    const QVector<QString>& amx_column,
    const CCMatrix::Pair& ccmPair,
    int clusterIndex,
    QVector<double>& results)
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
            if (QString::compare(amx_column.at(i), _base->_missing) == 0) {
                continue;
            }
            test_cluster_size++;
        }
    }

    // Don't do regression analysis if there are fewer than 5 samples.
    // This can occur if the annotation matrix has missing values
    // for samples in the cluster. The Degrees of Freedom for the error
    // is test_cluster_size - 4 so we need at least 5 samples. Odds
    // are these tests will get filterd out anyway for low power.
    if ( test_cluster_size < 5 )
    {
        results[0] = qQNaN();
        results[1] = qQNaN();
        return;
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

    // Useful for debugging purposes.
    // QString g1Name = geneX.toString();
    // QString g2Name = geneY.toString();

    // Look through all the samples in the mask.
    for ( int i = 0, j = 0; i < _base->_emx->sampleSize(); i++ )
    {
        // If the sample label matches with the given label.
        if ( ccmPair.at(clusterIndex, i) == 1 )
        {
            // Skip samples with a missing value in the annotation matrix
            if (QString::compare(amx_column.at(i), _base->_missing) == 0) {
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
            gsl_vector_set(Y, j, amx_column.at(i).toFloat());

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
    //double R2adj = 1.0 - ((double) test_cluster_size - 1) / DFE * (1 - R2);

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
        results[0] = qQNaN();
        results[1] = qQNaN();
    }
    else {               
        double power = pwr_f2_test(3, test_cluster_size - 4, R2, _base->_powerThresholdAlpha);
        if ( power >= _base->_powerThresholdPower ) {
            results[0] = pValue;
            results[1] = R2;
        }
        else {
            results[0] = std::numeric_limits<double>::quiet_NaN();
            results[1] = std::numeric_limits<double>::quiet_NaN();
        }
    }
}

/**
 * @brief ConditionalTest::Serial::pwr_f2_test
 * @param df1 degrees of freedom of the numerator: the number of coefficients in
 *        the model minus the intercept.
 * @param df2 degrees of freedom of the denominator: the number of error degrees
 *        of freedom which is n - df1 - 1
 * @param f2 effect size
 * @param sig_level signficance level
 * @return
 */
double ConditionalTest::Serial::pwr_f2_test(int df1, int df2, double f2, double sig_level)
{
    double lambda = f2 * (df1 + df2 + 1);
    double p = gsl_cdf_fdist_Qinv(sig_level, df1, df2);
    boost::math::non_central_f_distribution<> ncfd(df1, df2, lambda);
    double power = 1 - boost::math::cdf(ncfd, p);
    return power;
}
