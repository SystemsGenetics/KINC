#include <cblas.h>
#include <cusolverDn.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <lapacke.h>

#define MAJOR_VERSION KINC_MAJOR_VERSION
#define MINOR_VERSION KINC_MINOR_VERSION
#define REVISION KINC_REVISION
#include <ace/core/ace_settings.h>
#include <ace/core/cudaxx.h>

#include "rmt.h"
#include "rmt_input.h"



using namespace std;
using RawPair = CorrelationMatrix::RawPair;



#define CHECK_ERROR(condition, message)  \
    if ( !(condition) )                  \
    {                                    \
        E_MAKE_EXCEPTION(e);             \
        e.setTitle(tr("CUDA Error"));    \
        e.setDetails(tr(message));       \
        throw e;                         \
    }



#define CHECK_CUDA(ret) \
    CHECK_ERROR(ret == cudaSuccess, #ret)



#define CHECK_CUSOLVER(ret) \
    CHECK_ERROR(ret == CUSOLVER_STATUS_SUCCESS, #ret)



/*!
 * Return the total number of blocks this analytic must process as steps
 * or blocks of work.
 */
int RMT::size() const
{
    EDEBUG_FUNC(this);

    return 1;
}



/*!
 * Process the given index with a possible block of results if this analytic
 * produces work blocks. This analytic implementation has no work blocks.
 *
 * @param result
 */
void RMT::process(const EAbstractAnalyticBlock*)
{
    EDEBUG_FUNC(this);

    // initialize log text stream
    QTextStream stream(_logfile);

    // initialize helper variables
    float finalThreshold {0};
    float finalChi {numeric_limits<float>::infinity()};
    float maxChi {-numeric_limits<float>::infinity()};

    float threshold {_thresholdStart};

    // load raw correlation data, row-wise maximums
    std::vector<RawPair> pairs {_input->dumpRawData()};
    std::vector<float> maximums {computeMaximums(pairs)};

    // continue while max chi is less than final threshold
    while ( maxChi < _chiSquareThreshold2 )
    {
        qInfo("\n");
        qInfo("threshold: %0.3f", threshold);

        // compute pruned matrix based on threshold
        size_t size;
        std::vector<float> pruneMatrix {computePruneMatrix(pairs, maximums, threshold, &size)};

        qInfo("prune matrix: %lu", size);

        // make sure that pruned matrix is not empty
        float chi = -1;
        std::vector<float> eigens;

        if ( size > 0 )
        {
            // compute eigenvalues of pruned matrix
            eigens = computeEigenvalues(&pruneMatrix, size);

            qInfo("eigenvalues: %lu", eigens.size());

            // compute unique eigenvalues
            eigens = computeUnique(eigens);

            qInfo("unique eigenvalues: %lu", eigens.size());

            // compute chi-squared value from NNSD of eigenvalues
            chi = computeChiSquare(eigens);

            qInfo("chi-squared: %g", chi);
        }

        // make sure that chi-squared test succeeded
        if ( chi != -1 )
        {
            // save the most recent chi-squared value less than critical value
            if ( chi < _chiSquareThreshold1 )
            {
                finalChi = chi;
                finalThreshold = threshold;
            }

            // save the largest chi-squared value which occurs after finalChi
            if ( finalChi < _chiSquareThreshold1 && chi > finalChi )
            {
                maxChi = chi;
            }
        }

        // output to log file
        stream
            << QString::number(threshold, 'f', 3) << "\t"
            << size << "\t"
            << eigens.size() << "\t"
            << chi << "\n";

        // decrement threshold and fail if minimum threshold is reached
        threshold -= _thresholdStep;
        if ( threshold < _thresholdStop )
        {
            E_MAKE_EXCEPTION(e);
            e.setTitle(tr("RMT Threshold Error"));
            e.setDetails(tr("Could not find non-random threshold above stopping threshold."));
            throw e;
        }
    }

    // write threshold where chi was first above final threshold
    stream << finalThreshold << "\n";
}



/*!
 * Make a new input object and return its pointer.
 */
EAbstractAnalyticInput* RMT::makeInput()
{
    EDEBUG_FUNC(this);

    return new Input(this);
}



/*!
 * Initialize this analytic. This implementation checks to make sure the input
 * data object and output log file have been set, and that various integer
 * arguments are valid.
 */
void RMT::initialize()
{
    EDEBUG_FUNC(this);

    // make sure input and output were set properly
    if ( !_input || !_logfile )
    {
        E_MAKE_EXCEPTION(e);
        e.setTitle(tr("Invalid Argument"));
        e.setDetails(tr("Did not get valid input or logfile arguments."));
        throw e;
    }

    // make sure threshold arguments are valid
    if ( _thresholdStart <= _thresholdStop )
    {
        E_MAKE_EXCEPTION(e);
        e.setTitle(tr("Invalid Argument"));
        e.setDetails(tr("Starting threshold must be greater than stopping threshold."));
        throw e;
    }

    // make sure pace arguments are valid
    if ( _minSplinePace >= _maxSplinePace )
    {
        E_MAKE_EXCEPTION(e);
        e.setTitle(tr("Invalid Argument"));
        e.setDetails(tr("Minimum spline pace must be less than maximum spline pace."));
        throw e;
    }

    // set the number of openblas threads
    openblas_set_num_threads(_numThreads);
}



/*!
 * Compute the row-wise maximums of a correlation matrix.
 *
 * @param pairs
 */
std::vector<float> RMT::computeMaximums(const std::vector<RawPair>& pairs)
{
    EDEBUG_FUNC(this,&pairs);

    // initialize elements to minimum value
    std::vector<float> maximums(_input->geneSize(), 0);

    // compute maximum correlation of each row
    for ( auto& pair : pairs )
    {
        int i = pair.index.getX();

        for ( size_t k = 0; k < pair.correlations.size(); ++k )
        {
            float correlation = fabs(pair.correlations[k]);

            if ( maximums[i] < correlation )
            {
                maximums[i] = correlation;
            }
        }
    }

    // return row-wise maximums
    return maximums;
}



/*!
 * Compute the pruned matrix of a correlation matrix with a given threshold. This
 * function uses the pre-computed row-wise maximums for faster computation. The
 * returned matrix is equivalent to the correlation matrix with all correlations
 * below the given threshold removed, and all zero-columns removed. Additionally,
 * the number of rows in the pruned matrix is returned as a pointer argument.
 *
 * @param pairs
 * @param maximums
 * @param threshold
 * @param size
 */
std::vector<float> RMT::computePruneMatrix(const std::vector<RawPair>& pairs, const std::vector<float>& maximums, float threshold, size_t* size)
{
    EDEBUG_FUNC(this,&pairs,&maximums,threshold,size);

    // generate vector of row indices that have a correlation above threshold
    std::vector<int> indices(_input->geneSize(), -1);
    size_t pruneSize = 0;

    for ( size_t i = 0; i < maximums.size(); ++i )
    {
        if ( maximums[i] >= threshold )
        {
            indices[i] = pruneSize;
            pruneSize++;
        }
    }

    // extract pruned matrix from correlation matrix
    std::vector<float> pruneMatrix(pruneSize * pruneSize);

    // initialize diagonal
    for ( size_t i = 0; i < pruneSize; ++i )
    {
        pruneMatrix[i * pruneSize + i] = 1;
    }

    // iterate through all pairs
    for ( auto& pair : pairs )
    {
        // get indices into pruned matrix
        int i = indices[pair.index.getX()];
        int j = indices[pair.index.getY()];

        // skip pair if it was pruned
        if ( i == -1 || j == -1 )
        {
            continue;
        }

        // select correlation from pair using reduction method
        float correlation = 0;

        switch ( _reductionMethod )
        {
        case ReductionMethod::First:
        {
            correlation = pair.correlations[0];
            break;
        }
        case ReductionMethod::MaximumCorrelation:
        {
            for ( size_t k = 0; k < pair.correlations.size(); k++ )
            {
                float r = fabs(pair.correlations[k]);

                if ( correlation < r )
                {
                    correlation = r;
                }
            }
            break;
        }
        case ReductionMethod::MaximumSize:
        {
            E_MAKE_EXCEPTION(e);
            e.setTitle(tr("Unsupported Option"));
            e.setDetails(tr("Pairwise reduction by maximum size is not yet supported."));
            throw e;
        }
        case ReductionMethod::Random:
        {
            int k = qrand() % pair.correlations.size();
            correlation = pair.correlations[k];
            break;
        }
        };

        // save correlation if it is above threshold
        if ( fabs(correlation) >= threshold )
        {
            pruneMatrix[i * pruneSize + j] = correlation;
            pruneMatrix[j * pruneSize + i] = correlation;
        }
    }

    // save size of pruned matrix
    *size = pruneSize;

    // return pruned matrix
    return pruneMatrix;
}



/*!
 * Compute the eigenvalues of a correlation matrix.
 *
 * @param matrix
 * @param size
 */
std::vector<float> RMT::computeEigenvalues(std::vector<float>* matrix, size_t size)
{
    EDEBUG_FUNC(this,matrix,&size);

    // initialize helper variables
    int n = size;
    int lda = size;

    // determine whether a device is available
    if ( Ace::Settings::instance().cudaDevicePointer() )
    {
        // initialize cuSOLVER handle
        cusolverDnHandle_t cusolver_handle;

        CHECK_CUSOLVER(cusolverDnCreate(&cusolver_handle));

        // allocate device buffers
        ::CUDA::Buffer<float> matrixBuffer(n * n);
        ::CUDA::Buffer<float> eigensBuffer(n);
        ::CUDA::Buffer<int> info(1);

        // copy pruned matrix to device
        memcpy(matrixBuffer.hostData(), matrix->data(), matrix->size() * sizeof(float));
        matrixBuffer.write();

        // determine the size of the workspace
        int lwork;

        CHECK_CUSOLVER(cusolverDnSsyevd_bufferSize(
            cusolver_handle,
            CUSOLVER_EIG_MODE_NOVECTOR,
            CUBLAS_FILL_MODE_UPPER,
            n, (float *)matrixBuffer.deviceData(), lda,
            (float *)eigensBuffer.deviceData(),
            &lwork
        ));

        // allocate device buffer for workspace
        ::CUDA::Buffer<float> work(lwork);

        // compute eigenvalues
        CHECK_CUSOLVER(cusolverDnSsyevd(
            cusolver_handle,
            CUSOLVER_EIG_MODE_NOVECTOR,
            CUBLAS_FILL_MODE_UPPER,
            n, (float *)matrixBuffer.deviceData(), lda,
            (float *)eigensBuffer.deviceData(),
            (float *)work.deviceData(), lwork,
            (int *)info.deviceData()
        ));

        // destroy cuSOLVER handle
        CHECK_CUSOLVER(cusolverDnDestroy(cusolver_handle));

        // copy results to host
        std::vector<float> eigens(n);

        info.read().wait();
        eigensBuffer.read().wait();
        memcpy(eigens.data(), eigensBuffer.hostData(), eigens.size() * sizeof(float));

        // print warning if cuSOLVER returned error code
        if ( info.at(0) != 0 )
        {
            qInfo("warning: cuSOLVER ssyev returned %d", info.at(0));
        }

        // return eigenvalues
        return eigens;
    }
    else
    {
        // initialize eigenvalues and workspace
        std::vector<float> eigens(n);
        std::vector<float> work(5 * n);

        // compute eigenvalues
        int info = LAPACKE_ssyev_work(
            LAPACK_COL_MAJOR, 'N', 'U',
            n, matrix->data(), lda,
            eigens.data(),
            work.data(), work.size());

        // print warning if LAPACKE returned error code
        if ( info != 0 )
        {
            qInfo("warning: LAPACKE ssyev returned %d", info);
        }

        // return eigenvalues
        return eigens;
    }
}



/*!
 * Return the unique values of a sorted list of real numbers. Two real numbers
 * are unique if their absolute difference is greater than some small value
 * epsilon.
 *
 * @param values
 */
std::vector<float> RMT::computeUnique(const std::vector<float>& values)
{
    EDEBUG_FUNC(this,&values);

    std::vector<float> unique;

    for ( size_t i = 0; i < values.size(); ++i )
    {
        if ( unique.empty() || fabs(values.at(i) - unique.back()) > _uniqueEpsilon )
        {
            unique.push_back(values.at(i));
        }
    }

    return unique;
}



/*!
 * Compute the chi-squared test for the nearest-neighbor spacing distribution
 * (NNSD) of a list of eigenvalues. The list should be sorted and should contain
 * only unique values. If spline interpolation is enabled, the chi-squared value
 * is an average of several chi-squared tests, in which splines of varying pace
 * are applied to the eigenvalues. Otherwise, a single chi-squared test is
 * performed directly on the eigenvalues.
 *
 * @param eigens
 */
float RMT::computeChiSquare(const std::vector<float>& eigens)
{
    EDEBUG_FUNC(this,&eigens);

    // make sure there are enough eigenvalues
    if ( eigens.size() < (size_t) _minUniqueEigenvalues )
    {
        return -1;
    }

    // determine whether spline interpolation is enabled
    if ( _splineInterpolation )
    {
        // perform several chi-squared tests with spline interpolation by varying the pace
        float chi {0.0};
        int chiTestCount {0};

        for ( int pace = _minSplinePace; pace <= _maxSplinePace; ++pace )
        {
            // perform test only if there are enough eigenvalues for pace
            if ( eigens.size() / pace < 5 )
            {
                break;
            }

            // compute spline-interpolated eigenvalues
            std::vector<float> splineEigens {computeSpline(eigens, pace)};

            // compute chi-squared value
            float chiPace {computeChiSquareHelper(splineEigens)};

            qInfo("pace: %d, chi-squared: %g", pace, chiPace);

            // append chi-squared value to running sum
            chi += chiPace;
            ++chiTestCount;
        }

        // return average of chi-squared tests
        return chi / chiTestCount;
    }
    else
    {
        // perform a single chi-squared test without spline interpolation
        return computeChiSquareHelper(eigens);
    }
}



/*!
 * Compute the chi-squared test for the nearest-neighbor spacing distribution
 * (NNSD) of a list of values. The list should be sorted and should contain only
 * unique values.
 *
 * @param values
 */
float RMT::computeChiSquareHelper(const std::vector<float>& values)
{
    EDEBUG_FUNC(this,&values);

    // compute spacings
    std::vector<float> spacings {computeSpacings(values)};

    // compute histogram of spacings
    const float histogramMin {0};
    const float histogramMax {3};
    const float histogramBinWidth {(histogramMax - histogramMin) / _histogramBinSize};
    std::vector<float> histogram(_histogramBinSize);

    for ( auto& spacing : spacings )
    {
        if ( histogramMin <= spacing && spacing < histogramMax )
        {
            size_t index = static_cast<size_t>((spacing - histogramMin) / histogramBinWidth);

            ++histogram[index];
        }
    }

    // compute chi-squared value from the histogram
    float chi {0.0};

    for ( int i = 0; i < (int) histogram.size(); ++i )
    {
        // compute O_i, the number of elements in bin i
        float O_i {histogram[i]};

        // compute E_i, the expected value of Poisson distribution for bin i
        float E_i {(expf(-i * histogramBinWidth) - expf(-(i + 1) * histogramBinWidth)) * values.size()};

        // update chi-squared value based on difference between O_i and E_i
        chi += (O_i - E_i) * (O_i - E_i) / E_i;
    }

    return chi;
}



/*!
 * Compute a spline interpolation of a list of values using the given pace. The
 * list should be sorted and should contain only unique values. The pace determines
 * the ratio of values which are used as points to create the spline; for example,
 * a pace of 10 means that every 10th value is used to create the spline.
 *
 * @param values
 * @param pace
 */
std::vector<float> RMT::computeSpline(const std::vector<float>& values, int pace)
{
    EDEBUG_FUNC(this,&values,pace);

    // using declarations for gsl resource pointers
    using gsl_interp_accel_ptr = unique_ptr<gsl_interp_accel, decltype(&gsl_interp_accel_free)>;
    using gsl_spline_ptr = unique_ptr<gsl_spline, decltype(&gsl_spline_free)>;

    // extract values for spline based on pace
    int splineSize {(int) values.size() / pace};
    unique_ptr<double[]> x(new double[splineSize]);
    unique_ptr<double[]> y(new double[splineSize]);

    for ( int i = 0; i < splineSize; ++i )
    {
        x[i] = (double)values.at(i*pace);
        y[i] = (double)(i*pace + 1) / values.size();
    }
    x[splineSize - 1] = values.back();
    y[splineSize - 1] = 1.0;

    // initialize gsl spline
    gsl_interp_accel_ptr interp(gsl_interp_accel_alloc(), &gsl_interp_accel_free);
    gsl_spline_ptr spline(gsl_spline_alloc(gsl_interp_akima, splineSize), &gsl_spline_free);
    gsl_spline_init(spline.get(), x.get(), y.get(), splineSize);

    // extract interpolated values from spline
    std::vector<float> splineValues(values.size());

    splineValues[0] = 0.0;
    splineValues[values.size() - 1] = 1.0;

    for ( size_t i = 1; i < values.size() - 1; ++i )
    {
        splineValues[i] = static_cast<float>(gsl_spline_eval(spline.get(), values.at(i), interp.get()));
    }

    // return interpolated values
    return splineValues;
}



/*!
 * Compute the spacings of a list of values. The list should be sorted and should
 * contain only unique values.
 *
 * @param values
 */
std::vector<float> RMT::computeSpacings(const std::vector<float>& values)
{
    EDEBUG_FUNC(this,&values);

    std::vector<float> spacings(values.size() - 1);

    for ( size_t i = 0; i < spacings.size(); ++i )
    {
        spacings[i] = (values.at(i + 1) - values.at(i)) * values.size();
    }

    return spacings;
}
