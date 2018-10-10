#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

#include "rmt.h"
#include "rmt_input.h"
#include "correlationmatrix.h"



using namespace std;






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
void RMT::process(const EAbstractAnalytic::Block*)
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
   QVector<float> matrix {_input->dumpRawData()};
   QVector<float> maximums {computeMaximums(matrix)};

   // continue while max chi is less than final threshold
   while ( maxChi < _chiSquareThreshold2 )
   {
      qInfo("\n");
      qInfo("threshold: %g", threshold);

      // compute pruned matrix based on threshold
      int size;
      QVector<float> pruneMatrix {computePruneMatrix(matrix, maximums, threshold, &size)};

      qInfo("prune matrix: %d", size);

      // make sure that pruned matrix is not empty
      float chi = -1;

      if ( size > 0 )
      {
         // compute eigenvalues of pruned matrix
         QVector<float> eigens {computeEigenvalues(&pruneMatrix, size)};

         qInfo("eigenvalues: %d", eigens.size());

         // compute unique eigenvalues
         QVector<float> unique {computeUnique(eigens)};

         qInfo("unique eigenvalues: %d", unique.size());

         // compute chi-squared value from NNSD of eigenvalues
         chi = computeChiSquare(unique);

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
      stream << threshold << "\t" << size << "\t" << chi << "\n";

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
EAbstractAnalytic::Input* RMT::makeInput()
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
}






/*!
 * Compute the row-wise maximums of a correlation matrix.
 *
 * @param matrix
 */
QVector<float> RMT::computeMaximums(const QVector<float>& matrix)
{
   EDEBUG_FUNC(this,&matrix);

   const int N {_input->geneSize()};
   const int K {_input->maxClusterSize()};

   // initialize elements to minimum value
   QVector<float> maximums(N * K, 0);

   // compute maximum correlation of each row
   for ( int i = 0; i < N; ++i )
   {
      for ( int j = 0; j < N; ++j )
      {
         for ( int k = 0; k < K; ++k )
         {
            float correlation = fabs(matrix[i * N * K + j * K + k]);

            if ( maximums[i * K + k] < correlation )
            {
               maximums[i * K + k] = correlation;
            }
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
 * @param matrix
 * @param maximums
 * @param threshold
 * @param size
 */
QVector<float> RMT::computePruneMatrix(const QVector<float>& matrix, const QVector<float>& maximums, float threshold, int* size)
{
   EDEBUG_FUNC(this,&matrix,&maximums,threshold,size);

   const int N {_input->geneSize()};
   const int K {_input->maxClusterSize()};

   // generate vector of row indices that have a correlation above threshold
   QVector<int> indices;

   for ( int i = 0; i < maximums.size(); ++i )
   {
      if ( maximums[i] >= threshold )
      {
         indices.append(i);
      }
   }

   // extract pruned matrix from correlation matrix
   QVector<float> pruneMatrix(indices.size() * indices.size());

   for ( int i = 0; i < indices.size(); ++i )
   {
      for ( int j = 0; j < i; ++j )
      {
         // make sure that i and j refer to the same cluster number
         if ( indices[i] % K != indices[j] % K )
         {
            continue;
         }

         // save correlation if it is above threshold
         float correlation = matrix[indices[i]/K * N * K + indices[j]/K * K + indices[i] % K];

         if ( fabs(correlation) >= threshold )
         {
            pruneMatrix[i * indices.size() + j] = correlation;
         }
      }

      pruneMatrix[i * indices.size() + i] = 1;
   }

   // save size of pruned matrix
   *size = indices.size();

   return pruneMatrix;
}






/*!
 * Compute the eigenvalues of a correlation matrix.
 *
 * @param matrix
 * @param size
 */
QVector<float> RMT::computeEigenvalues(QVector<float>* matrix, int size)
{
   EDEBUG_FUNC(this,matrix,size);

   // using declarations for gsl resources
   using gsl_vector_ptr = unique_ptr<gsl_vector,decltype(&gsl_vector_free)>;
   using gsl_matrix_ptr = unique_ptr<gsl_matrix,decltype(&gsl_matrix_free)>;
   using gsl_eigen_symmv_workspace_ptr = unique_ptr<gsl_eigen_symmv_workspace
      ,decltype(&gsl_eigen_symmv_free)>;

   QVector<double> temp;
   for (auto val: *matrix) temp.append(val);

   // make and initialize gsl eigen resources
   gsl_matrix_view view = gsl_matrix_view_array(temp.data(),size,size);
   gsl_vector_ptr eval (gsl_vector_alloc(size),&gsl_vector_free);
   gsl_matrix_ptr evec (gsl_matrix_alloc(size,size),&gsl_matrix_free);
   gsl_eigen_symmv_workspace_ptr work (gsl_eigen_symmv_alloc(size),&gsl_eigen_symmv_free);

   // have gsl compute eigen values for the pruned matrix
   gsl_eigen_symmv(&view.matrix,eval.get(),evec.get(),work.get());
   gsl_eigen_symmv_sort(eval.get(),evec.get(),GSL_EIGEN_SORT_VAL_ASC);

   // create return vector and get eigen values from gsl
   QVector<float> ret(size);
   for (int i = 0; i < size ;i++)
   {
      ret[i] = gsl_vector_get(eval.get(),i);
   }

   // return eigen values vector
   return ret;
}






/*!
 * Return the unique values of a sorted list of real numbers. Two real numbers
 * are unique if their absolute difference is greater than some small value
 * epsilon.
 *
 * @param values
 */
QVector<float> RMT::computeUnique(const QVector<float>& values)
{
   EDEBUG_FUNC(this,&values);

   const float EPSILON {1e-6};
   QVector<float> unique;

   for ( int i = 1; i < values.size(); ++i )
   {
      if ( unique.isEmpty() || fabs(values.at(i) - unique.last()) > EPSILON )
      {
         unique.append(values.at(i));
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
float RMT::computeChiSquare(const QVector<float>& eigens)
{
   EDEBUG_FUNC(this,&eigens);

   // make sure there are enough eigenvalues
   if ( eigens.size() < _minEigenvalueSize )
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
         QVector<float> splineEigens {computeSpline(eigens, pace)};

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
float RMT::computeChiSquareHelper(const QVector<float>& values)
{
   EDEBUG_FUNC(this,&values);

   // compute spacings
   QVector<float> spacings {computeSpacings(values)};

   // compute histogram of spacings
   const float histogramMin {0};
   const float histogramMax {3};
   const float histogramBinWidth {(histogramMax - histogramMin) / _histogramBinSize};
   QVector<float> histogram(_histogramBinSize);

   for ( auto& spacing : spacings )
   {
      if ( histogramMin <= spacing && spacing < histogramMax )
      {
         ++histogram[(spacing - histogramMin) / histogramBinWidth];
      }
   }

   // compute chi-squared value from the histogram
   float chi {0.0};

   for ( int i = 0; i < histogram.size(); ++i )
   {
      // compute O_i, the number of elements in bin i
      float O_i {histogram[i]};

      // compute E_i, the expected value of Poisson distribution for bin i
      float E_i {(exp(-i * histogramBinWidth) - exp(-(i + 1) * histogramBinWidth)) * values.size()};

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
QVector<float> RMT::computeSpline(const QVector<float>& values, int pace)
{
   EDEBUG_FUNC(this,&values,pace);

   // using declarations for gsl resource pointers
   using gsl_interp_accel_ptr = unique_ptr<gsl_interp_accel, decltype(&gsl_interp_accel_free)>;
   using gsl_spline_ptr = unique_ptr<gsl_spline, decltype(&gsl_spline_free)>;

   // extract eigenvalues for spline based on pace
   int splineSize {values.size() / pace};
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

   // extract interpolated eigenvalues from spline
   QVector<float> splineValues(values.size());

   splineValues[0] = 0.0;
   splineValues[values.size() - 1] = 1.0;

   for ( int i = 1; i < values.size() - 1; ++i )
   {
      splineValues[i] = gsl_spline_eval(spline.get(), values.at(i), interp.get());
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
QVector<float> RMT::computeSpacings(const QVector<float>& values)
{
   EDEBUG_FUNC(this,&values);

   // compute spacings between interpolated eigenvalues
   QVector<float> spacings(values.size() - 1);

   for ( int i = 0; i < spacings.size(); ++i )
   {
      spacings[i] = (values.at(i + 1) - values.at(i)) * values.size();
   }

   return spacings;
}
