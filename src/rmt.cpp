#include <memory>
#include <random>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <lapacke.h>

#include "rmt.h"
#include "rmt_input.h"
#include "correlationmatrix.h"
#include "datafactory.h"



using namespace std;






int RMT::size() const
{
   return 1;
}






void RMT::process(const EAbstractAnalytic::Block* /*result*/)
{
   // initialize log text stream
   QTextStream stream(_logfile);

   // initialize helper variables
   float finalThreshold {0};
   float finalChi {INFINITY};
   float maxChi {-INFINITY};

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

         // compute chi-square value from NNSD of eigenvalues
         chi = computeChiSquare(eigens);

         qInfo("chi-square: %g", chi);
      }

      // make sure that chi-square test succeeded
      if ( chi != -1 )
      {
         // save the most recent chi-square value less than critical value
         if ( chi < _chiSquareThreshold1 )
         {
            finalChi = chi;
            finalThreshold = threshold;
         }

         // save the largest chi-square value which occurs after finalChi
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






EAbstractAnalytic::Input* RMT::makeInput()
{
   return new Input(this);
}






void RMT::initialize()
{
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
   if ( _minUnfoldingPace >= _maxUnfoldingPace )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Minimum unfolding pace must be less than maximum unfolding pace."));
      throw e;
   }
}






QVector<float> RMT::computeMaximums(const QVector<float>& matrix)
{
   const int N {_input->geneSize()};
   const int K {_input->maxClusterSize()};

   // initialize elements to minimum value
   QVector<float> maximums(N * K, 0);

   // compute maximum of each row/column
   for ( int i = 0; i < N; ++i )
   {
      for ( int j = 0; j < i; ++j )
      {
         for ( int k = 0; k < K; ++k )
         {
            float correlation = fabs(matrix[i * N * K + j * K + k]);

            if ( maximums[i * K + k] < correlation )
            {
               maximums[i * K + k] = correlation;
            }

            if ( maximums[j * K + k] < correlation )
            {
               maximums[j * K + k] = correlation;
            }
         }
      }
   }

   return maximums;
}






QVector<float> RMT::computePruneMatrix(const QVector<float>& matrix, const QVector<float>& maximums, float threshold, int* size)
{
   const int N {_input->geneSize()};
   const int K {_input->maxClusterSize()};

   // generate vector of row/column indices that have a correlation above threshold
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
         if ( indices[i] % K != indices[j] % K )
         {
            continue;
         }

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






QVector<float> RMT::computeEigenvalues(QVector<float>* pruneMatrix, int size)
{
   QVector<float> eigens(size);
   QVector<float> work(5 * size);

   int info = LAPACKE_ssyev_work(
      LAPACK_COL_MAJOR, 'N', 'U',
      size, pruneMatrix->data(), size,
      eigens.data(),
      work.data(), work.size()
   );

   if ( info != 0 )
   {
      qInfo("warning: ssyev returned %d", info);
   }

   return eigens;
}






float RMT::computeChiSquare(const QVector<float>& eigens)
{
   // compute unique eigenvalues
   QVector<float> unique {degenerate(eigens)};

   qInfo("unique eigenvalues: %d", unique.size());

   // make sure there are enough unique eigenvalues
   if ( unique.size() < _minEigenvalueSize )
   {
      return -1;
   }

   // perform several chi-square tests by varying the pace
   float chi {0.0};
   int chiTestCount {0};

   for ( int pace = _minUnfoldingPace; pace <= _maxUnfoldingPace; ++pace )
   {
      // perform test only if there are enough eigenvalues for pace
      if ( unique.size() / pace < 5 )
      {
         break;
      }

      chi += computePaceChiSquare(unique, pace);
      ++chiTestCount;
   }

   // compute average of chi-square tests
   chi /= chiTestCount;

   // return chi value
   return chi;
}






float RMT::computePaceChiSquare(const QVector<float>& eigens, int pace)
{
   // compute eigenvalue spacings
   QVector<float> spacings {unfold(eigens, pace)};

   // compute nearest-neighbor spacing distribution
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

   // compute chi-square value from nearest-neighbor spacing distribution
   float chi {0.0};

   for ( int i = 0; i < histogram.size(); ++i )
   {
      // compute O_i, the number of elements in bin i
      float O_i {histogram[i]};

      // compute E_i, the expected value of Poisson distribution for bin i
      float E_i {(exp(-i * histogramBinWidth) - exp(-(i + 1) * histogramBinWidth)) * eigens.size()};

      // update chi-square value based on difference between O_i and E_i
      chi += (O_i - E_i) * (O_i - E_i) / E_i;
   }

   qInfo("pace: %d, chi: %g", pace, chi);

   return chi;
}






QVector<float> RMT::degenerate(const QVector<float>& eigens)
{
   const float EPSILON {1e-6};
   QVector<float> unique;

   for ( int i = 1; i < eigens.size(); ++i )
   {
      if ( unique.isEmpty() || fabs(eigens.at(i) - unique.last()) > EPSILON )
      {
         unique.append(eigens.at(i));
      }
   }

   return unique;
}






QVector<float> RMT::unfold(const QVector<float>& eigens, int pace)
{
   // using declarations for gsl resource pointers
   using gsl_interp_accel_ptr = unique_ptr<gsl_interp_accel, decltype(&gsl_interp_accel_free)>;
   using gsl_spline_ptr = unique_ptr<gsl_spline, decltype(&gsl_spline_free)>;

   // extract eigenvalues for spline based on pace
   int splineSize {eigens.size() / pace};
   unique_ptr<double[]> x(new double[splineSize]);
   unique_ptr<double[]> y(new double[splineSize]);

   for ( int i = 0; i < splineSize; ++i )
   {
      x[i] = (double)eigens.at(i*pace);
      y[i] = (double)(i*pace + 1) / eigens.size();
   }
   x[splineSize - 1] = eigens.back();
   y[splineSize - 1] = 1.0;

   // initialize gsl spline
   gsl_interp_accel_ptr interp(gsl_interp_accel_alloc(), &gsl_interp_accel_free);
   gsl_spline_ptr spline(gsl_spline_alloc(gsl_interp_akima, splineSize), &gsl_spline_free);
   gsl_spline_init(spline.get(), x.get(), y.get(), splineSize);

   // extract interpolated eigenvalues from spline
   QVector<float> splineEigens(eigens.size());

   splineEigens[0] = 0.0;
   splineEigens[eigens.size() - 1] = 1.0;

   for ( int i = 1; i < eigens.size() - 1; ++i )
   {
      splineEigens[i] = gsl_spline_eval(spline.get(), eigens.at(i), interp.get());
   }

   // compute spacings between interpolated eigenvalues
   QVector<float> spacings(eigens.size() - 1);

   for ( int i = 0; i < spacings.size(); ++i )
   {
      spacings[i] = (splineEigens.at(i + 1) - splineEigens.at(i)) * eigens.size();
   }

   return spacings;
}
