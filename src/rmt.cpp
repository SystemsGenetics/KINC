#include <memory>
#include <random>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <ace/core/metadata.h>

#include "rmt.h"
#include "correlationmatrix.h"
#include "datafactory.h"
#include <iostream>



using namespace std;






EAbstractAnalytic::ArgumentType RMT::getArgumentData(int argument)
{
   // use type declaration
   using Type = EAbstractAnalytic::ArgumentType;

   // figure out which argument is being queried and return its type
   switch (argument)
   {
   case InputData: return Type::DataIn;
   case LogFile: return Type::FileOut;
   default: return Type::Bool;
   }
}






QVariant RMT::getArgumentData(int argument, EAbstractAnalytic::Role role)
{
   // use role declaration
   using Role = EAbstractAnalytic::Role;

   // figure out which role is being queried
   switch (role)
   {
   case Role::CommandLineName:
      // figure out which argument is being queried and return command line name
      switch (argument)
      {
      case InputData: return QString("input");
      case LogFile: return QString("log");
      default: return QVariant();
      }
   case Role::Title:
      // figure out which argument is being queried and return title
      switch (argument)
      {
      case InputData: return tr("Input:");
      case LogFile: return tr("Log File:");
      default: return QVariant();
      }
   case Role::WhatsThis:
      // figure out which argument is being queried and return "What's This?" text
      switch (argument)
      {
      case InputData: return tr("Input correlation matrix that will be used to output correlation"
                                " above a threshold determined by Random Matrix Theory.");
      case LogFile: return tr("Output text file that logs all chi trials.");
      default: return QVariant();
      }
   case Role::DataType:
      // see if this is input data and return type else return nothing
      switch (argument)
      {
      case InputData: return DataFactory::CorrelationMatrixType;
      default: return QVariant();
      }
   case Role::FileFilters:
      // figure out which argument is being queried and return file filters
      switch (argument)
      {
      case LogFile: return tr("Text File %1").arg("(*.txt)");
      default: return QVariant();
      }
   default:
      return QVariant();
   }
}






void RMT::setArgument(int argument, QFile* file)
{
   // figure out which argument is being queried and set file pointer
   switch (argument)
   {
   case LogFile:
      _logfile = file;
      break;
   }
}






void RMT::setArgument(int argument, EAbstractData* data)
{
   // figure out which argument is being queried and set data pointer
   switch (argument)
   {
   case InputData:
      _input = dynamic_cast<CorrelationMatrix*>(data);
      break;
   }
}






bool RMT::initialize()
{
   // make sure input and output were set properly
   if ( !_input || !_logfile )
   {
      // report argument error and fail
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Argument Error"));
      e.setDetails(QObject::tr("Did not get valid input or logfile arguments."));
      throw e;
   }

   // nothing to pre-allocate
   return false;
}






void RMT::runSerial()
{
   // generate list of maximum correlations for each gene
   computeGeneThresholds();

   // initialize log text stream
   QTextStream stream(_logfile);

   // initialize helper variables
   float finalThreshold {0};
   float finalChi {INFINITY};
   float maxChi {-INFINITY};

   float threshold {_thresholdStart};

   // initialize last percent, steps, and total steps
   int lastPercent {10};
   int steps {0};
   int totalSteps = (_thresholdStart - _thresholdMinimum) / _thresholdStep;

   // continue while max chi is less than final threshold
   while ( maxChi < _chiSquareThreshold2 )
   {
      // make sure interruption is not requested
      if ( isInterruptionRequested() )
      {
         return;
      }

      // compute pruned matrix based on threshold
      int size;
      QVector<double> pruneMatrix {computePruneMatrix(threshold, &size)};

      // make sure that pruned matrix is not empty
      float chi = -1;

      if ( size > 0 )
      {
         // make sure interruption is not requested
         if ( isInterruptionRequested() )
         {
            return;
         }

         // compute eigenvalues of pruned matrix
         QVector<float> eigens {computeEigenvalues(&pruneMatrix, size)};

         // make sure interruption is not requested
         if ( isInterruptionRequested() )
         {
            return;
         }

         // compute chi-square value from NNSD of eigenvalues
         chi = computeChiSquare(&eigens);
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
      if ( threshold < _thresholdMinimum )
      {
         E_MAKE_EXCEPTION(e);
         e.setTitle(QObject::tr("RMT Threshold Error"));
         e.setDetails(QObject::tr("Could not find non-random threshold above minimum."));
         throw e;
      }

      // determine new percent complete and check if it is new
      int newPercent {10 + 90*steps/totalSteps};
      if ( newPercent != lastPercent )
      {
         // update new percentage and emit progressed signal
         lastPercent = newPercent;
         emit progressed(lastPercent);
      }

      // increment steps taken
      ++steps;
   }

   // write threshold where chi was first above final threshold
   stream << finalThreshold << "\n";
}






void RMT::computeGeneThresholds()
{
   // initialize percent complete and steps
   int lastPercent {0};
   qint64 steps {0};
   qint64 totalSteps {_input->size()};

   // initialize elements to minimum value
   _maxCorrelations.fill(0, _input->geneSize() * _maxClusterSize);

   // iterate through all gene pairs
   CorrelationMatrix::Pair pair(_input);

   while ( pair.hasNext() )
   {
      // make sure interruption is not requested
      if ( isInterruptionRequested() )
      {
         return;
      }

      // read in next gene pair
      pair.readNext();

      // iterate through each cluster
      for ( int cluster = 0; cluster < pair.clusterSize(); ++cluster )
      {
         int index1 = pair.index().getX() * _maxClusterSize + cluster;
         int index2 = pair.index().getY() * _maxClusterSize + cluster;
         float correlation = fabs(pair.at(cluster, 0));

         // if value is greater than current max of gene x set it to new max
         if ( _maxCorrelations[index1] < correlation )
         {
            _maxCorrelations[index1] = correlation;
         }

         // if value is greater than current max of gene y set it to new max
         if ( _maxCorrelations[index2] < correlation )
         {
            _maxCorrelations[index2] = correlation;
         }
      }

      // increment steps and compute new percent complete
      ++steps;
      qint64 newPercent {10*steps/totalSteps};

      // check to see if percent complete changed
      if ( newPercent != lastPercent )
      {
         // update percent complete and emit progressed signal
         lastPercent = newPercent;
         emit progressed(lastPercent);
      }
   }
}






QVector<double> RMT::computePruneMatrix(float threshold, int* size)
{
   // generate vector of cluster indexes that have a correlation above threshold
   int numClusters = 0;
   QVector<int> indices(_input->geneSize() * _maxClusterSize, -1);

   for ( int i = 0; i < indices.size(); ++i )
   {
      if ( _maxCorrelations.at(i) >= threshold )
      {
         indices[i] = numClusters;
         ++numClusters;
      }
   }

   // allocate new pruned matrix with cluster size
   QVector<double> pruneMatrix(numClusters * numClusters);

   // initialize diagonal
   for ( int i = 0; i < numClusters; ++i )
   {
      pruneMatrix[i * numClusters + i] = 1;
   }

   // iterate through all gene pairs
   CorrelationMatrix::Pair pair(_input);

   while ( pair.hasNext() )
   {
      // read in next gene pair
      pair.readNext();

      // iterate through each cluster
      for ( int cluster = 0; cluster < pair.clusterSize(); ++cluster )
      {
         float correlation = pair.at(cluster, 0);

         if ( !isnan(correlation) && fabs(correlation) >= threshold )
         {
            int i = indices[pair.index().getX() * _maxClusterSize + cluster];
            int j = indices[pair.index().getY() * _maxClusterSize + cluster];

            pruneMatrix[i * numClusters + j] = correlation;
         }
      }
   }

   // set size to size of pruned matrix and return pointer
   *size = numClusters;
   return pruneMatrix;
}






QVector<float> RMT::computeEigenvalues(QVector<double>* pruneMatrix, int size)
{
   // using declarations for gsl resources
   using gsl_vector_ptr = unique_ptr<gsl_vector, decltype(&gsl_vector_free)>;
   using gsl_matrix_ptr = unique_ptr<gsl_matrix, decltype(&gsl_matrix_free)>;
   using gsl_eigen_symmv_workspace_ptr = unique_ptr<gsl_eigen_symmv_workspace, decltype(&gsl_eigen_symmv_free)>;

   // initialize gsl eigen resources
   gsl_matrix_view view = gsl_matrix_view_array(pruneMatrix->data(), size, size);
   gsl_vector_ptr eval (gsl_vector_alloc(size), &gsl_vector_free);
   gsl_matrix_ptr evec (gsl_matrix_alloc(size, size), &gsl_matrix_free);
   gsl_eigen_symmv_workspace_ptr work (gsl_eigen_symmv_alloc(size), &gsl_eigen_symmv_free);

   // compute eigenvalues for the pruned matrix
   gsl_eigen_symmv(&view.matrix, eval.get(), evec.get(), work.get());
   gsl_eigen_symmv_sort(eval.get(), evec.get(), GSL_EIGEN_SORT_ABS_ASC);

   // copy eigenvalues from gsl vector
   QVector<float> eigens(size);
   for (int i = 0; i < size; ++i )
   {
      eigens[i] = gsl_vector_get(eval.get(), i);
   }

   // return eigenvalues
   return eigens;
}






float RMT::computeChiSquare(QVector<float>* eigens)
{
   // remove duplicate eigenvalues
   degenerate(eigens);

   // make sure there are enough unique eigenvalues
   if ( eigens->size() < _minEigenvalueSize )
   {
      return -1;
   }

   // perform several chi-square tests by varying the pace
   float chi {0.0};
   int chiTestCount {0};

   for ( int pace = _minUnfoldingPace; pace <= _maxUnfoldingPace; ++pace )
   {
      // perform test only if there are enough eigenvalues for pace
      if ( eigens->size() / pace < 5 )
      {
         continue;
      }

      chi += computePaceChiSquare(*eigens, pace);
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
   QVector<float> histogram((int)((histogramMax - histogramMin) / _chiSquareBinSize))

   for ( auto& spacing : spacings )
   {
      if ( histogramMin <= spacing && spacing < histogramMax )
      {
         ++histogram[(spacing - histogramMin) / _chiSquareBinSize];
      }
   }

   // compute chi-square value from nearest-neighbor spacing distribution
   float chi {0.0};

   for ( int i = 0; i < histogram.size(); ++i )
   {
      // compute O_i, the number of elements in bin i
      float O_i {histogram[i]};

      // compute E_i, the expected value of Poisson distribution for bin i
      float E_i {(exp(-i * _chiSquareBinSize) - exp(-(i + 1) * _chiSquareBinSize)) * eigens.size()};

      // update chi-square value based on difference between O_i and E_i
      chi += (O_i - E_i) * (O_i - E_i) / E_i;
   }

   return chi;
}






void RMT::degenerate(QVector<float>* eigens)
{
   const float EPSILON {1e-6};

   // set extremely small eigenvalues to zero
   for ( auto& eigen : *eigens )
   {
      if ( fabs(eigen) < EPSILON )
      {
         eigen = 0.0;
      }
   }

   // sort all eigenvalues
   sort(eigens->begin(), eigens->end());

   // remove duplicate eigenvalues
   for ( int i = 0; i + 1 < eigens->size(); ++i )
   {
      if ( fabs(eigens->at(i) - eigens->at(i + 1)) < EPSILON )
      {
         eigens->removeAt(i);
         --i;
      }
   }
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
