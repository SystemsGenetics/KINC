#include <memory>
#include <random>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <lapacke.h>
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
   case ThresholdStart: return Type::Double;
   case ThresholdStep: return Type::Double;
   case ThresholdStop: return Type::Double;
   case MinUnfoldingPace: return Type::Integer;
   case MaxUnfoldingPace: return Type::Integer;
   case HistogramBinSize: return Type::Integer;
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
      case ThresholdStart: return QString("tstart");
      case ThresholdStep: return QString("tstep");
      case ThresholdStop: return QString("tstop");
      case MinUnfoldingPace: return QString("minpace");
      case MaxUnfoldingPace: return QString("maxpace");
      case HistogramBinSize: return QString("bins");
      default: return QVariant();
      }
   case Role::Title:
      // figure out which argument is being queried and return title
      switch (argument)
      {
      case InputData: return tr("Input:");
      case LogFile: return tr("Log File:");
      case ThresholdStart: return tr("Threshold Start:");
      case ThresholdStep: return tr("Threshold Step:");
      case ThresholdStop: return tr("Threshold Stop:");
      case MinUnfoldingPace: return tr("Minimum Unfolding Pace:");
      case MaxUnfoldingPace: return tr("Maximum Unfolding Pace:");
      case HistogramBinSize: return tr("Histogram Bin Size:");
      default: return QVariant();
      }
   case Role::WhatsThis:
      // figure out which argument is being queried and return "What's This?" text
      switch (argument)
      {
      case InputData: return tr("Correlation matrix for which an appropriate correlation threshold will be found.");
      case LogFile: return tr("Output text file that logs all results.");
      case ThresholdStart: return tr("Starting threshold.");
      case ThresholdStep: return tr("Threshold step size.");
      case ThresholdStop: return tr("Stopping threshold.");
      case MinUnfoldingPace: return tr("The minimum pace with which to perform unfolding.");
      case MaxUnfoldingPace: return tr("The maximum pace with which to perform unfolding.");
      case HistogramBinSize: return tr("The number of bins for the nearest-neighbor spacing histogram.");
      default: return QVariant();
      }
   case Role::DefaultValue:
      // figure out which argument is being queried and if applicable return default value else
      // return nothing
      switch (argument)
      {
      case ThresholdStart: return 0.99;
      case ThresholdStep: return 0.001;
      case ThresholdStop: return 0.5;
      case MinUnfoldingPace: return 10;
      case MaxUnfoldingPace: return 40;
      case HistogramBinSize: return 60;
      default: return QVariant();
      }
   case Role::Minimum:
      // figure out which argument is being queried and if applicable return minimum value else
      // return nothing
      switch (argument)
      {
      case ThresholdStart: return 0;
      case ThresholdStep: return 0;
      case ThresholdStop: return 0;
      case MinUnfoldingPace: return 1;
      case MaxUnfoldingPace: return 1;
      case HistogramBinSize: return 1;
      default: return QVariant();
      }
   case Role::Maximum:
      // figure out which argument is being queried and if applicable return maximum value else
      // return nothing
      switch (argument)
      {
      case ThresholdStart: return 1;
      case ThresholdStep: return 1;
      case ThresholdStop: return 1;
      case MinUnfoldingPace: return INT_MAX;
      case MaxUnfoldingPace: return INT_MAX;
      case HistogramBinSize: return INT_MAX;
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






void RMT::setArgument(int argument, QVariant value)
{
   // figure out which argument is being set and set it
   switch (argument)
   {
   case ThresholdStart:
      _thresholdStart = value.toDouble();
      break;
   case ThresholdStep:
      _thresholdStep = value.toDouble();
      break;
   case ThresholdStop:
      _thresholdStop = value.toDouble();
      break;
   case MinUnfoldingPace:
      _minUnfoldingPace = value.toInt();
      break;
   case MaxUnfoldingPace:
      _maxUnfoldingPace = value.toInt();
      break;
   case HistogramBinSize:
      _histogramBinSize = value.toInt();
      break;
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
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Argument Error"));
      e.setDetails(QObject::tr("Did not get valid input or logfile arguments."));
      throw e;
   }

   // make sure threshold arguments are valid
   if ( _thresholdStart <= _thresholdStop )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Argument Error"));
      e.setDetails(QObject::tr("Starting threshold must be greater than stopping threshold."));
      throw e;
   }

   // make sure pace arguments are valid
   if ( _minUnfoldingPace >= _maxUnfoldingPace )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Argument Error"));
      e.setDetails(QObject::tr("Minimum unfolding pace must be less than maximum unfolding pace."));
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
   int totalSteps = (_thresholdStart - _thresholdStop) / _thresholdStep;

   // continue while max chi is less than final threshold
   while ( maxChi < _chiSquareThreshold2 )
   {
      // make sure interruption is not requested
      if ( isInterruptionRequested() )
      {
         return;
      }

      qInfo("\n");
      qInfo("threshold: %g", threshold);

      // compute pruned matrix based on threshold
      int size;
      QVector<float> pruneMatrix {computePruneMatrix(threshold, &size)};

      qInfo("prune matrix: %d", size);

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

         qInfo("eigenvalues: %d", eigens.size());

         // make sure interruption is not requested
         if ( isInterruptionRequested() )
         {
            return;
         }

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
         e.setTitle(QObject::tr("RMT Threshold Error"));
         e.setDetails(QObject::tr("Could not find non-random threshold above stopping threshold."));
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
   _maxCorrelations.fill(0, _input->geneSize() * _input->maxClusterSize());

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
         int index1 = pair.index().getX() * _input->maxClusterSize() + cluster;
         int index2 = pair.index().getY() * _input->maxClusterSize() + cluster;
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






QVector<float> RMT::computePruneMatrix(float threshold, int* size)
{
   // generate vector of cluster indexes that have a correlation above threshold
   int numClusters = 0;
   QVector<int> indices(_input->geneSize() * _input->maxClusterSize(), -1);

   for ( int i = 0; i < indices.size(); ++i )
   {
      if ( _maxCorrelations.at(i) >= threshold )
      {
         indices[i] = numClusters;
         ++numClusters;
      }
   }

   // allocate new pruned matrix with cluster size
   QVector<float> pruneMatrix(numClusters * numClusters);

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
            int i = indices[pair.index().getX() * _input->maxClusterSize() + cluster];
            int j = indices[pair.index().getY() * _input->maxClusterSize() + cluster];

            pruneMatrix[i * numClusters + j] = correlation;
            pruneMatrix[j * numClusters + i] = correlation;
         }
      }
   }

   // set size to size of pruned matrix and return pointer
   *size = numClusters;
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
