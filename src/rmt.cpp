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
   int lastPercent {50};
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
      QVector<double> pruneMatrix {computePruneMatrix(threshold,&size)};

      // use chi-square value of zero if matrix is empty
      float chi = 0;

      if ( size > 0 )
      {
         // make sure interruption is not requested
         if ( isInterruptionRequested() )
         {
            return;
         }

         // compute eigenvalues of pruned matrix
         QVector<double> eigens {computeEigenvalues(&pruneMatrix,size)};

         // make sure interruption is not requested
         if ( isInterruptionRequested() )
         {
            return;
         }

         // compute chi-square value from NNSD of eigenvalues
         chi = computeChiSquare(&eigens);
      }

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

      // output to log file
      stream << threshold << "\t" << size << "\t" << chi << "\n";

      // decrement threshold by step and make sure minimum is not reached
      threshold -= _thresholdStep;
      if ( threshold < _thresholdMinimum )
      {
         // report no threshold could be found and fail
         E_MAKE_EXCEPTION(e);
         e.setTitle(QObject::tr("RMT Threshold Error"));
         e.setDetails(QObject::tr("Could not find non-random threshold above minimum."));
         throw e;
      }

      // determine new percent complete and check if it is new
      int newPercent {50 + 50*steps/totalSteps};
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
      qint64 newPercent {50*steps/totalSteps};

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
   // generate vector of cluster indexes that have max threshold above given threshold
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






QVector<double> RMT::computeEigenvalues(QVector<double>* pruneMatrix, int size)
{
   // using declarations for gsl resources
   using gsl_vector_ptr = unique_ptr<gsl_vector,decltype(&gsl_vector_free)>;
   using gsl_matrix_ptr = unique_ptr<gsl_matrix,decltype(&gsl_matrix_free)>;
   using gsl_eigen_symmv_workspace_ptr = unique_ptr<gsl_eigen_symmv_workspace
      ,decltype(&gsl_eigen_symmv_free)>;

   // make and initialize gsl eigen resources
   gsl_matrix_view view = gsl_matrix_view_array(pruneMatrix->data(),size,size);
   gsl_vector_ptr eval (gsl_vector_alloc(size),&gsl_vector_free);
   gsl_matrix_ptr evec (gsl_matrix_alloc(size,size),&gsl_matrix_free);
   gsl_eigen_symmv_workspace_ptr work (gsl_eigen_symmv_alloc(size),&gsl_eigen_symmv_free);

   // have gsl compute eigen values for the pruned matrix
   gsl_eigen_symmv(&view.matrix,eval.get(),evec.get(),work.get());
   gsl_eigen_symmv_sort(eval.get(),evec.get(),GSL_EIGEN_SORT_ABS_ASC);

   // create return vector and get eigen values from gsl
   QVector<double> ret(size);
   for (int i = 0; i < size ;i++)
   {
      ret[i] = gsl_vector_get(eval.get(),i);
   }

   // return eigen values vector
   return ret;
}






float RMT::computeChiSquare(QVector<double>* eigens)
{
   // initialize chi value and degenerate eigen vector
   float chi {0.0};
   degenerate(eigens);

   // make sure new size of eigen vector is minimum size
   if ( eigens->size() >= _minEigenvalueSize )
   {
      // initialize chi test count and iterate through all paces
      int chiTestCount {0};
      for (int pace = _minUnfoldingPace; pace <= _maxUnfoldingPace ;++pace)
      {
         // if ??? increment chi value with paced chi square value and increment test count
         if ( (eigens->size()/pace) >= 5 )
         {
            chi += computePaceChiSquare(*eigens,pace);
            ++chiTestCount;
         }
      }

      // divide chi value by amount of chi tests added to it
      chi /= chiTestCount;
   }

   // if chi is not a real number set it to zero
   if ( isnan(chi) || isinf(chi) )
   {
      chi = 0.0;
   }

   // return chi value
   return chi;
}






float RMT::computePaceChiSquare(const QVector<double>& eigens, int pace)
{
   // initialize chi value that will be returned and get unfolded vector
   float chi {0.0};
   QVector<double> unfolded {unfold(eigens,pace)};

   // iterate through each bin to count matches
   for (int i = 0; i < ((int)(3.0/_chiSquareBinSize) + 1) ;++i)
   {
      // initialize count and iterate through all unfolded values
      int count {0};
      for (const auto& unfold : unfolded)
      {
         // if unfold value is within bin range increment count
         if ( unfold < (i + 1)*_chiSquareBinSize && unfold > i*_chiSquareBinSize )
         {
            ++count;
         }
      }

      // calculate expected value and add to chi based off difference from expected count and actual
      float expect {(exp(-1*i*_chiSquareBinSize)
               - exp(-1*(i + 1)*_chiSquareBinSize))*eigens.size()};
      chi += ((double)count - expect)*((double)count - expect)/expect;
   }

   // return chi value
   return chi;
}






void RMT::degenerate(QVector<double>* eigens)
{
   // iterate through all eigen values
   for (auto& eigen : *eigens)
   {
      // if eigen value is less than minimum set it to zero
      if ( fabs(eigen) < _minEigenvalue )
      {
         eigen = 0.0;
      }
   }

   // sort all eigen values
   sort(eigens->begin(),eigens->end());

   // iterate through eigen values until second to last is reached
   int i {0};
   while ( (i + 1) < eigens->size() )
   {
      // if next eigen value equals current eigen value remove current value
      if ( eigens->at(i) == eigens->at(i + 1) )
      {
         eigens->removeAt(i);
      }

      // else iterate to next value
      else
      {
         ++i;
      }
   }
}






QVector<double> RMT::unfold(const QVector<double>& eigens, int pace)
{
   // using declarations for gsl resource pointers
   using gsl_interp_accel_ptr = unique_ptr<gsl_interp_accel,decltype(&gsl_interp_accel_free)>;
   using gsl_spline_ptr = unique_ptr<gsl_spline,decltype(&gsl_spline_free)>;

   // determine pace size and create x and y arrays
   int paceSize {eigens.size()/pace};
   unique_ptr<double[]> paceXs(new double[paceSize]);
   unique_ptr<double[]> paceYs(new double[paceSize]);

   // populate x and y arrays
   for (int i = 0; i < paceSize ;++i)
   {
      paceXs[i] = eigens.at(i*pace);
      paceYs[i] = (i*pace + 1)/(double)eigens.size();
   }
   paceXs[paceSize - 1] = eigens.back();
   paceYs[paceSize - 1] = 1.0;

   // create resources and initialize gsl spline
   gsl_interp_accel_ptr interp(gsl_interp_accel_alloc(),&gsl_interp_accel_free);
   gsl_spline_ptr spline(gsl_spline_alloc(gsl_interp_akima,paceSize),&gsl_spline_free);
   gsl_spline_init(spline.get(),paceXs.get(),paceYs.get(),paceSize);

   // create unfolding vector
   QVector<double> ret;
   ret.reserve(eigens.size());
   ret.push_back(0.0);

   // populate unfolding vector from gsl spline values
   for (int i = 1; i < (eigens.size() - 1) ;++i)
   {
      ret.push_back(gsl_spline_eval(spline.get(),eigens.at(i),interp.get()));
   }
   ret.push_back(1.0);

   // rewrite unfolding vector values as distance between values popping the very last value
   for (int i = 0; i < (eigens.size() - 1) ;++i)
   {
      ret[i] = (ret.at(i + 1) - ret.at(i))*eigens.size();
   }
   ret.pop_back();

   // sort unfolding vector and return it
   sort(ret.begin(),ret.end());
   return ret;
}
