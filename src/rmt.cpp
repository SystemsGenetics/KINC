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
   case OutputFile: return Type::FileOut;
   case LogFile: return Type::FileOut;
   case FilterSize: return Type::Integer;
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
      case OutputFile: return QString("output");
      case LogFile: return QString("log");
      case FilterSize: return QString("filter");
      default: return QVariant();
      }
   case Role::Title:
      // figure out which argument is being queried and return title
      switch (argument)
      {
      case InputData: return tr("Input:");
      case OutputFile: return tr("Output:");
      case LogFile: return tr("Log File:");
      case FilterSize: return tr("Filter Size:");
      default: return QVariant();
      }
   case Role::WhatsThis:
      // figure out which argument is being queried and return "What's This?" text
      switch (argument)
      {
      case InputData: return tr("Input correlation matrix that will be used to output correlation"
                                " above a threshold determined by Random Matrix Theory.");
      case OutputFile: return tr("Output text file that will hold listing of all gene correlations"
                                 " above a certain threshold.");
      case LogFile: return tr("Output text file that logs all chi trials.");
      case FilterSize: return tr("Sample size for low pass filter of chi results.");
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
      case OutputFile: return tr("Text File %1").arg("(*.txt)");
      case LogFile: return tr("Text File %1").arg("(*.txt)");
      default: return QVariant();
      }
   case Role::DefaultValue:
      // see if this is filter size and return default else return nothing
      switch (argument)
      {
      case FilterSize: return 10;
      default: return QVariant();
      }
   case Role::Minimum:
      // see if this is filter size and return minimum else return nothing
      switch (argument)
      {
      case FilterSize: return 5;
      default: return QVariant();
      }
   case Role::Maximum:
      // see if this is filter size and return maximum else return nothing
      switch (argument)
      {
      case FilterSize: return 200;
      default: return QVariant();
      }
   default:
      return QVariant();
   }
}






void RMT::setArgument(int argument, QVariant value)
{
   // set filter size if this is correct argument
   if ( argument == FilterSize )
   {
      _filterSize = value.toInt();
   }
}






void RMT::setArgument(int argument, QFile* file)
{
   // figure out which argument is being queried and set file pointer
   switch (argument)
   {
   case OutputFile:
      _output = file;
      break;
   case LogFile:
      _logfile = file;
      break;
   }
}






void RMT::setArgument(int argument, EAbstractData* data)
{
   // set input argument if this is correct argument
   if ( argument == InputData )
   {
      _input = dynamic_cast<CorrelationMatrix*>(data);
   }
}






bool RMT::initialize()
{
   // make sure input and output were set properly
   if ( !_input || !_output || !_logfile )
   {
      // report argument error and fail
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Argument Error"));
      e.setDetails(QObject::tr("Did not get valid input, output, or logfile arguments."));
      throw e;
   }

   // nothing to pre-allocate
   return false;
}






void RMT::runSerial()
{
   // determine threshold using RMT and make sure interruption is not requested
   float threshold {determineThreshold()};
   if ( isInterruptionRequested() )
   {
      return;
   }

   // initialize text stream for output and get gene names
   QTextStream stream(_output);
   const EMetadata::List* names {_input->geneNames().toArray()};

   // initialize pair iterator, xy gene indexes, and last percent complete
   CorrelationMatrix::Pair pair(_input);
   GenePair::Vector vector;
   int lastPercent {66};

   // iterate through all gene pairs until end is reached
   while ( vector.geneX() < _input->geneSize() )
   {
      // make sure interruption is not requested
      if ( isInterruptionRequested() )
      {
         return;
      }

      // read gene pair and check if it meets threshold
      pair.read(vector);
      if ( pair.at(0,0) >= threshold )
      {
         // write correlation to output
         stream << *(names->at(vector.geneX())->toString()) << "\t"
                << *(names->at(vector.geneY())->toString()) << "\tco\t" << pair.at(0,0) << "\n";
      }

      // increment to next gene pair
      ++vector;

      // determine new percentage complete and check if new
      int newPercent {66 + 34*vector.geneX()/_input->geneSize()};
      if ( newPercent != lastPercent )
      {
         // update to new percentage and emit progressed signal
         lastPercent = newPercent;
         emit progressed(lastPercent);
      }
   }
}






float RMT::determineThreshold()
{
   // generate list of maximum threshold for each gene and initialize log text stream
   generateGeneThresholds();
   QTextStream stream(_logfile);

   // initialize chi, threshold, and history of both
   float chi {0.0};
   float threshold {_initialThreshold};
   QQueue<float> previousChi;

   // initialize last percent, steps, and total steps
   int lastPercent {33};
   int steps {0};
   int totalSteps = (_initialThreshold - _thresholdMinimum)/_thresholdStep;

   // continue while chi is less than 200
   while ( chi < 100.0 )
   {
      // make sure interruption is not requested
      if ( isInterruptionRequested() )
      {
         return 0.0;
      }

      // add new raw chi value to queue
      int size;
      previousChi.enqueue(determineChi(threshold,&size));

      // check if queue has enough samples
      if ( previousChi.size() > _filterSize )
      {
         // remove last chi in queue and get sum of all chi samples
         previousChi.dequeue();
         float sum {0};
         for (const auto& chi : previousChi)
         {
            sum += chi;
         }

         // divide the sum by total number of samples and set that to filtered chi
         chi = sum/(float)_filterSize;
      }

      // output to log file
      stream << threshold << "\t" << previousChi.back() << "\t" << chi << "\t" << size << "\n";

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
      int newPercent {33 + 33*steps/totalSteps};
      if ( newPercent != lastPercent )
      {
         // update new percentage and emit progressed signal
         lastPercent = newPercent;
         emit progressed(lastPercent);
      }

      // increment steps taken
      ++steps;
   }

   // return threshold where chi was first above 100
   return threshold - _thresholdStep;
}






float RMT::determineChi(float threshold, int* size)
{
   // initialize size and generate pruned matrix based off threshold
   QVector<double> pruneMatrix {generatePruneMatrix(threshold,size)};

   // check and make sure matrix is not empty
   if ( *size > 0 )
   {
      // make sure interruption is not requested
      if ( isInterruptionRequested() )
      {
         return 0.0;
      }

      // generate eigen vector of pruned matrix and make sure interruption is not requested
      QVector<double> eigens {generateMatrixEigens(&pruneMatrix,*size)};
      if ( isInterruptionRequested() )
      {
         return 0.0;
      }

      // generate chi from eigen vector and return it
      return getChiSquare(&eigens);
   }

   // if no real chi was found return zero
   return 0.0;
}






void RMT::generateGeneThresholds()
{
   // resize thresholds matrix and initialize all to the minimum
   _geneThresholds.fill(-1.0,_input->geneSize());

   // initialize percent complete and steps
   int lastPercent {0};
   qint64 steps {0};
   qint64 totalSteps {_input->geneSize()*(_input->geneSize() - 1)/2};

   // xy gene indexes and iterator
   GenePair::Vector vector;
   CorrelationMatrix::Pair pair(_input);

   // iterate through all gene pairs
   while ( vector.geneX() < _input->geneSize() )
   {
      // make sure interruption is not requested
      if ( isInterruptionRequested() )
      {
         return;
      }

      // read in gene pair and check if it is a real number
      pair.read(vector);
      if ( !pair.isEmpty() && !isnan(pair.at(0,0)) )
      {
         // if value is greater than current max of gene x set it to new max
         if ( pair.at(0,0) > _geneThresholds.at(vector.geneX()) )
         {
            _geneThresholds[vector.geneX()] = pair.at(0,0);
         }

         // if value is greater than current max of gene y set it to new max
         if ( pair.at(0,0) > _geneThresholds.at(vector.geneY()) )
         {
            _geneThresholds[vector.geneY()] = pair.at(0,0);
         }
      }

      // increment to next gene pair, steps, and compute new percent complete
      ++vector;
      ++steps;
      qint64 newPercent {33*steps/totalSteps};

      // check to see if percent complete changed
      if ( newPercent != lastPercent )
      {
         // update percent complete and emit progressed signal
         lastPercent = newPercent;
         emit progressed(lastPercent);
      }
   }
}






QVector<double> RMT::generatePruneMatrix(float threshold, int* size)
{
   // make random number generator with normal distribution from -1 to 1
   default_random_engine generator;
   normal_distribution<float> distribution(0.0,0.3);

   // generate vector of gene indexes that have max threshold above given threshold
   QVector<int> genes;
   for (int i = 0; i < _input->geneSize() ;++i)
   {
      if ( _geneThresholds.at(i) >= threshold )
      {
         genes.push_back(i);
      }
   }

   // allocate new pruned matrix with gene vector size and initialize iterator
   QVector<double> pruneMatrix(genes.size()*genes.size());
   CorrelationMatrix::Pair pair(_input);

   // populate the pruned matrix with gene correlation values
   for (int i = 0; i < genes.size() ;++i)
   {
      for (int j = 0; j < genes.size() ;++j)
      {
         // make sure interruption is not requested
         if ( isInterruptionRequested() )
         {
            return QVector<double>();
         }

         // get indexes for both genes
         int g1 {genes[i]};
         int g2 {genes[j]};

         // if genes have same index set ij matrix value to 1
         if ( g1 == g2 )
         {
            pruneMatrix[i*genes.size() + j] = 1.0;
         }

         // else get correlation value
         else
         {
            // make sure first gene index is bigger and read gene pair
            if ( g2 > g1 )
            {
               swap(g1,g2);
            }
            pair.read({g1,g2});

            // if the correlation is a real number get value
            if ( !pair.isEmpty() && !isnan(pair.at(0,0)) && !isinf(pair.at(0,0)) )
            {
               pruneMatrix[i*genes.size() + j] = pair.at(0,0);
            }

            // else set value to random value
            else
            {
               pruneMatrix[i*genes.size() + j] = distribution(generator);
            }
         }
      }
   }

   // set size to size of pruned matrix and return pointer
   *size = genes.size();
   return pruneMatrix;
}






QVector<double> RMT::generateMatrixEigens(QVector<double>* pruneMatrix, int size)
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






float RMT::getChiSquare(QVector<double>* eigens)
{
   // initialize chi value and degenerate eigen vector
   float chi {0.0};
   degenerate(eigens);

   // make sure new size of eigen vector is minimum size
   if ( eigens->size() >= _minEigenVectorSize )
   {
      // initialize chi test count and iterate through all paces
      int chiTestCount {0};
      for (int pace = _minUnfoldingPace; pace < _maxUnfoldingPace ;++pace)
      {
         // if ??? increment chi value with paced chi square value and increment test count
         if ( (eigens->size()/pace) >= 5 )
         {
            chi += getPaceChiSquare(*eigens,pace);
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






float RMT::getPaceChiSquare(const QVector<double>& eigens, int pace)
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






void RMT::degenerate(QVector<double>* eigens)
{
   // iterate through all eigen values
   for (auto& eigen : *eigens)
   {
      // if eigen value is less than minimum set it to zero
      if ( eigen < _minimumEigenValue )
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
