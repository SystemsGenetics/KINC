#include <ace/core/metadata.h>
#include <gsl/gsl_matrix.h>

#include "gmm.h"
#include "datafactory.h"

using namespace std;






EAbstractAnalytic::ArgumentType GMM::getArgumentData(int argument)
{
   // use type declaration
   using Type = EAbstractAnalytic::ArgumentType;

   // figure out which argument is being queried and return its type
   switch (argument)
   {
   case InputData: return Type::DataIn;
   case OutputData: return Type::DataOut;
   case MinSamples: return Type::Integer;
   case MinClusters: return Type::Integer;
   case MaxClusters: return Type::Integer;
   default: return Type::Bool;
   }
}






QVariant GMM::getArgumentData(int argument, Role role)
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
      case OutputData: return QString("output");
      case MinSamples: return QString("min");
      case MinClusters: return QString("minclus");
      case MaxClusters: return QString("maxclus");
      default: return QVariant();
      }
   case Role::Title:
      // figure out which argument is being queried and return title
      switch (argument)
      {
      case InputData: return tr("Input:");
      case OutputData: return tr("Output:");
      case MinSamples: return tr("Minimum Sample Size:");
      case MinClusters: return tr("Minimum Clusters:");
      case MaxClusters: return tr("Maximum Clusters:");
      default: return QVariant();
      }
   case Role::WhatsThis:
      // figure out which argument is being queried and return "What's This?" text
      switch (argument)
      {
      case InputData: return tr("Input expression matrix that will be used to compute pair-wise clusters.");
      case OutputData: return tr("Output cluster matrix that will store pair-wise clusters.");
      case MinSamples: return tr("Minimum size of samples two genes must share to perform clustering.");
      case MinClusters: return tr("Minimum number of clusters to test.");
      case MaxClusters: return tr("Maximum number of clusters to test.");
      default: return QVariant();
      }
   case Role::DefaultValue:
      // figure out which argument is being queried and if applicable return default value else
      // return nothing
      switch (argument)
      {
      case MinSamples: return 30;
      case MinClusters: return 1;
      case MaxClusters: return 5;
      default: return QVariant();
      }
   case Role::Minimum:
      // figure out which argument is being queried and if applicable return minimum value else
      // return nothing
      switch (argument)
      {
      case MinSamples: return 1;
      case MinClusters: return 1;
      case MaxClusters: return 1;
      default: return QVariant();
      }
   case Role::Maximum:
      // figure out which argument is being queried and if applicable return maximum value else
      // return nothing
      switch (argument)
      {
      case MinSamples: return INT_MAX;
      case MinClusters: return INT_MAX;
      case MaxClusters: return INT_MAX;
      default: return QVariant();
      }
   case Role::DataType:
      // figure out which argument is being queried and if applicable return data type else
      // return nothing
      switch (argument)
      {
      case InputData: return DataFactory::ExpressionMatrixType;
      case OutputData: return DataFactory::CCMatrixType;
      default: return QVariant();
      }
   default:
      return QVariant();
   }
}






void GMM::setArgument(int argument, QVariant value)
{
   // figure out which argument is being set and set it
   switch (argument)
   {
   case MinSamples:
      _minSamples = value.toInt();
      break;
   case MinClusters:
      _minClusters = value.toInt();
      break;
   case MaxClusters:
      _maxClusters = value.toInt();
      break;
   }
}






void GMM::setArgument(int argument, EAbstractData *data)
{
   // figure out which argument is having its data set and if applicable set it
   switch (argument)
   {
   case InputData:
      _input = dynamic_cast<ExpressionMatrix*>(data);
      break;
   case OutputData:
      _output = dynamic_cast<CCMatrix*>(data);
      break;
   }
}






bool GMM::initialize()
{
   // make sure there is valid input and output
   if ( !_input || !_output )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Argument Error"));
      e.setDetails(QObject::tr("Did not get valid input and/or output arguments."));
      throw e;
   }

   // make sure minimum sample size is a legal value
   if ( _minSamples < 1 )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Argument Error"));
      e.setDetails(QObject::tr("Minimum sample size must be at least 1 or greater."));
      throw e;
   }

   // make sure cluster range is valid
   if ( _maxClusters < _minClusters )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Argument Error"));
      e.setDetails(QObject::tr("Minimum clusters must be less than or equal to maximum clusters."));
      throw e;
   }

   // initialize new cc matrix and return pre-allocation argument
   _output->initialize(_input->getGeneNames(),_input->getSampleNames());
   return false;
}






float GMM::computeBIC(const GenePair::GMM& gmm, int N, int D)
{
   int p = gmm.numClusters() * (1 + D + D * D);

   return log(N) * p - 2 * gmm.logLikelihood();
}






CCMatrix::Pair GMM::computePair(const float *X, int N, int D)
{
   // run each clustering model
   QVector<GenePair::GMM> models(_maxClusters - _minClusters + 1);

   for ( int K = _minClusters; K <= _maxClusters; ++K )
   {
      auto& gmm = models[K - _minClusters];

      gmm.fit(X, N, D, K);
   }

   // select the model with the lowest criterion value
   float bestValue = 0;
   auto bestModel = models.end();

   for ( auto iter = models.begin(); iter != models.end(); ++iter )
   {
      if ( !iter->success() )
      {
         continue;
      }

      float value = computeBIC(*iter, N, D);

      if ( bestModel == models.end() || value < bestValue )
      {
         bestValue = value;
         bestModel = iter;
      }
   }

   // compute cluster labels for gene pair
   CCMatrix::Pair pair(_output);

   if ( bestModel == models.end() )
   {
      fprintf(stderr, "warning: all models failed\n");
      return pair;
   }

   pair.addCluster(bestModel->numClusters());

   for ( int i = 0; i < N; ++i )
   {
      for ( int k = 0; k < bestModel->numClusters(); ++k )
      {
         pair.at(k, i) = (k == bestModel->labels()[i]);
      }
   }

   return pair;
}






void GMM::runSerial()
{
   // initialize percent complete and steps
   int lastPercent {0};
   qint64 steps {0};
   qint64 totalSteps {_output->geneSize()*(_output->geneSize() - 1)/2};

   // initialize arrays used for k-means clustering
   float X[_input->getSampleSize() * 2];

   // initialize expression genes for input/output
   ExpressionMatrix::Gene gene1(_input);
   ExpressionMatrix::Gene gene2(_input);

   // initialize xy gene indexes
   GenePair::Vector vector;

   // increment through all gene pairs
   while ( vector.geneX() < _output->geneSize() )
   {
      // make sure interruption is not requested
      if ( isInterruptionRequested() )
      {
         return;
      }

      // initialize sample size and read in gene expressions
      int size {0};
      gene1.read(vector.geneX());
      gene2.read(vector.geneY());

      // populate a and b arrays with shared expressions of gene x and y
      for ( auto i = 0; i < _input->getSampleSize(); ++i )
      {
         if ( !std::isnan(gene1.at(i)) && !std::isnan(gene2.at(i)) )
         {
            X[size * 2 + 0] = gene1.at(i);
            X[size * 2 + 1] = gene2.at(i);
            ++size;
         }
      }

      // perform clustering only if there are enough samples
      if ( size >= _minSamples )
      {
         CCMatrix::Pair pair = computePair(X, size, 2);

         pair.write(vector);
      }

      // increment to next pair
      ++vector;

      // increment steps and calculate percent complete
      ++steps;
      qint64 newPercent {100*steps/totalSteps};

      // check to see if percent has changed
      if ( newPercent != lastPercent )
      {
         // update percent complete and emit progressed signal
         lastPercent = newPercent;
         emit progressed(lastPercent);
      }
   }
}
