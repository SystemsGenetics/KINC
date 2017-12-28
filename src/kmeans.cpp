#include <ace/core/metadata.h>
#include <gsl/gsl_cblas.h>

#include "kmeans.h"
#include "datafactory.h"


using namespace std;






KMeans::~KMeans()
{
}






EAbstractAnalytic::ArgumentType KMeans::getArgumentData(int argument)
{
   // use type declaration
   using Type = EAbstractAnalytic::ArgumentType;

   // figure out which argument is being queried and return its type
   switch (argument)
   {
   case InputData: return Type::DataIn;
   case OutputData: return Type::DataOut;
   case Minimum: return Type::Integer;
   case MinClusters: return Type::Integer;
   case MaxClusters: return Type::Integer;
   default: return Type::Bool;
   }
}






QVariant KMeans::getArgumentData(int argument, Role role)
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
      case Minimum: return QString("min");
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
      case Minimum: return tr("Minimum Sample Size:");
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
      case Minimum: return tr("Minimum size of samples two genes must share to perform clustering.");
      case MinClusters: return tr("Minimum number of clusters to test.");
      case MaxClusters: return tr("Maximum number of clusters to test.");
      default: return QVariant();
      }
   case Role::DefaultValue:
      // figure out which argument is being queried and if applicable return default value else
      // return nothing
      switch (argument)
      {
      case Minimum: return 30;
      case MinClusters: return 1;
      case MaxClusters: return 5;
      default: return QVariant();
      }
   case Role::Minimum:
      // figure out which argument is being queried and if applicable return minimum value else
      // return nothing
      switch (argument)
      {
      case Minimum: return 1;
      case MinClusters: return 1;
      case MaxClusters: return 1;
      default: return QVariant();
      }
   case Role::Maximum:
      // figure out which argument is being queried and if applicable return maximum value else
      // return nothing
      switch (argument)
      {
      case Minimum: return INT_MAX;
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






void KMeans::setArgument(int argument, QVariant value)
{
   // figure out which argument is being set and set it
   switch (argument)
   {
   case Minimum:
      _minimum = value.toInt();
      break;
   case MinClusters:
      _minClusters = value.toInt();
      break;
   case MaxClusters:
      _maxClusters = value.toInt();
      break;
   }
}






void KMeans::setArgument(int argument, EAbstractData *data)
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






bool KMeans::initialize()
{
   // make sure there is valid input and output
   if ( !_input || !_output )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Argument Error"));
      e.setDetails(QObject::tr("Did not get valid input and/or output arguments."));
      throw e;
   }

   // make sure minimum is a legal value
   if ( _minimum < 1 )
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








float computeVectorDistance(const float *a, const float *b, int n)
{
   float dist = 0;
   for ( int i = 0; i < n; ++i )
   {
      float diff = a[i] - b[i];
      dist += diff * diff;
   }

   return sqrt(dist);
}

QVector<int> KMeans::computeKmeans(const float *X, int N, int D, float *Mu, int K)
{
   QVector<int> y;
   QVector<int> y_next(N);

   while ( true )
   {
      // E step
      for ( int i = 0; i < N; ++i )
      {
         // find k that minimizes norm(x_i - mu_k)
         int min_k = -1;
         float min_dist;

         for ( int k = 0; k < K; ++k )
         {
            const float *x_i = &X[i * D];
            const float *mu_k = &Mu[k * D];

            float dist = computeVectorDistance(x_i, mu_k, D);

            if ( min_k == -1 || dist < min_dist )
            {
               min_k = k;
               min_dist = dist;
            }
         }

         y_next[i] = min_k;
      }

      // check for convergence
      if ( y == y_next )
      {
         break;
      }

      y = y_next;

      // M step
      for ( int k = 0; k < K; ++k )
      {
         // compute mu_k = mean of all x_i in cluster k
         float *mu_k = &Mu[k * D];
         int n_k = 0;

         memset(mu_k, 0, D * sizeof(float));

         for ( int i = 0; i < N; ++i )
         {
            const float *x_i = &X[i * D];

            if ( y[i] == k )
            {
               cblas_saxpy(D, 1.0f, x_i, 1, mu_k, 1);
               n_k++;
            }
         }

         cblas_sscal(D, 1.0f / n_k, mu_k, 1);
      }
   }

   return y;
}








float computeLogLikelihood(const float *X, int N, int D, const float *Mu, int K, const QVector<int>& y)
{
   // compute within-class scatter
   float S = 0;

   for ( int k = 0; k < K; ++k )
   {
      for ( int i = 0; i < N; ++i )
      {
         if ( y[i] != k )
         {
            continue;
         }

         const float *mu_k = &Mu[k * D];
         const float *x_i = &X[i * D];

         float dist = computeVectorDistance(x_i, mu_k, D);

         S += dist * dist;
      }
   }

   return -S;
}

float computeBIC(const float *X, int N, int D, const float *Mu, int K, const QVector<int>& y)
{
   int p = K * D;
   float L = computeLogLikelihood(X, N, D, Mu, K, y);

   return log(N) * p - 2 * L;
}

CCMatrix::Pair KMeans::computePair(const float *X, int N, int D)
{
   // compute clustering model
   int bestK = 0;
   QVector<int> bestLabels;
   float bestValue = 0;

   float Mu[_maxClusters * D];

   for ( int K = _minClusters; K <= _maxClusters; ++K )
   {
      // run each clustering model
      QVector<int> labels = computeKmeans(X, N, D, Mu, K);

      // select the model with lowest criterion value
      float value = computeBIC(X, N, D, Mu, K, labels);

      if ( bestK == 0 || value < bestValue ) {
         bestK = K;
         bestLabels = labels;
         bestValue = value;
      }
   }

   // compute cluster labels for gene pair
   CCMatrix::Pair pair(_output);
   pair.addCluster(bestK);

   for ( int i = 0; i < bestLabels.size(); ++i )
   {
      for ( int k = 0; k < bestK; ++k )
      {
         pair.at(k, i) = (k == bestLabels[i]);
      }
   }

   return pair;
}








void KMeans::runSerial()
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
      if ( size >= _minimum )
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
