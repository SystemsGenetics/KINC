#include "similarity.h"
#include "datafactory.h"
#include "genepair_gmm.h"
#include "genepair_kmeans.h"
#include "genepair_pearson.h"
#include "genepair_spearman.h"



using namespace std;






const char* Similarity::GMM {QT_TR_NOOP("gmm")};
const char* Similarity::KMeans {QT_TR_NOOP("kmeans")};
const char* Similarity::Pearson {QT_TR_NOOP("pearson")};
const char* Similarity::Spearman {QT_TR_NOOP("spearman")};
const char* Similarity::BIC {QT_TR_NOOP("BIC")};
const char* Similarity::ICL {QT_TR_NOOP("ICL")};






EAbstractAnalytic::ArgumentType Similarity::getArgumentData(int argument)
{
   // use type declaration
   using Type = EAbstractAnalytic::ArgumentType;

   // figure out which argument is being queried and return its type
   switch (argument)
   {
   case InputData: return Type::DataIn;
   case ClusterData: return Type::DataOut;
   case CorrelationData: return Type::DataOut;
   case ClusterArg: return Type::Combo;
   case CorrelationArg: return Type::Combo;
   case MinExpression: return Type::Double;
   case MinSamples: return Type::Integer;
   case MinClusters: return Type::Integer;
   case MaxClusters: return Type::Integer;
   case CriterionArg: return Type::Combo;
   case RemovePreOutliers: return Type::Bool;
   case RemovePostOutliers: return Type::Bool;
   case MinCorrelation: return Type::Double;
   case MaxCorrelation: return Type::Double;
   default: return Type::Bool;
   }
}






QVariant Similarity::getArgumentData(int argument, Role role)
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
      case ClusterData: return QString("clus");
      case CorrelationData: return QString("corr");
      case ClusterArg: return QString("clusmethod");
      case CorrelationArg: return QString("corrmethod");
      case MinExpression: return QString("minexpr");
      case MinSamples: return QString("minsamp");
      case MinClusters: return QString("minclus");
      case MaxClusters: return QString("maxclus");
      case CriterionArg: return QString("crit");
      case RemovePreOutliers: return QString("preout");
      case RemovePostOutliers: return QString("postout");
      case MinCorrelation: return QString("mincorr");
      case MaxCorrelation: return QString("maxcorr");
      default: return QVariant();
      }
   case Role::Title:
      // figure out which argument is being queried and return title
      switch (argument)
      {
      case InputData: return tr("Expressiom Matrix:");
      case ClusterData: return tr("Cluster Matrix:");
      case CorrelationData: return tr("Correlation Matrix:");
      case ClusterArg: return tr("Clustering Method:");
      case CorrelationArg: return tr("Correlation Method:");
      case MinExpression: return tr("Minimum Expression:");
      case MinSamples: return tr("Minimum Sample Size:");
      case MinClusters: return tr("Minimum Clusters:");
      case MaxClusters: return tr("Maximum Clusters:");
      case CriterionArg: return tr("Criterion:");
      case RemovePreOutliers: return tr("Remove pre-clustering outliers:");
      case RemovePostOutliers: return tr("Remove post-clustering outliers:");
      case MinCorrelation: return tr("Minimum Correlation:");
      case MaxCorrelation: return tr("Maximum Correlation:");
      default: return QVariant();
      }
   case Role::WhatsThis:
      // figure out which argument is being queried and return "What's This?" text
      switch (argument)
      {
      case InputData: return tr("Input expression matrix.");
      case ClusterData: return tr("Output matrix that will contain gene pair clusters.");
      case CorrelationData: return tr("Output matrix that will contain gene pair correlations.");
      case ClusterArg: return tr("Clustering method to use for gene pairs.");
      case CorrelationArg: return tr("Correlation method to use for gene pairs.");
      case MinExpression: return tr("Minimum threshold for a sample to be included in a gene pair.");
      case MinSamples: return tr("Minimum number of shared samples for a gene pair to be processed.");
      case MinClusters: return tr("Minimum number of clusters to test.");
      case MaxClusters: return tr("Maximum number of clusters to test.");
      case CriterionArg: return tr("Criterion to select a clustering model.");
      case RemovePreOutliers: tr("Remove pre-clustering outliers.");
      case RemovePostOutliers: tr("Remove post-clustering outliers.");
      case MinCorrelation: return tr("Minimum threshold (absolute value) for a correlation to be saved.");
      case MaxCorrelation: return tr("Maximum threshold (absolute value) for a correlation to be saved.");
      default: return QVariant();
      }
   case Role::ComboValues:
      // if this is criterion argument return combo values else return nothing
      switch (argument)
      {
      case ClusterArg: return QStringList({ tr(GMM), tr(KMeans) });
      case CorrelationArg: return QStringList({ tr(Pearson), tr(Spearman) });
      case CriterionArg: return QStringList({ tr(BIC), tr(ICL) });
      default: return QStringList();
      }
   case Role::DefaultValue:
      // figure out which argument is being queried and if applicable return default value else
      // return nothing
      switch (argument)
      {
      case ClusterArg: return tr(GMM);
      case CorrelationArg: return tr(Pearson);
      case MinExpression: return -INFINITY;
      case MinSamples: return 30;
      case MinClusters: return 1;
      case MaxClusters: return 5;
      case CriterionArg: return tr(BIC);
      case RemovePreOutliers: return false;
      case RemovePostOutliers: return false;
      case MinCorrelation: return 0.5;
      case MaxCorrelation: return 1.0;
      default: return QVariant();
      }
   case Role::Minimum:
      // figure out which argument is being queried and if applicable return minimum value else
      // return nothing
      switch (argument)
      {
      case MinExpression: return -INFINITY;
      case MinSamples: return 1;
      case MinClusters: return 1;
      case MaxClusters: return 1;
      case MinCorrelation: return 0.0;
      case MaxCorrelation: return 0.0;
      default: return QVariant();
      }
   case Role::Maximum:
      // figure out which argument is being queried and if applicable return maximum value else
      // return nothing
      switch (argument)
      {
      case MinExpression: return +INFINITY;
      case MinSamples: return INT_MAX;
      case MinClusters: return GenePair::Vector::MAX_CLUSTER_SIZE;
      case MaxClusters: return GenePair::Vector::MAX_CLUSTER_SIZE;
      case MinCorrelation: return 1.0;
      case MaxCorrelation: return 1.0;
      default: return QVariant();
      }
   case Role::Decimals:
      switch (argument)
      {
      case MinCorrelation:
      case MaxCorrelation:
         return 6;
      default:
         return QVariant();
      }
   case Role::DataType:
      // figure out which argument is being queried and if applicable return data type else
      // return nothing
      switch (argument)
      {
      case InputData: return DataFactory::ExpressionMatrixType;
      case ClusterData: return DataFactory::CCMatrixType;
      case CorrelationData: return DataFactory::CorrelationMatrixType;
      default: return QVariant();
      }
   default:
      return QVariant();
   }
}






void Similarity::setArgument(int argument, QVariant value)
{
   // figure out which argument is being set and set it
   switch (argument)
   {
   case ClusterArg:
      {
         const QString option = value.toString();
         if ( option == tr(GMM) )
         {
            _clusMethod = new GenePair::GMM();
         }
         else if ( option == tr(KMeans) )
         {
            _clusMethod = new GenePair::KMeans();
         }
      }
      break;
   case CorrelationArg:
      {
         const QString option = value.toString();
         if ( option == tr(Pearson) )
         {
            _corrMethod = new GenePair::Pearson();
         }
         else if ( option == tr(Spearman) )
         {
            _corrMethod = new GenePair::Spearman();
         }
      }
      break;
   case MinExpression:
      _minExpression = value.toDouble();
      break;
   case MinSamples:
      _minSamples = value.toInt();
      break;
   case MinClusters:
      _minClusters = value.toInt();
      break;
   case MaxClusters:
      _maxClusters = value.toInt();
      break;
   case CriterionArg:
      {
         const QString option = value.toString();
         if ( option == tr(BIC) )
         {
            _criterion = GenePair::Criterion::BIC;
         }
         else if ( option == tr(ICL) )
         {
            _criterion = GenePair::Criterion::ICL;
         }
      }
      break;
   case RemovePreOutliers:
      _removePreOutliers = value.toBool();
      break;
   case RemovePostOutliers:
      _removePostOutliers = value.toBool();
      break;
   case MinCorrelation:
      _minCorrelation = value.toDouble();
      break;
   case MaxCorrelation:
      _maxCorrelation = value.toDouble();
      break;
   }
}






void Similarity::setArgument(int argument, EAbstractData *data)
{
   // figure out which argument is having its data set and if applicable set it
   switch (argument)
   {
   case InputData:
      _input = dynamic_cast<ExpressionMatrix*>(data);
      break;
   case ClusterData:
      _clusMatrix = dynamic_cast<CCMatrix*>(data);
      break;
   case CorrelationData:
      _corrMatrix = dynamic_cast<CorrelationMatrix*>(data);
      break;
   }
}






bool Similarity::initialize()
{
   // make sure there is valid input and output
   if ( !_input || !_clusMatrix || !_corrMatrix )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Argument Error"));
      e.setDetails(QObject::tr("Did not get valid input and/or output arguments."));
      throw e;
   }

   // make sure minimum is a legal value
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

   // initialize cluster matrix
   _clusMethod->initialize(_input, _clusMatrix);

   // initialize correlation matrix
   _corrMethod->initialize(_input, _corrMatrix);

   // return pre-allocation argument
   return false;
}






void Similarity::savePair(GenePair::Vector vector, qint8 K, const qint8 *labels, int N, const float *correlations)
{
   // save clusters whose correlations are within correlations
   if ( K > 1 )
   {
      CCMatrix::Pair clusPair(_clusMatrix);

      for ( qint8 k = 0; k < K; ++k )
      {
         float corr = correlations[k];

         if ( !isnan(corr) && _minCorrelation <= abs(corr) && abs(corr) <= _maxCorrelation )
         {
            clusPair.addCluster();

            for ( int i = 0; i < N; ++i )
            {
               clusPair.at(clusPair.clusterSize() - 1, i) = (labels[i] >= 0)
                  ? (k == labels[i])
                  : -labels[i];
            }
         }
      }

      if ( clusPair.clusterSize() > 0 )
      {
         clusPair.write(vector);
      }
   }

   // save correlations that are within thresholds
   if ( K > 0 )
   {
      CorrelationMatrix::Pair corrPair(_corrMatrix);

      for ( qint8 k = 0; k < K; ++k )
      {
         float corr = correlations[k];

         if ( !isnan(corr) && _minCorrelation <= abs(corr) && abs(corr) <= _maxCorrelation )
         {
            corrPair.addCluster();
            corrPair.at(corrPair.clusterSize() - 1, 0) = corr;
         }
      }

      if ( corrPair.clusterSize() > 0 )
      {
         corrPair.write(vector);
      }
   }
}






void Similarity::runSerial()
{
   // initialize percent complete and steps
   int lastPercent {0};
   qint64 steps {0};
   qint64 totalSteps {_input->getGeneSize()*(_input->getGeneSize() - 1)/2};

   // initialize gene pair index
   GenePair::Vector vector;

   // iterate through all gene pairs
   while ( vector.geneX() < _input->getGeneSize() )
   {
      // make sure interruption is not requested
      if ( isInterruptionRequested() )
      {
         return;
      }

      // compute clustering model
      _clusMethod->compute(
         vector,
         _minSamples,
         _minExpression,
         _minClusters,
         _maxClusters,
         _criterion,
         _removePreOutliers,
         _removePostOutliers
      );

      qint8 K {_clusMethod->clusterSize()};
      const QVector<qint8>& labels {_clusMethod->labels()};

      // compute correlation
      QVector<float> correlations(K);

      for ( qint8 k = 0; k < K; ++k )
      {
         correlations[k] = _corrMethod->compute(_clusMethod->dataMatrix(), labels, k, _minSamples);
      }

      // save gene pair data
      savePair(vector, K, labels.data(), labels.size(), correlations.data());

      // increment to next pair
      ++vector;

      // increment steps and calculate percent complete
      ++steps;
      qint64 newPercent {100*steps/totalSteps};

      // emit progressed signal when percent changes
      if ( newPercent != lastPercent )
      {
         lastPercent = newPercent;
         emit progressed(lastPercent);
      }
   }
}
