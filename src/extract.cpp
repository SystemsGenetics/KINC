#include "extract.h"
#include "datafactory.h"



using namespace std;






EAbstractAnalytic::ArgumentType Extract::getArgumentData(int argument)
{
   // use type declaration
   using Type = EAbstractAnalytic::ArgumentType;

   // figure out which argument is being queried and return its type
   switch (argument)
   {
   case ExpressionData: return Type::DataIn;
   case ClusterData: return Type::DataIn;
   case CorrelationData: return Type::DataIn;
   case OutputFile: return Type::FileOut;
   case GraphMLFile: return Type::FileOut;
   case MinCorrelation: return Type::Double;
   case MaxCorrelation: return Type::Double;
   default: return Type::Bool;
   }
}






QVariant Extract::getArgumentData(int argument, Role role)
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
      case ExpressionData: return QString("expr");
      case ClusterData: return QString("clus");
      case CorrelationData: return QString("corr");
      case OutputFile: return QString("output");
      case GraphMLFile: return QString("graphml");
      case MinCorrelation: return QString("mincorr");
      case MaxCorrelation: return QString("maxcorr");
      default: return QString();
      }
   case Role::Title:
      // figure out which argument is being queried and return title
      switch (argument)
      {
      case ExpressionData: return tr("Expression Matrix:");
      case ClusterData: return tr("Cluster Matrix:");
      case CorrelationData: return tr("Correlation Matrix:");
      case OutputFile: return tr("Output File:");
      case GraphMLFile: return tr("GraphML File:");
      case MinCorrelation: return tr("Minimum Correlation:");
      case MaxCorrelation: return tr("Maximum Correlation:");
      default: return QString();
      }
   case Role::WhatsThis:
      // figure out which argument is being queried and return "What's This?" text
      switch (argument)
      {
      case ExpressionData: return tr("Input expression matrix containing gene expression data.");
      case ClusterData: return tr("Input cluster matrix containing cluster composition data.");
      case CorrelationData: return tr("Input correlation matrix containing correlation data.");
      case OutputFile: return tr("Raw output text file that will contain network edges.");
      case GraphMLFile: return tr("Raw output text file that will contain network in GraphML format.");
      case MinCorrelation: return tr("Minimum (absolute) correlation threshold for gene pairs.");
      case MaxCorrelation: return tr("Maximum (absolute) correlation threshold for gene pairs.");
      default: return QString();
      }
   case Role::DefaultValue:
      // figure out which argument is being queried and if applicable return default value else
      // return nothing
      switch (argument)
      {
      case MinCorrelation: return 0.85;
      case MaxCorrelation: return 1.00;
      default: return QVariant();
      }
   case Role::Minimum:
      // figure out which argument is being queried and if applicable return minimum value else
      // return nothing
      switch (argument)
      {
      case MinCorrelation: return 0.0;
      case MaxCorrelation: return 0.0;
      default: return QVariant();
      }
   case Role::Maximum:
      // figure out which argument is being queried and if applicable return maximum value else
      // return nothing
      switch (argument)
      {
      case MinCorrelation: return 1.0;
      case MaxCorrelation: return 1.0;
      default: return QVariant();
      }
   case Role::DataType:
      // if this is input data argument return data type else return nothing
      switch (argument)
      {
      case ExpressionData: return DataFactory::ExpressionMatrixType;
      case ClusterData: return DataFactory::CCMatrixType;
      case CorrelationData: return DataFactory::CorrelationMatrixType;
      default: return QString();
      }
   case Role::FileFilters:
      // if this is output file argument return file filter else return nothing
      switch (argument)
      {
      case OutputFile: return tr("Raw text file %1").arg("(*.txt)");
      case GraphMLFile: return tr("GraphML file %1").arg("(*.graphml)");
      default: return QString();
      }
   default:
      return QVariant();
   }
}






void Extract::setArgument(int argument, EAbstractData* data)
{
   // if argument is input data set it
   if ( argument == ExpressionData )
   {
      _eMatrix = dynamic_cast<ExpressionMatrix*>(data);
   }
   else if ( argument == ClusterData )
   {
      _ccm = dynamic_cast<CCMatrix*>(data);
   }
   else if ( argument == CorrelationData )
   {
      _cmx = dynamic_cast<CorrelationMatrix*>(data);
   }
}






void Extract::setArgument(int argument, QFile* file)
{
   // figure out which argument is being set and set it
   if ( argument == OutputFile )
   {
      _output = file;
   }
   else if ( argument == GraphMLFile )
   {
      _graphml = file;
   }
}






void Extract::setArgument(int argument, QVariant value)
{
   // figure out which argument is being set and set it
   switch (argument)
   {
   case MinCorrelation:
      _minCorrelation = value.toFloat();
      break;
   case MaxCorrelation:
      _maxCorrelation = value.toFloat();
      break;
   }
}






bool Extract::initialize()
{
   // make sure input and output arguments were set
   if ( !_eMatrix || !_ccm || !_cmx || !_output )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Argument Error"));
      e.setDetails(QObject::tr("Did not get valid input and/or output arguments."));
      throw e;
   }

   // do not have data pre-allocate
   return false;
}






void Extract::runSerial()
{
   // initialize percent complete and steps
   int lastPercent {0};
   qint64 steps {0};
   qint64 totalSteps {_cmx->size() * 2};

   // initialize pair iterators
   CorrelationMatrix::Pair cmxPair(_cmx);
   CCMatrix::Pair ccmPair(_ccm);

   // get gene names
   const QList<EMetadata *> *geneNames = _cmx->geneNames().toArray();

   // initialize workspace
   QString sampleMask(_ccm->sampleSize(), '0');

   // create text stream to output file and write until end reached
   QTextStream stream(_output);
   stream.setRealNumberPrecision(6);

   // write header to file
   stream
      << "Source"
      << "\t" << "Target"
      << "\t" << "sc"
      << "\t" << "Interaction"
      << "\t" << "Cluster"
      << "\t" << "Num_Clusters"
      << "\t" << "Cluster_Samples"
      << "\t" << "Missing_Samples"
      << "\t" << "Cluster_Outliers"
      << "\t" << "Pair_Outliers"
      << "\t" << "Too_Low"
      << "\t" << "Samples"
      << "\n";

   // increment through all gene pairs
   while ( cmxPair.hasNext() )
   {
      // make sure interruption is not requested
      if ( isInterruptionRequested() )
      {
         return;
      }

      // read next gene pair
      cmxPair.readNext();

      if ( cmxPair.clusterSize() > 1 )
      {
         ccmPair.read(cmxPair.index());
      }

      // write gene pair data to output file
      for ( int cluster = 0; cluster < cmxPair.clusterSize(); cluster++ )
      {
         QString& source(*geneNames->at(cmxPair.index().getX())->toString());
         QString& target(*geneNames->at(cmxPair.index().getY())->toString());
         float correlation = cmxPair.at(cluster, 0);
         QString interaction("co");
         int numSamples = 0;
         int numMissing = 0;
         int numPostOutliers = 0;
         int numPreOutliers = 0;
         int numThreshold = 0;

         // exclude cluster if correlation is not within thresholds
         if ( fabs(correlation) < _minCorrelation || _maxCorrelation < fabs(correlation) )
         {
            continue;
         }

         // if there are multiple clusters then use cluster data
         if ( cmxPair.clusterSize() > 1 )
         {
            // compute summary statistics
            for ( int i = 0; i < _ccm->sampleSize(); i++ )
            {
               switch ( ccmPair.at(cluster, i) )
               {
               case 1:
                  numSamples++;
                  break;
               case 6:
                  numThreshold++;
                  break;
               case 7:
                  numPreOutliers++;
                  break;
               case 8:
                  numPostOutliers++;
                  break;
               case 9:
                  numMissing++;
                  break;
               }
            }

            // write sample mask to string
            for ( int i = 0; i < _ccm->sampleSize(); i++ )
            {
               sampleMask[i] = '0' + ccmPair.at(cluster, i);
            }
         }

         // otherwise use expression data
         else
         {
            // read in gene expressions
            ExpressionMatrix::Gene gene1(_eMatrix);
            ExpressionMatrix::Gene gene2(_eMatrix);

            gene1.read(cmxPair.index().getX());
            gene2.read(cmxPair.index().getY());

            // determine sample mask from expression data
            for ( int i = 0; i < _eMatrix->getSampleSize(); ++i )
            {
               if ( isnan(gene1.at(i)) || isnan(gene2.at(i)) )
               {
                  sampleMask[i] = '9';
               }
               else
               {
                  sampleMask[i] = '1';
               }
            }
         }

         // write cluster to output file
         stream
            << source
            << "\t" << target
            << "\t" << correlation
            << "\t" << interaction
            << "\t" << cluster
            << "\t" << cmxPair.clusterSize()
            << "\t" << numSamples
            << "\t" << numMissing
            << "\t" << numPostOutliers
            << "\t" << numPreOutliers
            << "\t" << numThreshold
            << "\t" << sampleMask
            << "\n";
      }

      // increment steps and calculate percent complete
      ++steps;
      qint64 newPercent {50*steps/totalSteps};

      // check to see if percent has changed
      if ( newPercent != lastPercent )
      {
         // update percent complete and emit progressed signal
         lastPercent = newPercent;
         emit progressed(lastPercent);
      }
   }

   // make sure writing output file worked
   if ( stream.status() != QTextStream::Ok )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Qt Text Stream encountered an unknown error."));
      throw e;
   }

   // reset gene pair iterator
   cmxPair.reset();

   // create text stream to graphml file and write until end reached
   stream.setDevice(_graphml);

   // write header to file
   stream
      << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
      << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"\n"
      << "    xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
      << "    xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n"
      << "  <graph id=\"G\" edgedefault=\"undirected\">\n";

   // write each node to file
   for ( int i = 0; i < _cmx->geneSize(); i++ )
   {
      QString& id = *geneNames->at(i)->toString();

      stream << "    <node id=\"" << id << "\"/>\n";
   }

   // increment through all gene pairs
   while ( cmxPair.hasNext() )
   {
      // make sure interruption is not requested
      if ( isInterruptionRequested() )
      {
         return;
      }

      // read next gene pair
      cmxPair.readNext();

      if ( cmxPair.clusterSize() > 1 )
      {
         ccmPair.read(cmxPair.index());
      }

      // write gene pair edges to file
      for ( int cluster = 0; cluster < cmxPair.clusterSize(); cluster++ )
      {
         QString& source(*geneNames->at(cmxPair.index().getX())->toString());
         QString& target(*geneNames->at(cmxPair.index().getY())->toString());
         float correlation = cmxPair.at(cluster, 0);

         // exclude edge if correlation is not within thresholds
         if ( fabs(correlation) < _minCorrelation || _maxCorrelation < fabs(correlation) )
         {
            continue;
         }

         // if there are multiple clusters then use cluster data
         if ( cmxPair.clusterSize() > 1 )
         {
            // write sample mask to string
            for ( int i = 0; i < _ccm->sampleSize(); i++ )
            {
               sampleMask[i] = '0' + ccmPair.at(cluster, i);
            }
         }

         // otherwise use expression data
         else
         {
            // read in gene expressions
            ExpressionMatrix::Gene gene1(_eMatrix);
            ExpressionMatrix::Gene gene2(_eMatrix);

            gene1.read(cmxPair.index().getX());
            gene2.read(cmxPair.index().getY());

            // determine sample mask from expression data
            for ( int i = 0; i < _eMatrix->getSampleSize(); ++i )
            {
               if ( isnan(gene1.at(i)) || isnan(gene2.at(i)) )
               {
                  sampleMask[i] = '9';
               }
               else
               {
                  sampleMask[i] = '1';
               }
            }
         }

         // write edge to file
         stream
            << "    <edge"
            << " source=\"" << source << "\""
            << " target=\"" << target << "\""
            << " samples=\"" << sampleMask << "\""
            << "/>\n";
      }

      // increment steps and calculate percent complete
      ++steps;
      qint64 newPercent {50 + 50*steps/totalSteps};

      // check to see if percent has changed
      if ( newPercent != lastPercent )
      {
         // update percent complete and emit progressed signal
         lastPercent = newPercent;
         emit progressed(lastPercent);
      }
   }

   // write footer to file
   stream
      << "  </graph>\n"
      << "</graphml>\n";

   // make sure writing graphml file worked
   if ( stream.status() != QTextStream::Ok )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Qt Text Stream encountered an unknown error."));
      throw e;
   }
}
