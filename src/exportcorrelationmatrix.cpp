#include "exportcorrelationmatrix.h"
#include "datafactory.h"



EAbstractAnalytic::ArgumentType ExportCorrelationMatrix::getArgumentData(int argument)
{
   // use type declaration
   using Type = EAbstractAnalytic::ArgumentType;

   // figure out which argument is being queried and return its type
   switch (argument)
   {
   case ClusterData: return Type::DataIn;
   case CorrelationData: return Type::DataIn;
   case OutputFile: return Type::FileOut;
   default: return Type::Bool;
   }
}






QVariant ExportCorrelationMatrix::getArgumentData(int argument, Role role)
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
      case ClusterData: return QString("clus");
      case CorrelationData: return QString("corr");
      case OutputFile: return QString("output");
      default: return QString();
      }
   case Role::Title:
      // figure out which argument is being queried and return title
      switch (argument)
      {
      case ClusterData: return tr("Cluster Matrix:");
      case CorrelationData: return tr("Correlation Matrix:");
      case OutputFile: return tr("Output File:");
      default: return QString();
      }
   case Role::WhatsThis:
      // figure out which argument is being queried and return "What's This?" text
      switch (argument)
      {
      case ClusterData: return tr("Input cluster matrix containing cluster composition data.");
      case CorrelationData: return tr("Input correlation matrix containing correlation data.");
      case OutputFile: return tr("Raw output text file that will contain gene pair correlation data.");
      default: return QString();
      }
   case Role::DataType:
      // if this is input data argument return data type else return nothing
      switch (argument)
      {
      case ClusterData: return DataFactory::CCMatrixType;
      case CorrelationData: return DataFactory::CorrelationMatrixType;
      default: return QString();
      }
   case Role::FileFilters:
      // if this is output file argument return file filter else return nothing
      switch (argument)
      {
      case OutputFile: return tr("Raw text file %1").arg("(*.txt)");
      default: return QString();
      }
   default:
      return QVariant();
   }
}






void ExportCorrelationMatrix::setArgument(int argument, EAbstractData* data)
{
   // if argument is input data set it
   if ( argument == ClusterData )
   {
      _ccMatrix = dynamic_cast<CCMatrix*>(data);
   }
   else if ( argument == CorrelationData )
   {
      _cMatrix = dynamic_cast<CorrelationMatrix*>(data);
   }
}






void ExportCorrelationMatrix::setArgument(int argument, QFile* file)
{
   // if argument is output file set it
   if ( argument == OutputFile )
   {
      _output = file;
   }
}






bool ExportCorrelationMatrix::initialize()
{
   // make sure input and output arguments were set
   if ( !_ccMatrix || !_cMatrix || !_output )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Argument Error"));
      e.setDetails(QObject::tr("Did not get valid input and/or output arguments."));
      throw e;
   }

   // do not have data pre-allocate
   return false;
}






void ExportCorrelationMatrix::runSerial()
{
   // initialize percent complete and steps
   int lastPercent {0};
   qint64 steps {0};
   qint64 totalSteps {_ccMatrix->geneSize()*(_ccMatrix->geneSize() - 1)/2};

   // initialize xy gene indexes
   GenePair::Vector vector;
   int cluster = 0;
   CCMatrix::Pair ccPair(_ccMatrix);
   CorrelationMatrix::Pair cPair(_cMatrix);

   // create text stream to output file and write until end reached
   QTextStream stream(_output);
   stream.setRealNumberPrecision(6);

   // increment through all gene pairs
   while ( vector.geneX() < _ccMatrix->geneSize() )
   {
      // make sure interruption is not requested
      if ( isInterruptionRequested() )
      {
         return;
      }

      // read cluster and correlation data for gene pair
      if ( cluster == 0 )
      {
         ccPair.read(vector);
         cPair.read(vector);
      }

      // write gene pair data to output file
      if ( cPair.clusterSize() > 0 )
      {
         int numClean = 0;
         int numMissing = 0;
         int numThreshold = 0;
         int numPreOutliers = 0;
         int numPostOutliers = 0;
         QString sampleMask(_ccMatrix->sampleSize(), '0');

         if ( ccPair.clusterSize() > 1 )
         {
            for ( int i = 0; i < _ccMatrix->sampleSize(); i++ )
            {
               sampleMask[i] = '0' + ccPair.at(cluster, i);
            }
         }

         stream
            << vector.geneX()
            << "\t" << vector.geneY()
            << "\t" << cluster
            << "\t" << cPair.clusterSize()
            << "\t" << numClean
            << "\t" << numMissing
            << "\t" << numThreshold
            << "\t" << numPreOutliers
            << "\t" << numPostOutliers
            << "\t" << cPair.at(cluster, 0)
            << "\t" << sampleMask
            << "\n";
      }

      // increment to next cluster
      cluster = (cluster + 1) % cPair.clusterSize();

      if ( cluster == 0 )
      {
         ++vector;
         ++steps;
      }

      // increment steps and calculate percent complete
      qint64 newPercent {100*steps/totalSteps};

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
}
