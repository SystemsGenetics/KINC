#include "importcorrelationmatrix.h"
#include "datafactory.h"
#include "pairwise_index.h"



EAbstractAnalytic::ArgumentType ImportCorrelationMatrix::getArgumentData(int argument)
{
   // use type declaration
   using Type = EAbstractAnalytic::ArgumentType;

   // figure out which argument is being queried and return its type
   switch (argument)
   {
   case InputFile: return Type::FileIn;
   case ClusterData: return Type::DataOut;
   case CorrelationData: return Type::DataOut;
   case GeneSize: return Type::Integer;
   case SampleSize: return Type::Integer;
   case CorrelationName: return Type::String;
   default: return Type::Bool;
   }
}






QVariant ImportCorrelationMatrix::getArgumentData(int argument, Role role)
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
      case InputFile: return QString("input");
      case ClusterData: return QString("clus");
      case CorrelationData: return QString("corr");
      case GeneSize: return QString("genes");
      case SampleSize: return QString("samples");
      case CorrelationName: return QString("corrname");
      default: return QString();
      }
   case Role::Title:
      // figure out which argument is being queried and return title
      switch (argument)
      {
      case InputFile: return tr("Input File:");
      case ClusterData: return tr("Cluster Matrix:");
      case CorrelationData: return tr("Correlation Matrix:");
      case GeneSize: return tr("Gene Size:");
      case SampleSize: return tr("Sample Size:");
      case CorrelationName: return tr("Correlation Name:");
      default: return QString();
      }
   case Role::WhatsThis:
      // figure out which argument is being queried and return "What's This?" text
      switch (argument)
      {
      case InputFile: return tr("Raw input text file containing gene pair correlation data.");
      case ClusterData: return tr("Output cluster matrix that will contain cluster composition data.");
      case CorrelationData: return tr("Output correlation matrix that will contain correlation data.");
      case GeneSize: return tr("Number of genes.");
      case SampleSize: return tr("Number of samples.");
      case CorrelationName: return tr("Name of correlation method.");
      default: return QString();
      }
   case Role::Minimum:
      // if this is sample size argument return minimum else return nothing
      switch (argument)
      {
      case GeneSize: return 1;
      case SampleSize: return 1;
      default: return QVariant();
      }
   case Role::Maximum:
      // if this is sample size argument return maximum else return nothing
      switch (argument)
      {
      case GeneSize: return INT_MAX;
      case SampleSize: return INT_MAX;
      default: return QVariant();
      }
   case Role::FileFilters:
      // if this is input file argument return file filter else return nothing
      switch (argument)
      {
      case InputFile: return tr("Raw text file %1").arg("(*.txt)");
      default: return QString();
      }
   case Role::DataType:
      // if this is output data argument return data type else return nothing
      switch (argument)
      {
      case ClusterData: return DataFactory::CCMatrixType;
      case CorrelationData: return DataFactory::CorrelationMatrixType;
      default: return QString();
      }
   default:
      return QVariant();
   }
}






void ImportCorrelationMatrix::setArgument(int argument, QVariant value)
{
   // figure out which argument is being set and set it
   switch (argument)
   {
   case GeneSize:
      _geneSize = value.toInt();
      break;
   case SampleSize:
      _sampleSize = value.toInt();
      break;
   case CorrelationName:
      _correlationName = value.toString();
      break;
   }
}






void ImportCorrelationMatrix::setArgument(int argument, QFile* file)
{
   // if argument is input file set it
   if ( argument == InputFile )
   {
      _input = file;
   }
}






void ImportCorrelationMatrix::setArgument(int argument, EAbstractData* data)
{
   // if argument is output data set it
   if ( argument == ClusterData )
   {
      _ccm = dynamic_cast<CCMatrix*>(data);
   }
   else if ( argument == CorrelationData )
   {
      _cmx = dynamic_cast<CorrelationMatrix*>(data);
   }
}






bool ImportCorrelationMatrix::initialize()
{
   // make sure input and output arguments were set
   if ( !_input || !_ccm || !_cmx )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Argument Error"));
      e.setDetails(QObject::tr("Did not get valid input and/or output arguments."));
      throw e;
   }

   // make sure gene size and sample size are positive
   if ( _geneSize <= 0 || _sampleSize <= 0 )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Argument Error"));
      e.setDetails(QObject::tr("Gene size and sample size must be positive."));
      throw e;
   }

   // make sure correlation name is not empty
   if ( _correlationName.isEmpty() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Argument Error"));
      e.setDetails(QObject::tr("Correlation name is required."));
      throw e;
   }

   // do not have data pre-allocate
   return false;
}






void ImportCorrelationMatrix::runSerial()
{
   // initialize variable used to track percent complete
   int lastPercent {0};

   // initialize metadata arrays
   EMetadata metaGeneNames(EMetadata::Array);
   EMetadata metaSampleNames(EMetadata::Array);
   EMetadata metaCorrelationNames(EMetadata::Array);

   // build gene name metadata
   for (int i = 0; i < _geneSize ;++i)
   {
      EMetadata* name {new EMetadata(EMetadata::String)};
      *(name->toString()) = QString::number(i);
      metaGeneNames.toArray()->append(name);
   }

   // build sample name metadata
   for (int i = 0; i < _sampleSize ;++i)
   {
      EMetadata* name {new EMetadata(EMetadata::String)};
      *(name->toString()) = QString::number(i);
      metaSampleNames.toArray()->append(name);
   }

   // build correlation name metadata
   EMetadata* name {new EMetadata(EMetadata::String)};
   *(name->toString()) = _correlationName;
   metaCorrelationNames.toArray()->append(name);

   // initialize output data
   _ccm->initialize(metaGeneNames, metaSampleNames);
   _cmx->initialize(metaGeneNames, metaCorrelationNames);

   Pairwise::Index index;
   CCMatrix::Pair ccmPair(_ccm);
   CorrelationMatrix::Pair cmxPair(_cmx);

   // create text stream from input file and read until end reached
   QTextStream stream(_input);
   while ( !stream.atEnd() )
   {
      // if interruption is requested exit immediately
      if ( isInterruptionRequested() )
      {
         return;
      }

      // read a single line from stream and split it into words using any whitespace characters as
      // delimiters
      QString line = stream.readLine();
      auto words = line.splitRef(QRegExp("\\s+"), QString::SkipEmptyParts);

      // check to see if there are any words from read line
      if ( words.size() == 11 )
      {
         int geneX = words[0].toInt();
         int geneY = words[1].toInt();
         float correlation = words[9].toFloat();
         QStringRef sampleMask = words[10];

         // make sure sample mask has correct length
         if ( sampleMask.size() != _sampleSize )
         {
            E_MAKE_EXCEPTION(e);
            e.setTitle(tr("Parsing Error"));
            e.setDetails(tr("Encountered sample mask with invalid length %1. "
                            "Sample size is %2.")
                         .arg(sampleMask.size()).arg(_sampleSize));
            throw e;
         }

         // save previous pair when new pair is read
         Pairwise::Index nextIndex(geneX, geneY);

         if ( index != nextIndex )
         {
            // save pairs
            if ( ccmPair.clusterSize() > 1 )
            {
               ccmPair.write(index);
            }

            if ( cmxPair.clusterSize() > 0 )
            {
               cmxPair.write(index);
            }

            // reset pairs
            ccmPair.clearClusters();
            cmxPair.clearClusters();

            // update index
            index = nextIndex;
         }

         // append cluster and correlation to gene pair
         int cluster = ccmPair.clusterSize();

         ccmPair.addCluster();
         cmxPair.addCluster();

         for ( int i = 0; i < sampleMask.size(); ++i )
         {
            ccmPair.at(cluster, i) = sampleMask[i].digitValue();
         }

         cmxPair.at(cluster, 0) = correlation;
      }
      else if ( words.size() != 1 && words.size() != 0 )
      {
         E_MAKE_EXCEPTION(e);
         e.setTitle(tr("Parsing Error"));
         e.setDetails(tr("Encountered line with incorrect amount of fields. "
                         "Read %1 fields when there should have been %2.")
                      .arg(words.size()).arg(11));
         throw e;
      }

      // figure out percent complete
      int newPercent = 100 * _input->pos() / _input->size();

      // if percent complete has changed update it and emit progressed
      if ( newPercent != lastPercent )
      {
         lastPercent = newPercent;
         emit progressed(lastPercent);
      }
   }

   // make sure reading input file worked
   if ( stream.status() != QTextStream::Ok )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Qt Text Stream encountered an unknown error."));
      throw e;
   }
}
