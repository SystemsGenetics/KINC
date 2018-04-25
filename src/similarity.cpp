#include <ace/core/ace_qmpi.h>

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






Similarity::~Similarity()
{
   // check if blocks are allocated
   if ( _blocks )
   {
      // iterate through all blocks and delete them
      for (int i = 0; i < _blockSize; ++i)
      {
         delete _blocks[i];
      }

      // delete pointer list
      delete[] _blocks;
   }

   // delete program, kernels, and buffers
   delete _program;
   delete _kernel1;
   delete _kernel2;
   delete _expressions;
}






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
   case ClusteringArg: return Type::Combo;
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
   case BlockSize: return Type::Integer;
   case KernelSize: return Type::Integer;
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
      case ClusteringArg: return QString("clusmethod");
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
      case BlockSize: return QString("bsize");
      case KernelSize: return QString("ksize");
      default: return QVariant();
      }
   case Role::Title:
      // figure out which argument is being queried and return title
      switch (argument)
      {
      case InputData: return tr("Expressiom Matrix:");
      case ClusterData: return tr("Cluster Matrix:");
      case CorrelationData: return tr("Correlation Matrix:");
      case ClusteringArg: return tr("Clustering Method:");
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
      case BlockSize: return tr("Block Size:");
      case KernelSize: return tr("Kernel Size:");
      default: return QVariant();
      }
   case Role::WhatsThis:
      // figure out which argument is being queried and return "What's This?" text
      switch (argument)
      {
      case InputData: return tr("Input expression matrix.");
      case ClusterData: return tr("Output matrix that will contain gene pair clusters.");
      case CorrelationData: return tr("Output matrix that will contain gene pair correlations.");
      case ClusteringArg: return tr("Clustering method to use for gene pairs.");
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
      case BlockSize: return tr("(OpenCL) Total number of blocks to run.");
      case KernelSize: return tr("(OpenCL) Total number of kernels to run per block.");
      default: return QVariant();
      }
   case Role::ComboValues:
      // if this is criterion argument return combo values else return nothing
      switch (argument)
      {
      case ClusteringArg: return QStringList({ tr(GMM), tr(KMeans) });
      case CorrelationArg: return QStringList({ tr(Pearson), tr(Spearman) });
      case CriterionArg: return QStringList({ tr(BIC), tr(ICL) });
      default: return QStringList();
      }
   case Role::DefaultValue:
      // figure out which argument is being queried and if applicable return default value else
      // return nothing
      switch (argument)
      {
      case ClusteringArg: return tr(GMM);
      case CorrelationArg: return tr(Pearson);
      case MinExpression: return -INFINITY;
      case MinSamples: return 30;
      case MinClusters: return 1;
      case MaxClusters: return 5;
      case CriterionArg: return tr(ICL);
      case RemovePreOutliers: return false;
      case RemovePostOutliers: return false;
      case MinCorrelation: return 0.5;
      case MaxCorrelation: return 1.0;
      case BlockSize: return 4;
      case KernelSize: return 4096;
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
      case BlockSize: return 1;
      case KernelSize: return 1;
      default: return QVariant();
      }
   case Role::Maximum:
      // figure out which argument is being queried and if applicable return maximum value else
      // return nothing
      switch (argument)
      {
      case MinExpression: return +INFINITY;
      case MinSamples: return INT_MAX;
      case MinClusters: return GenePair::Index::MAX_CLUSTER_SIZE;
      case MaxClusters: return GenePair::Index::MAX_CLUSTER_SIZE;
      case MinCorrelation: return 1.0;
      case MaxCorrelation: return 1.0;
      case BlockSize: return INT_MAX;
      case KernelSize: return INT_MAX;
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
   case ClusteringArg:
      {
         const QString option = value.toString();
         if ( option == tr(GMM) )
         {
            _clusMethod = ClusteringMethod::GMM;
            _clusModel = new GenePair::GMM();
         }
         else if ( option == tr(KMeans) )
         {
            _clusMethod = ClusteringMethod::KMeans;
            _clusModel = new GenePair::KMeans();
         }
      }
      break;
   case CorrelationArg:
      {
         const QString option = value.toString();
         if ( option == tr(Pearson) )
         {
            _corrMethod = CorrelationMethod::Pearson;
            _corrModel = new GenePair::Pearson();
         }
         else if ( option == tr(Spearman) )
         {
            _corrMethod = CorrelationMethod::Spearman;
            _corrModel = new GenePair::Spearman();
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
   case BlockSize:
      _blockSize = value.toInt();
      break;
   case KernelSize:
      _kernelSize = value.toInt();
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
   _clusMatrix->initialize(_input->getGeneNames(), _input->getSampleNames());

   // initialize correlation matrix
   EMetadata correlations(EMetadata::Array);
   EMetadata* name {new EMetadata(EMetadata::String)};
   *(name->toString()) = _corrModel->getName();
   correlations.toArray()->append(name);

   _corrMatrix->initialize(_input->getGeneNames(), correlations);

   // initialize total steps
   _totalSteps = (qint64) _input->getGeneSize() * (_input->getGeneSize() - 1) / 2;

   // return pre-allocation argument
   return false;
}






int Similarity::fetchPair(GenePair::Index index, QVector<GenePair::Vector2>& X, QVector<qint8>& labels)
{
   // read in gene expressions
   ExpressionMatrix::Gene gene1(_input);
   ExpressionMatrix::Gene gene2(_input);

   gene1.read(index.getX());
   gene2.read(index.getY());

   // populate X with shared expressions of gene pair
   int numSamples = 0;

   for ( int i = 0; i < _input->getSampleSize(); ++i )
   {
      if ( std::isnan(gene1.at(i)) || std::isnan(gene2.at(i)) )
      {
         labels[i] = -9;
      }
      else if ( gene1.at(i) < _minExpression || gene2.at(i) < _minExpression )
      {
         labels[i] = -6;
      }
      else
      {
         X[numSamples] = { gene1.at(i), gene2.at(i) };
         numSamples++;

         labels[i] = 0;
      }
   }

   // return size of X
   return numSamples;
}






void Similarity::savePair(GenePair::Index index, qint8 K, const qint8 *labels, int N, const float *correlations)
{
   // get MPI instance
   Ace::QMPI& mpi {Ace::QMPI::initialize()};

   if ( mpi.isMaster() )
   {
      // save clusters whose correlations are within thresholds
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
            clusPair.write(index);
         }
      }

      // save correlations that are within thresholds
      if ( K > 0 )
      {
         // save correlation data to output file
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
            corrPair.write(index);
         }
      }
   }
   else
   {
      // send cluster size to master node
      (*_mpiOut) << K;

      // send cluster data to master node
      if ( K > 1 )
      {
         _mpiOut->writeRawData(reinterpret_cast<const char *>(labels), N);
      }

      // send correlation data to master node
      if ( K > 0 )
      {
         _mpiOut->writeRawData(reinterpret_cast<const char *>(correlations), K * sizeof(float));
      }
   }

}






void Similarity::runSerial()
{
   // initialize clustering model
   _clusModel->initialize(_input);

   // initialize correlation model
   _corrModel->initialize(_input);

   // initialize workspace
   QVector<GenePair::Vector2> X(_input->getSampleSize());
   QVector<qint8> labels(_input->getSampleSize());

   // iterate through all gene pairs
   while ( _stepsComplete < _totalSteps )
   {
      // make sure interruption is not requested
      if ( isInterruptionRequested() )
      {
         return;
      }

      // fetch pairwise data
      int numSamples = fetchPair(_index, X, labels);

      // compute clusters
      qint8 K = _clusModel->compute(
         X,
         numSamples,
         labels,
         _minSamples,
         _minClusters,
         _maxClusters,
         _criterion,
         _removePreOutliers,
         _removePostOutliers
      );

      // compute correlation
      QVector<float> correlations = _corrModel->compute(
         X,
         K,
         labels,
         _minSamples
      );

      // save gene pair data
      savePair(_index, K, labels.data(), labels.size(), correlations.data());

      // increment to next pair
      ++_index;

      // increment steps and calculate percent complete
      ++_stepsComplete;
      qint64 newPercent {100 * _stepsComplete / _totalSteps};

      // emit progressed signal when percent changes
      if ( newPercent != _lastPercent )
      {
         _lastPercent = newPercent;
         emit progressed(_lastPercent);
      }
   }
}






QByteArray Similarity::buildMPIBlock()
{
   const qint64 MPI_BLOCK_SIZE { 32 * 1024 };

   // initialize output data
   QByteArray data;
   QDataStream stream(&data, QIODevice::WriteOnly);

   // check if there are more gene pairs to compute
   if ( _stepsStarted < _totalSteps )
   {
      // determine block size
      qint64 steps = min(_totalSteps - _stepsStarted, MPI_BLOCK_SIZE);

      // write block start and block size
      stream << _stepsStarted << steps;

      // update steps started
      _stepsStarted += steps;
   }

   // send data to worker
   return data;
}






bool Similarity::readMPIBlock(const QByteArray& block)
{
   // read block start and block size from worker
   QDataStream stream(block);

   qint64 blockStart;
   qint64 blockSize;
   stream >> blockStart >> blockSize;

   // defer this block if it is not next in line
   if ( blockStart != _stepsComplete )
   {
      return false;
   }

   // iterate through gene pairs in this block
   QVector<qint8> labels(_input->getSampleSize());
   QVector<float> correlations(_maxClusters);

   for ( int i = 0; i < blockSize; ++i )
   {
      // read cluster size
      qint8 K;
      stream >> K;

      // read cluster data
      if ( K > 1 )
      {
         stream.readRawData(reinterpret_cast<char *>(labels.data()), labels.size());
      }

      // read correlation data
      if ( K > 0 )
      {
         stream.readRawData(reinterpret_cast<char *>(correlations.data()), K * sizeof(float));
      }

      // save results
      savePair(_index, K, labels.data(), labels.size(), correlations.data());

      // increment gene pair index
      ++_index;
      ++_stepsComplete;
   }

   // calculate percent complete
   qint64 newPercent {100*_stepsComplete/_totalSteps};

   // emit progressed signal when percent changes
   if ( newPercent != _lastPercent )
   {
      _lastPercent = newPercent;
      emit progressed(_lastPercent);
   }

   return true;
}






QByteArray Similarity::processMPIBlock(const QByteArray& block)
{
   // get MPI instance
   Ace::QMPI& mpi {Ace::QMPI::initialize()};

   // read block start and block size from master
   QDataStream stream(block);

   qint64 blockStart;
   qint64 blockSize;
   stream >> blockStart >> blockSize;

   qInfo("%d: starting block %12lld", mpi.rank(), blockStart);

   // initialize output data
   QByteArray data;

   _mpiOut = new QDataStream(&data, QIODevice::WriteOnly);

   // write block start and block size
   (*_mpiOut) << blockStart << blockSize;

   // initialize pairwise index and total steps
   _index = GenePair::Index(blockStart);
   _nextIndex = _index;
   _totalSteps = blockSize;
   _stepsStarted = 0;
   _stepsComplete = 0;

   // check to see if analytic can run OpenCL and there is a device to use
   if ( getCapabilities()&Capabilities::OpenCL
        && EOpenCLDevice::getInstance().getStatus() == EOpenCLDevice::Ok )
   {
      // initialize block info
      int blockSize {getBlockSize()};
      int done {0};
      bool blocks[blockSize];
      for (int i = 0; i < blockSize ;++i)
      {
         blocks[i] = true;
      }

      // begin block while loop
      while ( done < blockSize )
      {
         // if interruption is requested exit
         if ( isInterruptionRequested() )
         {
            break;
         }

         for (int i = 0; i < blockSize ;++i)
         {
            // make sure block is still alive
            if ( blocks[i] )
            {
               // if block is still alive run it
               if ( !runBlock(i) )
               {
                  // block is done, remove it from active list
                  blocks[i] = false;
                  ++done;
               }
            }
         }
      }
   }

   // else just run serial if possible
   else if ( getCapabilities()&Capabilities::Serial )
   {
      runSerial();
   }

   // else analytic cannot run and throw failure
   else
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Failed Analytic Execution."));
      e.setDetails(tr("Could not execute analytic because it lacks any applicable capability."));
      throw e;
   }

   // cleanup
   delete _mpiOut;

   // send data to master node
   return data;
}






int Similarity::getBlockSize()
{
   // make sure block and kernel size are greater than zero
   if ( _blockSize < 1 || _kernelSize < 1 )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails((tr("Block size and/or kernel size are set to values less than 1.")));
      throw e;
   }

   // initialize all opencl components
   if ( !_blocks )
   {
      initializeOpenCL();
   }

   // reset blocks to initial state
   for ( int i = 0; i < _blockSize; ++i )
   {
      _blocks[i]->state = Block::Start;
   }

   // return block size
   return _blockSize;
}






bool Similarity::runBlock(int index)
{
   // figure out what state this block is in and execute it
   switch (_blocks[index]->state)
   {
   case Block::Start:
      runStartBlock(*_blocks[index]);
      break;
   case Block::Load:
      runLoadBlock(*_blocks[index]);
      break;
   case Block::Execute1:
      runExecute1Block(*_blocks[index]);
      break;
   case Block::Execute2:
      runExecute2Block(*_blocks[index]);
      break;
   case Block::Read:
      runReadBlock(*_blocks[index]);
      break;
   case Block::Done:
   default:
      // if state is done or unknown signal this block is done
      return false;
   }

   // figure out new percent complete
   qint64 newPercent {100 * _stepsComplete / _totalSteps};

   // emit progressed signal when percent changes
   if ( newPercent != _lastPercent )
   {
      _lastPercent = newPercent;
      emit progressed(_lastPercent);
   }

   // signal block is still running
   return true;
}






void Similarity::initializeOpenCL()
{
   // get opencl device
   EOpenCLDevice& device {EOpenCLDevice::getInstance()};

   // make program making sure it worked
   _program = device.makeProgram().release();
   if ( !device )
   {
      E_MAKE_EXCEPTION(e);
      device.fillException(e);
      throw e;
   }

   // add opencl c code
   _program->addFile(":/opencl/linalg.cl");
   _program->addFile(":/opencl/fetchpair.cl");
   _program->addFile(":/opencl/sort.cl");
   _program->addFile(":/opencl/outlier.cl");

   if ( _clusMethod == ClusteringMethod::GMM )
   {
      _program->addFile(":/opencl/gmm.cl");
   }
   else if ( _clusMethod == ClusteringMethod::KMeans )
   {
      _program->addFile(":/opencl/kmeans.cl");
   }

   if ( _corrMethod == CorrelationMethod::Pearson )
   {
      _program->addFile(":/opencl/pearson.cl");
   }
   else if ( _corrMethod == CorrelationMethod::Spearman )
   {
      _program->addFile(":/opencl/spearman.cl");
   }

   // compile program making sure it worked
   if ( !_program->compile() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("OpenCL Compile Error"));
      e.setDetails(tr("OpenCL program failed to compile:\n\n%1")
         .arg(_program->getBuildError())
         .replace('\n',"<br/>"));

      throw e;
   }

   // get clustering kernel from compiled code making sure it worked
   QString clusKernel;

   if ( _clusMethod == ClusteringMethod::GMM )
   {
      clusKernel = "computeGMMBlock";
   }
   else if ( _clusMethod == ClusteringMethod::KMeans )
   {
      clusKernel = "computeKmeansBlock";
   }

   _kernel1 = _program->makeKernel(clusKernel).release();
   if ( !*_program )
   {
      E_MAKE_EXCEPTION(e);
      _program->fillException(e);
      throw e;
   }

   // get correlation kernel from compiled code making sure it worked
   QString corrKernel;

   if ( _corrMethod == CorrelationMethod::Pearson )
   {
      corrKernel = "computePearsonBlock";
   }
   else if ( _corrMethod == CorrelationMethod::Spearman )
   {
      corrKernel = "computeSpearmanBlock";
   }

   _kernel2 = _program->makeKernel(corrKernel).release();
   if ( !*_program )
   {
      E_MAKE_EXCEPTION(e);
      _program->fillException(e);
      throw e;
   }

   // make new opencl buffer for expressions making sure it worked
   _expressions = device.makeBuffer<cl_float>(_input->getRawSize()).release();
   if ( !device )
   {
      E_MAKE_EXCEPTION(e);
      device.fillException(e);
      throw e;
   }

   // get raw expression data from input
   unique_ptr<ExpressionMatrix::Expression> rawData(_input->dumpRawData());
   ExpressionMatrix::Expression* rawDataRef {rawData.get()};

   // copy expression data to opencl buffer
   for (int i = 0; i < _input->getRawSize(); ++i)
   {
      (*_expressions)[i] = rawDataRef[i];
   }

   // write opencl expression buffer to device making sure it worked
   EOpenCLEvent event = _expressions->write();
   if ( !*_expressions )
   {
      E_MAKE_EXCEPTION(e);
      _expressions->fillException(e);
      throw e;
   }

   // wait for write to finish making sure it worked
   event.wait();
   if ( !event )
   {
      E_MAKE_EXCEPTION(e);
      event.fillException(e);
      throw e;
   }

   // get work size (sample size rounded up to next power of 2)
   int N_pow2 {2};
   while ( N_pow2 < _input->getSampleSize() )
   {
      N_pow2 *= 2;
   }

   // increase kernel size to next power of 2
   int pow2 {2};
   while ( pow2 < _kernelSize )
   {
      pow2 *= 2;
   }
   _kernelSize = pow2;

   // figure out workgroup size for clustering kernel
   int workgroupSize1 {_kernelSize};
   while ( workgroupSize1 > (int)_kernel1->getMaxWorkgroupSize() )
   {
      workgroupSize1 /= 2;
   }

   // set all static kernel arguments
   _kernel1->setDimensionCount(1);
   _kernel1->setGlobalSize(0, _kernelSize);
   _kernel1->setWorkgroupSize(0, workgroupSize1);

   if ( _clusMethod == ClusteringMethod::GMM )
   {
      _kernel1->setBuffer(0, _expressions);
      _kernel1->setArgument(1, (cl_int)_input->getSampleSize());
      _kernel1->setArgument(3, (cl_int)_minSamples);
      _kernel1->setArgument(4, (cl_int)_minExpression);
      _kernel1->setArgument(5, (cl_char)_minClusters);
      _kernel1->setArgument(6, (cl_char)_maxClusters);
      _kernel1->setArgument(7, (cl_int)_criterion);
      _kernel1->setArgument(8, (cl_int)_removePreOutliers);
      _kernel1->setArgument(9, (cl_int)_removePostOutliers);
   }
   else if ( _clusMethod == ClusteringMethod::KMeans )
   {
      _kernel1->setBuffer(0, _expressions);
      _kernel1->setArgument(1, (cl_int)_input->getSampleSize());
      _kernel1->setArgument(3, (cl_int)_minSamples);
      _kernel1->setArgument(4, (cl_int)_minExpression);
      _kernel1->setArgument(5, (cl_char)_minClusters);
      _kernel1->setArgument(6, (cl_char)_maxClusters);
      _kernel1->setArgument(7, (cl_int)_removePreOutliers);
      _kernel1->setArgument(8, (cl_int)_removePostOutliers);
   }

   // make sure everything with kernel worked
   if ( !*_kernel1 )
   {
      E_MAKE_EXCEPTION(e);
      _kernel1->fillException(e);
      throw e;
   }

   // figure out workgroup size for correlation kernel
   int workgroupSize2 {_kernelSize};
   while ( workgroupSize2 > (int)_kernel2->getMaxWorkgroupSize() )
   {
      workgroupSize2 /= 2;
   }

   // set all static kernel arguments
   _kernel2->setDimensionCount(1);
   _kernel2->setGlobalSize(0, _kernelSize);
   _kernel2->setWorkgroupSize(0, workgroupSize2);

   if ( _corrMethod == CorrelationMethod::Pearson )
   {
      _kernel2->setArgument(1, (cl_char)_maxClusters);
      _kernel2->setArgument(3, (cl_int)_input->getSampleSize());
      _kernel2->setArgument(4, (cl_int)_minSamples);
   }
   else if ( _corrMethod == CorrelationMethod::Spearman )
   {
      _kernel2->setArgument(1, (cl_char)_maxClusters);
      _kernel2->setArgument(3, (cl_int)_input->getSampleSize());
      _kernel2->setArgument(4, (cl_int)_minSamples);
   }

   // make sure everything with kernel worked
   if ( !*_kernel2 )
   {
      E_MAKE_EXCEPTION(e);
      _kernel2->fillException(e);
      throw e;
   }

   // initialize blocks
   _blocks = new Block*[_blockSize];
   for (int i = 0; i < _blockSize; ++i)
   {
      _blocks[i] = new Block(device, _input->getSampleSize(), N_pow2, _maxClusters, _kernelSize);
   }
}






void Similarity::runStartBlock(Block& block)
{
   // check if there are more gene pairs to compute
   if ( _stepsStarted < _totalSteps )
   {
      // set block xy to beginning of read comparisons
      block.index = _index;

      // copy list of pairs to be done on this run
      int index {0};
      while ( _stepsStarted + index < _totalSteps && index < _kernelSize )
      {
         (*block.pairs)[index] = { _index.getX(), _index.getY() };
         ++_index;
         ++index;
      }

      // update steps started
      _stepsStarted += index;

      // copy any remaining and unused pairs to zero
      while ( index < _kernelSize )
      {
         (*block.pairs)[index] = { 0, 0 };
         ++index;
      }

      // write pairs to device making sure it worked
      block.events.append(block.pairs->write());
      if ( !*block.pairs )
      {
         E_MAKE_EXCEPTION(e);
         block.pairs->fillException(e);
         throw e;
      }

      // change block state to load
      block.state = Block::Load;
   }

   // else all pairs are complete and this block is done
   else
   {
      block.state = Block::Done;
   }
}






void Similarity::runLoadBlock(Block& block)
{
   // check to see if reference loading is complete
   if ( !block.isWaiting() )
   {
      // make sure opencl events worked
      block.checkAllEvents();

      // set kernel arguments and execute it
      if ( _clusMethod == ClusteringMethod::GMM )
      {
         _kernel1->setBuffer(2, block.pairs);
         _kernel1->setBuffer(10, block.work_X);
         _kernel1->setBuffer(11, block.work_labels);
         _kernel1->setBuffer(12, block.work_components);
         _kernel1->setBuffer(13, block.work_MP);
         _kernel1->setBuffer(14, block.work_counts);
         _kernel1->setBuffer(15, block.work_logpi);
         _kernel1->setBuffer(16, block.work_loggamma);
         _kernel1->setBuffer(17, block.work_logGamma);
         _kernel1->setBuffer(18, block.out_K);
         _kernel1->setBuffer(19, block.out_labels);
      }
      else if ( _clusMethod == ClusteringMethod::KMeans )
      {
         _kernel1->setBuffer(2, block.pairs);
         _kernel1->setBuffer(9, block.work_X);
         _kernel1->setBuffer(10, block.work_loggamma);
         _kernel1->setBuffer(11, block.work_labels);
         _kernel1->setBuffer(12, block.work_MP);
         _kernel1->setBuffer(13, block.out_K);
         _kernel1->setBuffer(14, block.out_labels);
      }

      block.events.append(_kernel1->execute());

      // make sure kernel worked
      if ( !*_kernel1 )
      {
         E_MAKE_EXCEPTION(e);
         _kernel1->fillException(e);
         throw e;
      }

      // change block state to execute
      block.state = Block::Execute1;
   }
}






void Similarity::runExecute1Block(Block& block)
{
   // check to see if reference loading is complete
   if ( !block.isWaiting() )
   {
      // make sure opencl events worked
      block.checkAllEvents();

      // set kernel arguments and execute it
      if ( _corrMethod == CorrelationMethod::Pearson )
      {
         _kernel2->setBuffer(0, block.work_X);
         _kernel2->setBuffer(2, block.out_labels);
         _kernel2->setBuffer(5, block.out_correlations);
      }
      else if ( _corrMethod == CorrelationMethod::Spearman )
      {
         _kernel2->setBuffer(0, block.work_X);
         _kernel2->setBuffer(2, block.out_labels);
         _kernel2->setBuffer(5, block.work_x);
         _kernel2->setBuffer(6, block.work_y);
         _kernel2->setBuffer(7, block.work_rank);
         _kernel2->setBuffer(8, block.out_correlations);
      }

      block.events.append(_kernel2->execute());

      // make sure kernel worked
      if ( !*_kernel2 )
      {
         E_MAKE_EXCEPTION(e);
         _kernel2->fillException(e);
         throw e;
      }

      // change block state to execute
      block.state = Block::Execute2;
   }
}






void Similarity::runExecute2Block(Block& block)
{
   // check to see if kernel execution is complete
   if ( !block.isWaiting() )
   {
      // make sure opencl events worked
      block.checkAllEvents();

      // read results from device and make sure it worked
      block.events.append(block.out_K->read());
      if ( !*block.out_K )
      {
         E_MAKE_EXCEPTION(e);
         block.out_K->fillException(e);
         throw e;
      }

      block.events.append(block.out_labels->read());
      if ( !*block.out_labels )
      {
         E_MAKE_EXCEPTION(e);
         block.out_labels->fillException(e);
         throw e;
      }

      block.events.append(block.out_correlations->read());
      if ( !*block.out_correlations )
      {
         E_MAKE_EXCEPTION(e);
         block.out_correlations->fillException(e);
         throw e;
      }

      // change block state to read
      block.state = Block::Read;
   }
}






void Similarity::runReadBlock(Block& block)
{
   // check if read is complete and we are next in line
   if ( !block.isWaiting() && _nextIndex == block.index )
   {
      // make sure opencl events worked
      block.checkAllEvents();

      // save each valid clustering result to cluster matrix
      int index {0};
      while ( _stepsComplete < _totalSteps && index < _kernelSize )
      {
         // read results
         int N = _input->getSampleSize();
         qint8 K = (*block.out_K)[index];
         qint8 *labels = &(*block.out_labels)[index * N];
         float *correlations = &(*block.out_correlations)[index * _maxClusters];

         // save gene pair
         savePair(block.index, K, labels, N, correlations);

         // increment indices
         ++block.index;
         ++index;
         ++_stepsComplete;
      }

      // clear all opencl events
      block.events.clear();

      // update next index and change block state to start
      _nextIndex = block.index;
      block.state = Block::Start;
   }
}
