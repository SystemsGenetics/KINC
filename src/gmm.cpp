#include <ace/core/ace_qmpi.h>
#include <ace/core/metadata.h>

#include "gmm.h"
#include "datafactory.h"



using namespace std;






const char* GMM::BIC {QT_TR_NOOP("BIC")};
const char* GMM::ICL {QT_TR_NOOP("ICL")};






GMM::~GMM()
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

   // delete kernel and program
   delete _kernel;
   delete _program;
}






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
   case MinExpression: return Type::Double;
   case MinClusters: return Type::Integer;
   case MaxClusters: return Type::Integer;
   case CriterionArg: return Type::Combo;
   case RemovePreOutliers: return Type::Bool;
   case RemovePostOutliers: return Type::Bool;
   case BlockSize: return Type::Integer;
   case KernelSize: return Type::Integer;
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
      case MinSamples: return QString("minsamp");
      case MinExpression: return QString("minexpr");
      case MinClusters: return QString("minclus");
      case MaxClusters: return QString("maxclus");
      case CriterionArg: return QString("crit");
      case RemovePreOutliers: return QString("preout");
      case RemovePostOutliers: return QString("postout");
      case BlockSize: return QString("bsize");
      case KernelSize: return QString("ksize");
      default: return QVariant();
      }
   case Role::Title:
      // figure out which argument is being queried and return title
      switch (argument)
      {
      case InputData: return tr("Input:");
      case OutputData: return tr("Output:");
      case MinSamples: return tr("Minimum Sample Size:");
      case MinExpression: return tr("Minimum Expression:");
      case MinClusters: return tr("Minimum Clusters:");
      case MaxClusters: return tr("Maximum Clusters:");
      case CriterionArg: return tr("Criterion:");
      case RemovePreOutliers: return tr("Remove pre-clustering outliers:");
      case RemovePostOutliers: return tr("Remove post-clustering outliers:");
      case BlockSize: return tr("Block Size:");
      case KernelSize: return tr("Kernel Size:");
      default: return QVariant();
      }
   case Role::WhatsThis:
      // figure out which argument is being queried and return "What's This?" text
      switch (argument)
      {
      case InputData: return tr("Input expression matrix that will be used to compute pair-wise clusters.");
      case OutputData: return tr("Output cluster matrix that will store pair-wise clusters.");
      case MinSamples: return tr("Minimum size of samples two genes must share to perform clustering.");
      case MinExpression: return tr("Minimum threshold for a gene expression to be included in clustering.");
      case MinClusters: return tr("Minimum number of clusters to test.");
      case MaxClusters: return tr("Maximum number of clusters to test.");
      case CriterionArg: return tr("Criterion to determine the number of clusters.");
      case RemovePreOutliers: tr("Remove statistical outliers before clustering.");
      case RemovePostOutliers: tr("Remove statistical outliers after clustering.");
      case BlockSize: return tr("This option only applies if OpenCL is used. Total number of blocks"
                                " to run for execution.");
      case KernelSize: return tr("This option only applies if OpenCL is used. Total number of"
                                 " kernels to run per block of execution.");
      default: return QVariant();
      }
   case Role::ComboValues:
      // if this is criterion argument return combo values else return nothing
      switch (argument)
      {
      case CriterionArg: return QStringList() << tr(BIC) << tr(ICL);
      default: return QStringList();
      }
   case Role::DefaultValue:
      // figure out which argument is being queried and if applicable return default value else
      // return nothing
      switch (argument)
      {
      case MinSamples: return 30;
      case MinExpression: return -INFINITY;
      case MinClusters: return 1;
      case MaxClusters: return 5;
      case CriterionArg: return tr(BIC);
      case RemovePreOutliers: return false;
      case RemovePostOutliers: return false;
      case BlockSize: return 4;
      case KernelSize: return 4096;
      default: return QVariant();
      }
   case Role::Minimum:
      // figure out which argument is being queried and if applicable return minimum value else
      // return nothing
      switch (argument)
      {
      case MinSamples: return 1;
      case MinExpression: return -INFINITY;
      case MinClusters: return 1;
      case MaxClusters: return 1;
      case BlockSize: return 1;
      case KernelSize: return 1;
      default: return QVariant();
      }
   case Role::Maximum:
      // figure out which argument is being queried and if applicable return maximum value else
      // return nothing
      switch (argument)
      {
      case MinSamples: return INT_MAX;
      case MinExpression: return +INFINITY;
      case MinClusters: return GenePair::Vector::MAX_CLUSTER_SIZE;
      case MaxClusters: return GenePair::Vector::MAX_CLUSTER_SIZE;
      case BlockSize: return INT_MAX;
      case KernelSize: return INT_MAX;
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
   case MinExpression:
      _minExpression = value.toDouble();
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
   case BlockSize:
      _blockSize = value.toInt();
      break;
   case KernelSize:
      _kernelSize = value.toInt();
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

   // initialize new cc matrix
   _output->initialize(_input->getGeneNames(),_input->getSampleNames());

   // initialize total steps
   _totalSteps = (qint64)_output->geneSize()*(_output->geneSize() - 1)/2;

   // return pre-allocation argument
   return false;
}






void GMM::savePair(GenePair::Vector vector, qint8 K, const qint8 *labels, int N)
{
   // get MPI instance
   Ace::QMPI& mpi {Ace::QMPI::initialize()};

   if ( mpi.isMaster() )
   {
      // save cluster data to output file
      CCMatrix::Pair pair(_output);
      pair.addCluster(K);

      for ( int i = 0; i < N; ++i )
      {
         for ( qint8 k = 0; k < K; ++k )
         {
            pair.at(k, i) = (labels[i] >= 0)
               ? (k == labels[i])
               : -labels[i];
         }
      }

      pair.write(vector);
   }
   else
   {
      // send cluster data to master node
      _mpiOut->writeRawData(reinterpret_cast<const char *>(labels), N);
   }
}






void GMM::runSerial()
{
   // get MPI instance
   Ace::QMPI& mpi {Ace::QMPI::initialize()};

   // initialize clustering model
   GenePair::GMM clusteringModel;

   // iterate through gene pairs
   while ( _stepsComplete < _totalSteps )
   {
      // make sure interruption is not requested
      if ( isInterruptionRequested() )
      {
         return;
      }

      // compute clustering model
      clusteringModel.compute(
         _input,
         _vector,
         _minSamples,
         _minExpression,
         _minClusters,
         _maxClusters,
         _criterion,
         _removePreOutliers,
         _removePostOutliers
      );

      qint8 K {clusteringModel.clusterSize()};
      const QVector<qint8>& labels {clusteringModel.labels()};

      // send cluster size to master node
      if ( !mpi.isMaster() )
      {
         (*_mpiOut) << K;
      }

      // save cluster pair if multiple clusters are found
      if ( K > 1 )
      {
         savePair(_vector, K, labels.data(), labels.size());
      }

      // increment to next pair
      ++_vector;

      // increment steps and calculate percent complete
      ++_stepsComplete;
      qint64 newPercent {100*_stepsComplete/_totalSteps};

      // emit progressed signal when percent changes
      if ( newPercent != _lastPercent )
      {
         _lastPercent = newPercent;
         emit progressed(_lastPercent);
      }
   }
}






QByteArray GMM::buildMPIBlock()
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






bool GMM::readMPIBlock(const QByteArray& block)
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

   qInfo("%d: reading block (%d KB)...", 0, block.size() / 1024);

   // iterate through gene pairs in this block
   QVector<qint8> labels(_input->getSampleSize());

   for ( int i = 0; i < blockSize; ++i )
   {
      // read results
      qint8 K;
      stream >> K;

      // save results if necessary
      if ( K > 1 )
      {
         stream.readRawData(reinterpret_cast<char *>(labels.data()), labels.size());

         savePair(_vector, K, labels.data(), labels.size());
      }

      // increment gene pair index
      ++_vector;
      ++_stepsComplete;
   }

   qInfo("%d: finished reading block", 0);

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






QByteArray GMM::processMPIBlock(const QByteArray& block)
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

   // initialize gene pair vector and total steps
   _vector = GenePair::Vector(blockStart);
   _nextVector = _vector;
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

   qInfo("%d: finished block", mpi.rank());

   // send data to master node
   return data;
}






int GMM::getBlockSize()
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
      initializeKernel();
      initializeBlockExpressions();
      initializeKernelArguments();
   }

   // reset blocks to initial state
   for ( int i = 0; i < _blockSize; ++i )
   {
      _blocks[i]->state = Block::Start;
   }

   // return block size
   return _blockSize;
}






bool GMM::runBlock(int index)
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
   case Block::Execute:
      runExecuteBlock(*_blocks[index]);
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
   int newPercent = _stepsComplete*100/_totalSteps;

   // if percent complete has changed update it and emit progressed
   if ( newPercent != _lastPercent )
   {
      _lastPercent = newPercent;
      emit progressed(_lastPercent);
   }

   // signal block is still running
   return true;
}






void GMM::initializeKernel()
{
   // get opencl device and make program
   EOpenCLDevice& device {EOpenCLDevice::getInstance()};
   _program = device.makeProgram().release();

   // make sure it worked
   if ( !device )
   {
      E_MAKE_EXCEPTION(e);
      device.fillException(e);
      throw e;
   }

   // add opencl c code and compile it making sure it worked
   _program->addFile(":/opencl/gmm.cl");
   if ( !_program->compile() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("OpenCL Compile Error"));
      e.setDetails(tr("OpenCL program failed to compile:\n\n%1")
         .arg(_program->getBuildError())
         .replace('\n',"<br/>"));

      throw e;
   }

   // get kernel from compiled code making sure it worked
   _kernel = _program->makeKernel("computeGMMBlock").release();
   if ( !*_program )
   {
      E_MAKE_EXCEPTION(e);
      _program->fillException(e);
      throw e;
   }
}






void GMM::initializeBlockExpressions()
{
   // get opencl device
   EOpenCLDevice& device {EOpenCLDevice::getInstance()};

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
}






void GMM::initializeKernelArguments()
{
   // get opencl device
   EOpenCLDevice& device {EOpenCLDevice::getInstance()};

   // increase kernel size to first power of 2 reached
   int pow2 {2};
   while ( pow2 < _kernelSize )
   {
      pow2 *= 2;
   }
   _kernelSize = pow2;

   // initialize blocks
   _blocks = new Block*[_blockSize];
   for (int i = 0; i < _blockSize; ++i)
   {
      _blocks[i] = new Block(device, _input->getSampleSize(), _maxClusters, _kernelSize);
   }

   // figure out workgroup size for opencl kernels
   int workgroupSize {_kernelSize};
   while ( workgroupSize > (int)_kernel->getMaxWorkgroupSize() )
   {
      workgroupSize /= 2;
   }

   // set all static kernel arguments
   _kernel->setBuffer(0, _expressions);
   _kernel->setArgument(1, (cl_int)_input->getSampleSize());
   _kernel->setArgument(3, (cl_int)_minSamples);
   _kernel->setArgument(4, (cl_int)_minExpression);
   _kernel->setArgument(5, (cl_char)_minClusters);
   _kernel->setArgument(6, (cl_char)_maxClusters);
   _kernel->setArgument(7, (cl_int)_criterion);
   _kernel->setArgument(8, (cl_int)_removePreOutliers);
   _kernel->setArgument(9, (cl_int)_removePostOutliers);
   _kernel->setDimensionCount(1);
   _kernel->setGlobalSize(0, _kernelSize);
   _kernel->setWorkgroupSize(0, workgroupSize);

   // make sure everything with kernel worked
   if ( !*_kernel )
   {
      E_MAKE_EXCEPTION(e);
      _kernel->fillException(e);
      throw e;
   }
}






void GMM::runStartBlock(Block& block)
{
   Ace::QMPI& mpi {Ace::QMPI::initialize()};

   // check if there are more gene pairs to compute
   if ( _stepsStarted < _totalSteps )
   {
      qInfo("%d: loading %12lld / %lld...", mpi.rank(), _stepsStarted, _totalSteps);

      // set block xy to beginning of read comparisons
      block.vector = _vector;

      // copy list of pairs to be done on this run
      int index {0};
      while ( _stepsStarted + index < _totalSteps && index < _kernelSize )
      {
         (*block.pairs)[index] = { _vector.geneX(), _vector.geneY() };
         ++_vector;
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






void GMM::runLoadBlock(Block& block)
{
   Ace::QMPI& mpi {Ace::QMPI::initialize()};

   // check to see if reference loading is complete
   if ( !block.isWaiting() )
   {
      qInfo("%d: executing kernel...", mpi.rank());

      // make sure opencl events worked
      block.checkAllEvents();

      // set kernel arguments and execute it
      _kernel->setBuffer(2, block.pairs);
      _kernel->setBuffer(10, block.work_X);
      _kernel->setBuffer(11, block.work_labels);
      _kernel->setBuffer(12, block.work_components);
      _kernel->setBuffer(13, block.work_MP);
      _kernel->setBuffer(14, block.work_counts);
      _kernel->setBuffer(15, block.work_logpi);
      _kernel->setBuffer(16, block.work_loggamma);
      _kernel->setBuffer(17, block.work_logGamma);
      _kernel->setBuffer(18, block.result_K);
      _kernel->setBuffer(19, block.result_labels);
      block.events.append(_kernel->execute());

      // make sure kernel worked
      if ( !*_kernel )
      {
         E_MAKE_EXCEPTION(e);
         _kernel->fillException(e);
         throw e;
      }

      // change block state to execute
      block.state = Block::Execute;
   }
}






void GMM::runExecuteBlock(Block& block)
{
   Ace::QMPI& mpi {Ace::QMPI::initialize()};

   // check to see if kernel execution is complete
   if ( !block.isWaiting() )
   {
      qInfo("%d: reading data...", mpi.rank());

      // make sure opencl events worked
      block.checkAllEvents();

      // read clustering results from device and make sure it worked
      block.events.append(block.result_K->read());
      if ( !*block.result_K )
      {
         E_MAKE_EXCEPTION(e);
         block.result_K->fillException(e);
         throw e;
      }

      block.events.append(block.result_labels->read());
      if ( !*block.result_labels )
      {
         E_MAKE_EXCEPTION(e);
         block.result_labels->fillException(e);
         throw e;
      }

      // change block state to read
      block.state = Block::Read;
   }
}






void GMM::runReadBlock(Block& block)
{
   // get MPI instance
   Ace::QMPI& mpi {Ace::QMPI::initialize()};

   // check if read is complete and we are next in line
   if ( !block.isWaiting() && _nextVector == block.vector )
   {
      qInfo("%d: reading %12lld / %lld...", mpi.rank(), _stepsComplete, _totalSteps);

      // make sure opencl events worked
      block.checkAllEvents();

      // save each valid clustering result to cluster matrix
      int index {0};
      while ( _stepsComplete < _totalSteps && index < _kernelSize )
      {
         // read results
         int N = _input->getSampleSize();
         qint8 K = (*block.result_K)[index];
         qint8 *labels = &(*block.result_labels)[index * N];

         if ( !mpi.isMaster() )
         {
            (*_mpiOut) << K;
         }

         // save cluster pair if multiple clusters are found
         if ( K > 1 )
         {
            savePair(block.vector, K, labels, N);
         }

         // increment indices
         ++block.vector;
         ++index;
         ++_stepsComplete;
      }

      // update next vector and change block state to start
      _nextVector = block.vector;
      block.state = Block::Start;

      // clear all opencl events
      block.events.clear();
   }
}
