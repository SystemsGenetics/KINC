#include <ace/core/metadata.h>

#include "kmeans.h"
#include "datafactory.h"



using namespace std;






KMeans::~KMeans()
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






EAbstractAnalytic::ArgumentType KMeans::getArgumentData(int argument)
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
   case BlockSize: return Type::Integer;
   case KernelSize: return Type::Integer;
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
      case MinSamples: return QString("min");
      case MinClusters: return QString("minclus");
      case MaxClusters: return QString("maxclus");
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
      case MinClusters: return tr("Minimum Clusters:");
      case MaxClusters: return tr("Maximum Clusters:");
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
      case MinClusters: return tr("Minimum number of clusters to test.");
      case MaxClusters: return tr("Maximum number of clusters to test.");
      case BlockSize: return tr("This option only applies if OpenCL is used. Total number of blocks"
                                " to run for execution.");
      case KernelSize: return tr("This option only applies if OpenCL is used. Total number of"
                                 " kernels to run per block of execution.");
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
      case MinClusters: return INT_MAX;
      case MaxClusters: return INT_MAX;
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






void KMeans::setArgument(int argument, QVariant value)
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
   case BlockSize:
      _blockSize = value.toInt();
      break;
   case KernelSize:
      _kernelSize = value.toInt();
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






float KMeans::computeBIC(int K, float logL, int N, int D)
{
   int p = K * D;

   return log(N) * p - 2 * logL;
}






void KMeans::computeModel(const QVector<GenePair::Vector2>& X, int& bestK, QVector<int>& bestLabels)
{
   float bestValue = INFINITY;

   for ( int K = _minClusters; K <= _maxClusters; ++K )
   {
      // run each clustering model
      GenePair::KMeans model;

      model.fit(X, K);

      // evaluate model
      float value = computeBIC(K, model.logLikelihood(), X.size(), 2);

      // save the best model
      if ( value < bestValue )
      {
         bestK = K;
         bestValue = value;

         for ( int i = 0, j = 0; i < X.size(); ++i )
         {
            if ( bestLabels[i] != -1 )
            {
               bestLabels[i] = model.labels()[j];
               ++j;
            }
         }
      }
   }
}






void KMeans::savePair(const GenePair::Vector& vector, int K, const QVector<int>& labels)
{
   CCMatrix::Pair pair(_output);
   pair.addCluster(K);

   for ( int i = 0; i < labels.size(); ++i )
   {
      for ( int k = 0; k < K; ++k )
      {
         pair.at(k, i) = (k == labels[i]);
      }
   }

   pair.write(vector);
}






void KMeans::runSerial()
{
   // initialize percent complete and steps
   int lastPercent {0};
   qint64 steps {0};
   qint64 totalSteps {_output->geneSize()*(_output->geneSize() - 1)/2};

   // initialize arrays used for k-means clustering
   QVector<int> labels(_input->getSampleSize());
   QVector<GenePair::Vector2> X;

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
      gene1.read(vector.geneX());
      gene2.read(vector.geneY());

      // populate X with shared expressions of gene x and y
      X.clear();
      X.reserve(_input->getSampleSize());

      for ( auto i = 0; i < _input->getSampleSize(); ++i )
      {
         if ( !std::isnan(gene1.at(i)) && !std::isnan(gene2.at(i)) )
         {
            X.append({ gene1.at(i), gene2.at(i) });

            labels[i] = 0;
         }
         else
         {
            labels[i] = -1;
         }
      }

      // perform clustering only if there are enough samples
      int bestK = 0;

      if ( X.size() >= _minSamples )
      {
         computeModel(X, bestK, labels);
      }

      // save cluster pair if multiple clusters are found
      if ( bestK > 1 )
      {
         savePair(vector, bestK, labels);
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






int KMeans::getBlockSize()
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
   initializeKernel();
   initializeBlockExpressions();
   initializeKernelArguments();


   // calculate total number of calculations that will be done and return block size
   qint64 geneSize {_output->geneSize()};
   _totalPairs = geneSize*(geneSize-1)/2;
   return _blockSize;
}






bool KMeans::runBlock(int index)
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
   int newPercent = _pairsComplete*100/_totalPairs;

   // if percent complete has changed update it and emit progressed
   if ( newPercent != _lastPercent )
   {
      _lastPercent = newPercent;
      emit progressed(_lastPercent);
   }

   // signal block is still running
   return true;
}






void KMeans::initializeKernel()
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
   _program->addFile(":/opencl/kmeans.cl");
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
   _kernel = _program->makeKernel("computeKmeansBlock").release();
   if ( !*_program )
   {
      E_MAKE_EXCEPTION(e);
      _program->fillException(e);
      throw e;
   }
}






void KMeans::initializeBlockExpressions()
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






void KMeans::initializeKernelArguments()
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
   _kernel->setArgument(4, (cl_int)_minClusters);
   _kernel->setArgument(5, (cl_int)_maxClusters);
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






void KMeans::runStartBlock(Block& block)
{
   // check if there are more gene pairs to compute
   if ( _vector.geneX() < _output->geneSize() )
   {
      // set block xy to beginning of read comparisons
      block.vector = _vector;

      // copy list of pairs to be done on this run
      int index {0};
      while ( _vector.geneX() < _output->geneSize() && index < _kernelSize )
      {
         (*block.pairs)[index] = { _vector.geneX(), _vector.geneY() };
         ++_vector;
         ++index;
      }

      // copy any remaining and unused pairs to zero
      while ( index < _kernelSize )
      {
         (*block.pairs)[index] = { 0, 0 };
         ++index;
      }

      // write pairs to device making sure it worked
      block.event = block.pairs->write();
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






void KMeans::runLoadBlock(Block& block)
{
   // check to see if reference loading is complete
   if ( block.event.isDone() )
   {
      // make sure opencl event worked
      if ( !block.event )
      {
         E_MAKE_EXCEPTION(e);
         block.event.fillException(e);
         throw e;
      }

      // set kernel arguments and execute it
      _kernel->setBuffer(2, block.pairs);
      _kernel->setBuffer(6, block.work_X);
      _kernel->setBuffer(7, block.work_y);
      _kernel->setBuffer(8, block.work_Mu);
      _kernel->setBuffer(9, block.work_ynext);
      _kernel->setBuffer(10, block.result_K);
      _kernel->setBuffer(11, block.result_labels);
      block.event = _kernel->execute();

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






void KMeans::runExecuteBlock(Block& block)
{
   // check to see if kernel execution is complete
   if ( block.event.isDone() )
   {
      // make sure opencl event worked
      if ( !block.event )
      {
         E_MAKE_EXCEPTION(e);
         block.event.fillException(e);
         throw e;
      }

      // read clustering results from device and make sure it worked
      block.event = block.result_K->read();
      if ( !*block.result_K )
      {
         E_MAKE_EXCEPTION(e);
         block.result_K->fillException(e);
         throw e;
      }

      block.event = block.result_labels->read();
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






void KMeans::runReadBlock(Block& block)
{
   // check if read is complete and we are next in line
   if ( block.event.isDone() && _nextVector == block.vector )
   {
      // make sure opencl event worked
      if ( !block.event )
      {
         E_MAKE_EXCEPTION(e);
         block.event.fillException(e);
         throw e;
      }

      // save each valid clustering result to cluster matrix
      int index {0};
      while ( block.vector.geneX() < _output->geneSize() && index < _kernelSize )
      {
         // read results
         int N = _input->getSampleSize();
         int bestK = (*block.result_K)[index];
         int *bestLabels = &(*block.result_labels)[index * N];

         // save cluster pair if multiple clusters are found
         if ( bestK > 1 )
         {
            CCMatrix::Pair pair(_output);
            pair.addCluster(bestK);

            for ( int i = 0; i < N; ++i )
            {
               for ( int k = 0; k < bestK; ++k )
               {
                  pair.at(k, i) = (k == bestLabels[i]);
               }
            }

            pair.write(block.vector);
         }

         // increment indices
         ++block.vector;
         ++index;
         ++_pairsComplete;
      }

      // update next vector and change block state to start
      _nextVector = block.vector;
      block.state = Block::Start;
   }
}
