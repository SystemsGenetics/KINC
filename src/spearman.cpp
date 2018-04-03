#include <ace/core/metadata.h>

#include "spearman.h"
#include "datafactory.h"



using namespace std;






Spearman::~Spearman()
{
   // check if blocks are allocated
   if ( _blocks )
   {
      // iterate through all blocks and delete them
      for (int i = 0; i < _blockSize ;++i)
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






EAbstractAnalytic::ArgumentType Spearman::getArgumentData(int argument)
{
   // use type declaration
   using Type = EAbstractAnalytic::ArgumentType;

   // figure out which argument is being queried and return its type
   switch (argument)
   {
   case InputData: return Type::DataIn;
   case ClusterData: return Type::DataIn;
   case OutputData: return Type::DataOut;
   case MinSamples: return Type::Integer;
   case MinExpression: return Type::Double;
   case MinCorrelation: return Type::Double;
   case MaxCorrelation: return Type::Double;
   case BlockSize: return Type::Integer;
   case KernelSize: return Type::Integer;
   default: return Type::Bool;
   }
}






QVariant Spearman::getArgumentData(int argument, Role role)
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
      case ClusterData: return QString("cmatrix");
      case OutputData: return QString("output");
      case MinSamples: return QString("minsamp");
      case MinExpression: return QString("minexpr");
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
      case InputData: return tr("Input:");
      case ClusterData: return tr("Cluster Matrix:");
      case OutputData: return tr("Output:");
      case MinSamples: return tr("Minimum Sample Size:");
      case MinExpression: return tr("Minimum Expression:");
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
      case InputData: return tr("Input expression matrix that will be used to compute spearman"
                                " coefficients.");
      case ClusterData: return tr("Cluster matrix to compute correlations for separate clusters.");
      case OutputData: return tr("Output correlation matrixx that will store spearman coefficient"
                                 " results.");
      case MinSamples: return tr("Minimum size of samples two genes must share to generate a spearman"
                              " coefficient.");
      case MinExpression: return tr("Minimum threshold for a gene expression to be included in correlation.");
      case MinCorrelation: return tr("Minimum threshold (absolute value) for a correlation to be saved.");
      case MaxCorrelation: return tr("Maximum threshold (absolute value) for a correlation to be saved.");
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
      case MinExpression: return -INFINITY;
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
      case MinSamples: return 1;
      case MinExpression: return -INFINITY;
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
      case MinSamples: return INT_MAX;
      case MinExpression: return +INFINITY;
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
      case OutputData: return DataFactory::CorrelationMatrixType;
      default: return QVariant();
      }
   default:
      return QVariant();
   }
}






void Spearman::setArgument(int argument, QVariant value)
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






void Spearman::setArgument(int argument, EAbstractData *data)
{
   // figure out which argument is having its data set and if applicable set it
   switch (argument)
   {
   case InputData:
      _input = dynamic_cast<ExpressionMatrix*>(data);
      break;
   case ClusterData:
      _cMatrix = dynamic_cast<CCMatrix*>(data);
      break;
   case OutputData:
      _output = dynamic_cast<CorrelationMatrix*>(data);
      break;
   }
}






bool Spearman::initialize()
{
   // make sure there is valid input and output
   if ( !_input || !_cMatrix || !_output )
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

   // initialize new correlation matrix and return pre-allocation argument
   EMetadata correlations(EMetadata::Array);
   EMetadata* name {new EMetadata(EMetadata::String)};
   *(name->toString()) = "spearman";
   correlations.toArray()->append(name);
   _output->initialize(_input->getGeneNames(),correlations);
   return false;
}






void Spearman::runSerial()
{
   // initialize percent complete and steps
   int lastPercent {0};
   qint64 steps {0};
   qint64 totalSteps {_output->geneSize()*(_output->geneSize() - 1)/2};

   // initialize correlation model
   GenePair::Spearman model;

   // initialize input/output pairs
   CCMatrix::Pair inPair(_cMatrix);
   CorrelationMatrix::Pair outPair(_output);

   // initialize gene pair index
   GenePair::Vector vector;
   int cluster {0};

   // iterate through all gene pairs
   while ( vector.geneX() < _output->geneSize() )
   {
      // make sure interruption is not requested
      if ( isInterruptionRequested() )
      {
         return;
      }

      // read next cluster pair
      if ( cluster == 0 )
      {
         inPair.read(vector);
         outPair.clearClusters();
         outPair.addCluster(max(1, inPair.clusterSize()));
      }

      // compute correlation
      float result = model.compute(_input, vector, inPair, cluster, _minSamples, _minExpression);

      // save correlation if within threshold limits
      if ( !isnan(result) && _minCorrelation <= abs(result) && abs(result) <= _maxCorrelation )
      {
         outPair.at(cluster, 0) = result;

         if ( cluster == outPair.clusterSize() - 1 )
         {
            outPair.write(vector);
         }
      }

      // increment to next cluster
      cluster = (cluster + 1) % outPair.clusterSize();

      if ( cluster == 0 )
      {
         ++vector;
      }

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






int Spearman::getBlockSize()
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

   // initialize input/output pairs
   _inPair = CCMatrix::Pair(_cMatrix);
   _outPair = CorrelationMatrix::Pair(_output);

   // calculate total number of calculations that will be done and return block size
   qint64 geneSize {_output->geneSize()};
   _totalPairs = geneSize*(geneSize-1)/2;
   return _blockSize;
}






bool Spearman::runBlock(int index)
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






void Spearman::initializeKernel()
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
   _program->addFile(":/opencl/spearman.cl");
   if ( !_program->compile() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("OpenCL Compile Error"));
      e.setDetails(tr("OpenCL program failed to compile:\n\n%1").arg(_program->getBuildError())
                   .replace('\n',"<br/>"));
      throw e;
   }

   // get kernel from compiled code making sure it worked
   _kernel = _program->makeKernel("calculateSpearmanBlock").release();
   if ( !*_program )
   {
      E_MAKE_EXCEPTION(e);
      _program->fillException(e);
      throw e;
   }
}






void Spearman::initializeBlockExpressions()
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
   for (int i = 0; i < _input->getRawSize() ;++i)
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






void Spearman::initializeKernelArguments()
{
   // get opencl device
   EOpenCLDevice& device {EOpenCLDevice::getInstance()};

   // get work size (sample size rounded up to next power of 2)
   int workSize {2};
   while ( workSize < _input->getSampleSize() )
   {
      workSize *= 2;
   }

   // increase kernel size to first power of 2 reached
   int pow2 {2};
   while ( pow2 < _kernelSize )
   {
      pow2 *= 2;
   }
   _kernelSize = pow2;

   // initialize blocks
   _blocks = new Block*[_blockSize];
   for (int i = 0; i < _blockSize ;++i)
   {
      _blocks[i] = new Block(device,_input->getSampleSize(),workSize,_kernelSize);
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
   _kernel->setArgument(2, (cl_int)workSize);
   _kernel->setArgument(5, (cl_int)_minSamples);
   _kernel->setArgument(6, (cl_int)_minExpression);
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






void Spearman::runStartBlock(Block& block)
{
   // check if there are more gene comparisons to compute
   if ( _vector.geneX() < _output->geneSize() )
   {
      // set block xy to beginning of read comparisons
      block.vector = _vector;

      // copy gene pairs and sample masks for this block
      int index {0};
      int cluster {0};
      while ( _vector.geneX() < _output->geneSize() && index < _kernelSize )
      {
         // read next cluster pair
         if ( cluster == 0 )
         {
            _inPair.read(_vector);
         }

         // stop if gene pair won't fit in this block
         if ( index + _inPair.clusterSize() > _kernelSize )
         {
            break;
         }

         // copy gene pair
         (*block.pairs)[index] = { _vector.geneX(), _vector.geneY() };

         // copy sample mask if there is one
         int N = _input->getSampleSize();
         cl_char *sampleMask = &(*block.sampleMasks)[index * N];

         if ( _inPair.clusterSize() > 0 )
         {
            for ( int i = 0; i < N; ++i )
            {
               sampleMask[i] = _inPair.at(cluster, i);
            }
         }
         else
         {
            sampleMask[0] = -1;
         }

         // increment to next cluster
         cluster = (cluster + 1) % max(1, _inPair.clusterSize());

         if ( cluster == 0 )
         {
            ++_vector;
         }

         ++index;
      }

      // copy any remaining and unused comparisons to zero
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

      // write sample masks to device making sure it worked
      block.events.append(block.sampleMasks->write());
      if ( !*block.sampleMasks )
      {
         E_MAKE_EXCEPTION(e);
         block.sampleMasks->fillException(e);
         throw e;
      }

      // change block state to load
      block.state = Block::Load;
   }

   // else all comparisons are complete and this block is done
   else
   {
      block.state = Block::Done;
   }
}






void Spearman::runLoadBlock(Block& block)
{
   // check to see if reference loading is complete
   if ( !block.isWaiting() )
   {
      // make sure opencl events worked
      block.checkAllEvents();

      // set kernel arguments and execute it
      _kernel->setBuffer(3, block.pairs);
      _kernel->setBuffer(4, block.sampleMasks);
      _kernel->setBuffer(7, block.workBuffer);
      _kernel->setBuffer(8, block.rankBuffer);
      _kernel->setBuffer(9, block.results);
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






void Spearman::runExecuteBlock(Block& block)
{
   // check to see if kernel execution is complete
   if ( !block.isWaiting() )
   {
      // make sure opencl events worked
      block.checkAllEvents();

      // read spearman results from device and make sure it worked
      block.events.append(block.results->read());
      if ( !*block.results )
      {
         E_MAKE_EXCEPTION(e);
         block.results->fillException(e);
         throw e;
      }

      // change block state to read
      block.state = Block::Read;
   }
}






void Spearman::runReadBlock(Block& block)
{
   // check if read is complete and we are next in line
   if ( !block.isWaiting() && _nextVector == block.vector )
   {
      // make sure opencl events worked
      block.checkAllEvents();

      // save each valid correlation result to correlation matrix
      int index {0};
      int cluster {0};
      while ( block.vector.geneX() < _output->geneSize() && index < _kernelSize )
      {
         // stop if remaining pairs are empty
         if ( (*block.pairs)[index].x == 0 )
         {
            break;
         }

         // read next cluster pair
         if ( cluster == 0 )
         {
            _inPair.read(block.vector);
            _outPair.clearClusters();
            _outPair.addCluster(max(1, _inPair.clusterSize()));
         }

         // read result
         float result = (*block.results)[index];

         // save correlation if within threshold limits
         if ( !isnan(result) && _minCorrelation <= abs(result) && abs(result) <= _maxCorrelation )
         {
            _outPair.at(cluster, 0) = result;

            if ( cluster == _outPair.clusterSize() - 1 )
            {
               _outPair.write(block.vector);
            }
         }

         // increment to next cluster
         cluster = (cluster + 1) % _outPair.clusterSize();

         if ( cluster == 0 )
         {
            ++block.vector;
            ++_pairsComplete;
         }

         ++index;
      }

      // update next vector and change block state to start
      _nextVector = block.vector;
      block.state = Block::Start;

      // clear all opencl events
      block.events.clear();
   }
}
