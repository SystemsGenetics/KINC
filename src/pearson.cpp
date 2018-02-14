#include <ace/core/metadata.h>

#include "pearson.h"
#include "datafactory.h"



using namespace std;






Pearson::~Pearson()
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






EAbstractAnalytic::ArgumentType Pearson::getArgumentData(int argument)
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
   case MinThreshold: return Type::Double;
   case MaxThreshold: return Type::Double;
   case BlockSize: return Type::Integer;
   case KernelSize: return Type::Integer;
   default: return Type::Bool;
   }
}






QVariant Pearson::getArgumentData(int argument, Role role)
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
      case MinSamples: return QString("min");
      case MinThreshold: return QString("minthresh");
      case MaxThreshold: return QString("maxthresh");
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
      case MinThreshold: return tr("Minimum Threshold:");
      case MaxThreshold: return tr("Maximum Threshold:");
      case BlockSize: return tr("Block Size:");
      case KernelSize: return tr("Kernel Size:");
      default: return QVariant();
      }
   case Role::WhatsThis:
      // figure out which argument is being queried and return "What's This?" text
      switch (argument)
      {
      case InputData: return tr("Input expression matrix that will be used to compute pearson"
                                " coefficients.");
      case ClusterData: return tr("Cluster matrix to compute correlations for separate clusters.");
      case OutputData: return tr("Output correlation matrixx that will store pearson coefficient"
                                 " results.");
      case MinSamples: return tr("Minimum size of samples two genes must share to generate a spearman"
                              " coefficient.");
      case MinThreshold: return tr("Minimum threshold that a correlation value must be equal to or"
                                   " greater than to be added to the correlation matrix.");
      case MaxThreshold: return tr("Maximum threshold that a correlation value must be equal to or"
                                   " lesser than to be added to the correlation matrix.");
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
      case MinThreshold: return 0.5;
      case MaxThreshold: return 1.0;
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
      case MinThreshold: return -1.0;
      case MaxThreshold: return -1.0;
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
      case MinThreshold: return 1.0;
      case MaxThreshold: return 1.0;
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
      case ClusterData: return DataFactory::CCMatrixType;
      case OutputData: return DataFactory::CorrelationMatrixType;
      default: return QVariant();
      }
   default:
      return QVariant();
   }
}






void Pearson::setArgument(int argument, QVariant value)
{
   // figure out which argument is being set and set it
   switch (argument)
   {
   case MinSamples:
      _minSamples = value.toInt();
      break;
   case MinThreshold:
      _minThreshold = value.toDouble();
      break;
   case MaxThreshold:
      _maxThreshold = value.toDouble();
      break;
   case BlockSize:
      _blockSize = value.toInt();
      break;
   case KernelSize:
      _kernelSize = value.toInt();
      break;
   }
}






void Pearson::setArgument(int argument, EAbstractData *data)
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






bool Pearson::initialize()
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
   *(name->toString()) = "pearson";
   correlations.toArray()->append(name);
   _output->initialize(_input->getGeneNames(),correlations);
   return false;
}






int Pearson::fetchData(const GenePair::Vector& vector, const CCMatrix::Pair& pair, int k, float *x, float *y)
{
   // read in gene expressions
   ExpressionMatrix::Gene gene1(_input);
   ExpressionMatrix::Gene gene2(_input);

   gene1.read(vector.geneX());
   gene2.read(vector.geneY());

   // populate x and y with shared expressions of gene pair
   int numSamples = 0;

   if ( pair.clusterSize() > 0 )
   {
      // add samples that are in the cluster
      for ( int i = 0; i < _input->getSampleSize(); ++i )
      {
         if ( pair.at(k, i) == 1 )
         {
            x[numSamples] = gene1.at(i);
            y[numSamples] = gene2.at(i);
            ++numSamples;
         }
      }
   }
   else
   {
      // add samples that are valid
      for ( int i = 0; i < _input->getSampleSize(); ++i )
      {
         if ( !isnan(gene1.at(i)) && !isnan(gene2.at(i)) )
         {
            x[numSamples] = gene1.at(i);
            y[numSamples] = gene2.at(i);
            ++numSamples;
         }
      }
   }

   return numSamples;
}






void Pearson::runSerial()
{
   // initialize percent complete and steps
   int lastPercent {0};
   qint64 steps {0};
   qint64 totalSteps {_output->geneSize()*(_output->geneSize() - 1)/2};

   // initialize work arrays
   float x[_input->getSampleSize()];
   float y[_input->getSampleSize()];

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

      // fetch a and b arrays from expression matrix
      int n = fetchData(vector, inPair, cluster, x, y);

      // compute correlation only if there are enough samples
      float result = NAN;

      if ( n >= _minSamples )
      {
         // compute intermediate sums
         float sumx = 0;
         float sumy = 0;
         float sumx2 = 0;
         float sumy2 = 0;
         float sumxy = 0;

         for ( int i = 0; i < n; ++i )
         {
            sumx += x[i];
            sumy += y[i];
            sumx2 += x[i] * x[i];
            sumy2 += y[i] * y[i];
            sumxy += x[i] * y[i];
         }

         // compute Pearson correlation coefficient
         result = (n*sumxy - sumx*sumy) / sqrt((n*sumx2 - sumx*sumx) * (n*sumy2 - sumy*sumy));
      }

      // save correlation if within threshold limits
      if ( !isnan(result) && _minThreshold <= result && result <= _maxThreshold )
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






int Pearson::getBlockSize()
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






bool Pearson::runBlock(int index)
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






void Pearson::initializeKernel()
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
   _program->addFile(":/opencl/pearson.cl");
   if ( !_program->compile() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("OpenCL Compile Error"));
      e.setDetails(tr("OpenCL program failed to compile:\n\n%1").arg(_program->getBuildError())
                   .replace('\n',"<br/>"));
      throw e;
   }

   // get kernel from compiled code making sure it worked
   _kernel = _program->makeKernel("calculatePearsonBlock").release();
   if ( !*_program )
   {
      E_MAKE_EXCEPTION(e);
      _program->fillException(e);
      throw e;
   }
}






void Pearson::initializeBlockExpressions()
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






void Pearson::initializeKernelArguments()
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
   for (int i = 0; i < _blockSize ;++i)
   {
      _blocks[i] = new Block(device,_input->getSampleSize(),_kernelSize);
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
   _kernel->setArgument(4, (cl_int)_minSamples);
   _kernel->setDimensionCount(1);
   _kernel->setGlobalSize(0,_kernelSize);
   _kernel->setWorkgroupSize(0,workgroupSize);

   // make sure everything with kernel worked
   if ( !*_kernel )
   {
      E_MAKE_EXCEPTION(e);
      _kernel->fillException(e);
      throw e;
   }
}






void Pearson::runStartBlock(Block& block)
{
   // check if there are more gene comparisons to compute
   if ( _vector.geneX() < _output->geneSize() )
   {
      // set block xy to beginning of read comparisons
      block.vector = _vector;
      block.cluster = _cluster;

      // load first input pair
      if ( _cluster != 0 )
      {
         _inPair.read(_vector);
      }

      // copy list of comparisons to be done on this run
      int index {0};
      while ( _vector.geneX() < _output->geneSize() && index < _kernelSize )
      {
         // copy gene pair
         (*block.pairs)[index] = { _vector.geneX(), _vector.geneY() };

         // read next cluster pair
         if ( _cluster == 0 )
         {
            _inPair.read(_vector);
         }

         // copy sample mask if there is one
         int N = _input->getSampleSize();
         cl_char *sampleMask = &(*block.sampleMasks)[index * N];

         if ( _inPair.clusterSize() > 0 )
         {
            for ( int i = 0; i < N; ++i )
            {
               sampleMask[i] = _inPair.at(_cluster, i);
            }
         }
         else
         {
            sampleMask[0] = -1;
         }

         // increment to next cluster
         _cluster = (_cluster + 1) % max(1, _inPair.clusterSize());

         if ( _cluster == 0 )
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

      // write comparison references to device making sure it worked
      block.event = block.pairs->write();
      if ( !*block.pairs )
      {
         E_MAKE_EXCEPTION(e);
         block.pairs->fillException(e);
         throw e;
      }

      // write sample masks to device making sure it worked
      block.event = block.sampleMasks->write();
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






void Pearson::runLoadBlock(Block& block)
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
      _kernel->setBuffer(3, block.sampleMasks);
      _kernel->setBuffer(5, block.workBuffer);
      _kernel->setBuffer(6, block.results);
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






void Pearson::runExecuteBlock(Block& block)
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

      // read spearman answers from device and make sure it worked
      block.event = block.results->read();
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






void Pearson::runReadBlock(Block& block)
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

      // load first input/output pair
      if ( block.cluster != 0 )
      {
         _inPair.read(block.vector);
         _outPair.read(block.vector);
      }

      // save each valid correlation result to correlation matrix
      int index {0};
      while ( block.vector.geneX() < _output->geneSize() && index < _kernelSize )
      {
         // read next cluster pair
         if ( block.cluster == 0 )
         {
            _inPair.read(block.vector);
            _outPair.clearClusters();
            _outPair.addCluster(max(1, _inPair.clusterSize()));
         }

         // read result
         float result = (*block.results)[index];

         // save correlation if within threshold limits
         if ( !isnan(result) && _minThreshold <= result && result <= _maxThreshold )
         {
            _outPair.at(block.cluster, 0) = result;

            if ( block.cluster == _outPair.clusterSize() - 1 )
            {
               _outPair.write(block.vector);
            }
         }

         // increment to next cluster
         block.cluster = (block.cluster + 1) % _outPair.clusterSize();

         if ( block.cluster == 0 )
         {
            ++block.vector;
         }

         ++index;
         ++_pairsComplete;
      }

      // save last output pair
      if ( block.cluster != 0 )
      {
         _outPair.write(block.vector);
      }

      // update next vector and change block state to start
      _nextVector = block.vector;
      block.state = Block::Start;
   }
}
