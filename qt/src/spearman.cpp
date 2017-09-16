#include <ace/core/metadata.h>

#include "spearman.h"
#include "datafactory.h"
#include "expressionmatrix.h"
#include "correlationmatrix.h"



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
   case OutputData: return Type::DataOut;
   case Minimum: return Type::Integer;
   case BlockSize: return Type::Integer;
   case KernelSize: return Type::Integer;
   case PreAllocate: return Type::Bool;
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
      case OutputData: return QString("output");
      case Minimum: return QString("min");
      case BlockSize: return QString("bsize");
      case KernelSize: return QString("ksize");
      case PreAllocate: return QString("prealloc");
      default: return QVariant();
      }
   case Role::Title:
      // figure out which argument is being queried and return title
      switch (argument)
      {
      case InputData: return tr("Input:");
      case OutputData: return tr("Output:");
      case Minimum: return tr("Minimum Sample Size:");
      case BlockSize: return tr("Block Size:");
      case KernelSize: return tr("Kernel Size:");
      case PreAllocate: return tr("Pre-Allocate Output?");
      default: return QVariant();
      }
   case Role::WhatsThis:
      // figure out which argument is being queried and return "What's This?" text
      switch (argument)
      {
      case InputData: return tr("Input expression matrix that will be used to compute spearman"
                                " coefficients.");
      case OutputData: return tr("Output correlation matrixx that will store spearman coefficient"
                                 " results.");
      case Minimum: return tr("Minimum size of samples two genes must share to generate a spearman"
                              " coefficient.");
      case BlockSize: return tr("This option only applies if OpenCL is used. Total number of blocks"
                                " to run for execution.");
      case KernelSize: return tr("This option only applies if OpenCL is used. Total number of"
                                 " kernels to run per block of execution.");
      case PreAllocate: return tr("Should the output correlation matrix have file space"
                                  " pre-allocated? WARNING this only works in linux systems.");
      default: return QVariant();
      }
   case Role::DefaultValue:
      // figure out which argument is being queried and if applicable return default value else
      // return nothing
      switch (argument)
      {
      case Minimum: return 30;
      case BlockSize: return 4;
      case KernelSize: return 4096;
      default: return QVariant();
      }
   case Role::Minimum:
      // figure out which argument is being queried and if applicable return minimum value else
      // return nothing
      switch (argument)
      {
      case Minimum: return 1;
      case BlockSize: return 1;
      case KernelSize: return 1;
      default: return QVariant();
      }
   case Role::Maximum:
      // figure out which argument is being queried and if applicable return maximum value else
      // return nothing
      switch (argument)
      {
      case Minimum: return INT_MAX;
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
   case Minimum:
      _minimum = value.toInt();
      break;
   case BlockSize:
      _blockSize = value.toInt();
      break;
   case KernelSize:
      _kernelSize = value.toInt();
      break;
   case PreAllocate:
      _allocate = value.toBool();
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
   case OutputData:
      _output = dynamic_cast<CorrelationMatrix*>(data);
      break;
   }
}






bool Spearman::initialize()
{
   // make sure there is valid input and output
   if ( !_input || !_output )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Argument Error"));
      e.setDetails(QObject::tr("Did not get valid input and/or output arguments."));
      throw e;
   }

   // initialize new correlation matrix and return pre-allocation argument
   _output->initialize(_input->getGeneNames(),_input->getSampleSize(),1,1);
   return _allocate;
}






void Spearman::runSerial()
{
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


   // calculate total number of calculations that will be done and return block size
   qint64 geneSize {_output->getGeneSize()};
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

   // get first power of 2 number greater or equal to sample size
   int pow2Size {2};
   while ( pow2Size < _output->getSampleSize() )
   {
      pow2Size *= 2;
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
      _blocks[i] = new Block(device,pow2Size,_kernelSize);
   }

   // figure out workgroup size for opencl kernels
   int workgroupSize {2};
   while ( workgroupSize*2 <= _kernelSize && workgroupSize < (int)_kernel->getMaxWorkgroupSize() )
   {
      workgroupSize *= 2;
   }

   // set all static kernel arguments
   _kernel->setArgument(0,(cl_int)_output->getSampleSize());
   _kernel->setArgument(1,(cl_int)pow2Size);
   _kernel->setArgument(2,(cl_int)_minimum);
   _kernel->setBuffer(4,_expressions);
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






void Spearman::runStartBlock(Block& block)
{
   // check if there are more gene comparisons to compute
   if ( _x < _output->getGeneSize() )
   {
      // set block xy to beginning of read comparisons
      block.x = _x;
      block.y = _y;

      // copy list of comparisons to be done on this run
      int index {0};
      while ( _x < _output->getGeneSize() && index < _kernelSize )
      {
         (*block.references)[index*2] = _x;
         (*block.references)[(index*2)+1] = _y;
         CorrelationMatrix::increment(_x,_y);
         ++index;
      }

      // copy any remaining and unused comparisons to zero
      while ( index < _kernelSize )
      {
         (*block.references)[index*2] = 0;
         (*block.references)[(index*2)+1] = 0;
         ++index;
      }

      // write comparison references to device making sure it worked
      block.event = block.references->write();
      if ( !*block.references )
      {
         E_MAKE_EXCEPTION(e);
         block.references->fillException(e);
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
      _kernel->setBuffer(3,block.references);
      _kernel->setBuffer(5,block.workBuffer);
      _kernel->setBuffer(6,block.rankBuffer);
      _kernel->setBuffer(7,block.answers);
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






void Spearman::runExecuteBlock(Block& block)
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
      block.event = block.answers->read();
      if ( !*block.answers )
      {
         E_MAKE_EXCEPTION(e);
         block.answers->fillException(e);
         throw e;
      }

      // change block state to read
      block.state = Block::Read;
   }
}






void Spearman::runReadBlock(Block& block)
{
   // use pair declaration and check if read is complete
   using Pair = CorrelationMatrix::Pair;
   if ( block.event.isDone() )
   {
      // make sure opencl event worked
      if ( !block.event )
      {
         E_MAKE_EXCEPTION(e);
         block.event.fillException(e);
         throw e;
      }

      // make new correlation pair
      Pair pair(_output);

      // initialize mode size to one and set mask to all ones
      pair.setModeSize(1);
      for (int i = 0; i < _output->getSampleSize() ;++i)
      {
         pair.mode(0,i) = 1;
      }

      // iterate through all valid spearman answers and write each one to correlation matrix
      int index {0};
      while ( block.x < _output->getGeneSize() && index < _kernelSize )
      {
         pair.at(0,0) = (*block.answers)[index];
         pair.write(block.x,block.y);
         CorrelationMatrix::increment(block.x,block.y);
         ++index;
         ++_pairsComplete;
      }

      // change block state to start
      block.state = Block::Start;
   }
}
