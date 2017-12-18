#include <ace/core/metadata.h>
#include <gsl/gsl_statistics.h>

#include "pearson.h"
#include "datafactory.h"
#include "expressionmatrix.h"
#include "correlationmatrix.h"



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
   case OutputData: return Type::DataOut;
   case Minimum: return Type::Integer;
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
      case OutputData: return QString("output");
      case Minimum: return QString("min");
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
      case OutputData: return tr("Output:");
      case Minimum: return tr("Minimum Sample Size:");
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
      case OutputData: return tr("Output correlation matrixx that will store pearson coefficient"
                                 " results.");
      case Minimum: return tr("Minimum size of samples two genes must share to generate a spearman"
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
      case Minimum: return 30;
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
      case Minimum: return 1;
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
      case Minimum: return INT_MAX;
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
   case Minimum:
      _minimum = value.toInt();
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
   case OutputData:
      _output = dynamic_cast<CorrelationMatrix*>(data);
      break;
   }
}






bool Pearson::initialize()
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
   if ( _minimum < 1 )
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


   // calculate total number of calculations that will be done and return block size
   qint64 geneSize {_output->geneSize()};
   _totalPairs = geneSize*(geneSize-1)/2;
   return _blockSize;
}





/////////////////////////////////////
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
   // TODO: Connor never changed the kernel function name!!??
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


void Pearson::runSerial()
{
    // Initialize correlation gene pair and expression genes for input/output
    CorrelationMatrix::Pair pair(_output);
    pair.addCluster();
    ExpressionMatrix::Gene gene1(_input);
    ExpressionMatrix::Gene gene2(_input);

    // Initialize the first gene comparision
    GenePair::Vector vector;

    // Iterate through the gene pairs
    // _output->getGeneSize() represents all of the possible comparisons that can
    // be done limited by the gene with the least amount of samples
    while(vector.geneX() < _output->geneSize())
    {
        // Initialize the running sum variables
        double sumx {0.0}, sumy {0.0}, sumx2 {0.0}, sumy2 {0.0}, sumxy {0.0};
        // A trap to kill this function if there is an interrupt i.e. the user wants to kill
        // the process
        if ( isInterruptionRequested() )
        {
           return;
        }

        // Initialize sample size and read in gene expressions
        // Keep track of the count of viable samples that ARE NOT nans
        int size {0};
        gene1.read(vector.geneX());
        gene2.read(vector.geneY());

        // Iterate through the gene1 and gene2 sample row
        for (auto i = 0; i < _input->getSampleSize(); i++)
        {
            if(!std::isnan(gene1.at(i)) && !std::isnan(gene2.at(i)))
            {
                // Multiply the current scores together
                //  - keep a running sum of this
                sumxy += gene1.at(i) * gene2.at(i);
                //  - keep a running sum the x and y values
                sumx += gene1.at(i);
                sumy += gene2.at(i);
                //  - keep a running sum of the square of x and y values
                sumx2 += pow(gene1.at(i), 2.0);
                sumy2 += pow(gene2.at(i), 2.0);
                size++;
            }
        }

        // Once we're out of this loop we have the summation values needed to
        // calculate the Pearson Coefficient for the selected gene comparison
        double numerator {(size * sumxy) - (sumx * sumy)};
        double denominator {sqrt(((size*sumx2) - pow(sumx, 2.0)))*sqrt(((size*sumy2) - pow(sumy, 2.0)))};

        // Set result to correlation
        pair.at(0,0) = numerator / denominator;

        // Save this gene pair and increment to the next comparison
        if ( !isnan(pair.at(0,0)) && pair.at(0,0) >= _minThreshold
             && pair.at(0,0) <= _maxThreshold )
        {
         pair.write(vector);
        }
        ++vector;
    }
}



void Pearson::initializeKernelArguments()
{
   // get opencl device
   EOpenCLDevice& device {EOpenCLDevice::getInstance()};

   // get first power of 2 number greater or equal to sample size
   int pow2Size {2};
   while ( pow2Size < _input->getSampleSize() )
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
   int workgroupSize {_kernelSize};
   while ( workgroupSize > (int)_kernel->getMaxWorkgroupSize() )
   {
      workgroupSize /= 2;
   }

   // set all static kernel arguments
   _kernel->setArgument(0,(cl_int)_input->getSampleSize());
   _kernel->setArgument(1,(cl_int)_minimum);
   _kernel->setBuffer(2,_expressions);
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

      // copy list of comparisons to be done on this run
      int index {0};
      while ( _vector.geneX() < _output->geneSize() && index < _kernelSize )
      {
         (*block.references)[index*2] = _vector.geneX();
         (*block.references)[(index*2)+1] = _vector.geneY();
         ++_vector;
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
      _kernel->setBuffer(3,block.workBuffer);
      _kernel->setBuffer(4,block.references);
      _kernel->setBuffer(5,block.answers);
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






void Pearson::runReadBlock(Block& block)
{
   // use pair declaration and check if read is complete and we are next in line
   using Pair = CorrelationMatrix::Pair;
   if ( block.event.isDone() && _nextVector == block.vector )
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
      pair.addCluster();

      // iterate through all valid spearman answers and write each one to correlation matrix that is
      // not a NaN and is within threshold limits
      int index {0};
      while ( block.vector.geneX() < _output->geneSize() && index < _kernelSize )
      {
         pair.at(0,0) = (*block.answers)[index];
         if ( !isnan(pair.at(0,0)) && pair.at(0,0) >= _minThreshold
              && pair.at(0,0) <= _maxThreshold )
         {
            pair.write(block.vector);
         }
         ++(block.vector);
         ++index;
         ++_pairsComplete;
      }

      // update next vector and change block state to start
      _nextVector = block.vector;
      block.state = Block::Start;
   }
}
