#include "importexpressionmatrix.h"
#include "datafactory.h"



const char* ImportExpressionMatrix::NONE {QT_TR_NOOP("none")};
const char* ImportExpressionMatrix::NLOG {QT_TR_NOOP("natural logarithm")};
const char* ImportExpressionMatrix::LOG2 {QT_TR_NOOP("logarithm base 2")};
const char* ImportExpressionMatrix::LOG10 {QT_TR_NOOP("logarithm base 10")};






EAbstractAnalytic::ArgumentType ImportExpressionMatrix::getArgumentData(int argument)
{
   // use type declaration
   using Type = EAbstractAnalytic::ArgumentType;

   // figure out which argument is being queried and return its type
   switch (argument)
   {
   case InputFile: return Type::FileIn;
   case OutputData: return Type::DataOut;
   case NoSampleToken: return Type::String;
   case SampleSize: return Type::Integer;
   case TransformArg: return Type::Combo;
   default: return Type::Bool;
   }
}






QVariant ImportExpressionMatrix::getArgumentData(int argument, EAbstractAnalytic::Role role)
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
      case OutputData: return QString("output");
      case NoSampleToken: return QString("nan");
      case SampleSize: return QString("size");
      case TransformArg: return QString("transform");
      default: return QString();
      }
   case Role::Title:
      // figure out which argument is being queried and return title
      switch (argument)
      {
      case InputFile: return tr("Input:");
      case OutputData: return tr("Output:");
      case NoSampleToken: return tr("No Sample Token:");
      case SampleSize: return tr("Sample Size:");
      case TransformArg: return tr("Transform:");
      default: return QString();
      }
   case Role::WhatsThis:
      // figure out which argument is being queried and return "What's This?" text
      switch (argument)
      {
      case InputFile: return tr("Raw input text file containing space/tab divided gene expression"
                                " data.");
      case OutputData: return tr("Output expression matrix that will contain expression data.");
      case NoSampleToken: return tr("Expected token for expressions that have no value.");
      case SampleSize: return tr("Total number of samples per gene. 0 indicates the text file "
                                 "contains a header of sample names to be read to determine size.");
      case TransformArg: return tr("What type of transformation, if any, should be done to the raw"
                                   " epxression values when imported to the expression matrix.");
      default: return QString();
      }
   case Role::ComboValues:
      // if this is transform argument return combo values else return nothing
      switch (argument)
      {
      case TransformArg: return QStringList() << tr(NONE) << tr(NLOG) << tr(LOG2) << tr(LOG10);
      default: return QStringList();
      }
   case Role::Minimum:
      // if this is sample size argument return minimum else return nothing
      switch (argument)
      {
      case SampleSize: return 0;
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
      case OutputData: return DataFactory::ExpressionMatrixType;
      default: return QString();
      }
   default:
      return QVariant();
   }
}






void ImportExpressionMatrix::setArgument(int argument, QVariant value)
{
   // figure out which argument is being set and set it
   switch (argument)
   {
   case SampleSize:
      _sampleSize = value.toInt();
      break;
   case NoSampleToken:
      _noSampleToken = value.toString();
      break;
   case TransformArg:
      {
         // get option and map it to correct transform type
         const QString option = value.toString();
         if ( option == tr(NONE) )
         {
            _transform = Transform::None;
         }
         else if ( option == tr(NLOG) )
         {
            _transform = Transform::NLog;
         }
         else if ( option == tr(LOG2) )
         {
            _transform = Transform::Log2;
         }
         else if ( option == tr(LOG10) )
         {
            _transform = Transform::Log10;
         }
      }
      break;
   }
}






void ImportExpressionMatrix::setArgument(int argument, QFile* file)
{
   // if argument is input file set it
   if ( argument == InputFile )
   {
      _input = file;
   }
}






void ImportExpressionMatrix::setArgument(int argument, EAbstractData* data)
{
   // if argument is output data set it
   if ( argument == OutputData )
   {
      _output = dynamic_cast<ExpressionMatrix*>(data);
   }
}






bool ImportExpressionMatrix::initialize()
{
   // make sure input and output arguments were set
   if ( !_input || !_output )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Argument Error"));
      e.setDetails(QObject::tr("Did not get valid input and/or output arguments."));
      throw e;
   }

   // make sure sample size is not negative
   if ( _sampleSize < 0 )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Argument Error"));
      e.setDetails(QObject::tr("Sample size must be zero or positive."));
      throw e;
   }

   // do not have data pre-allocate
   return false;
}






bool ImportExpressionMatrix::runBlock(int block)
{
   // do not use block and use expression and gene declarations
   Q_UNUSED(block);
   using Expression = ExpressionMatrix::Expression;
   using EGene = ExpressionMatrix::Gene;

   // custom single linked list structure to use for building list of gene expressions
   struct Gene
   {
      Gene(int size)
      {
         expressions = new Expression[size];
      }
      ~Gene()
      {
         delete next;
         delete[] expressions;
      }

      Expression* expressions;
      Gene* next {nullptr};
   };

   // initialize variable used to track percent complete
   int lastPercent {0};

   // initialize gene expression linked list
   std::unique_ptr<Gene> geneRoot;
   Gene* geneTail {nullptr};

   // initialize gene and sample name lists
   QStringList geneNames;
   QStringList sampleNames;

   // if sample size is not zero then build sample name list
   if ( _sampleSize != 0 )
   {
      for (int i = 0; i < _sampleSize ;++i)
      {
         sampleNames.append(QString::number(i));
      }
   }

   // create text stream from input file and read until end reached
   QTextStream stream(_input);
   while ( !stream.atEnd() )
   {
      // read a single line from stream and split it into words using any whitespace characters as
      // delimiters
      QString line = stream.readLine();
      auto words = line.splitRef(QRegExp("\\s+"),QString::SkipEmptyParts);

      // check to see if there are any words from read line
      if ( words.size() > 0 )
      {
         // check to see if sample size has not yet been built
         if ( _sampleSize == 0 )
         {
            // set sample size and build sample name list
            _sampleSize = words.size();
            for (int i = 0; i < words.size() ;++i)
            {
               sampleNames.append(words.at(i).toString());
            }
         }

         // else this is a normal gene expression line to read in
         else
         {
            // make sure the number of words matches expected sample size
            if ( words.size() != (_sampleSize+1) )
            {
               E_MAKE_EXCEPTION(e);
               e.setTitle(tr("Parsing Error"));
               e.setDetails(tr("Encountered gene expression line with incorrect amount of fields. "
                               "Read in %1 fields when it should have been %2. Gene name is %3.")
                            .arg(words.size()-1).arg(_sampleSize).arg(words.at(0).toString()));
               throw e;
            }

            // make new linked list gene node
            Gene* gene {new Gene(_sampleSize)};

            // if linked list is empty make new node head of list
            if ( !geneTail )
            {
               geneTail = gene;
               geneRoot.reset(geneTail);
            }

            // else add to last node in last and make new node the tail
            else
            {
               geneTail->next = gene;
               geneTail = gene;
            }

            // add gene name to list
            geneNames.append(words.at(0).toString());

            // iterate through all gene expressions
            for (int i = 1; i < words.size() ;++i)
            {
               // if word matches no sample token string set it as such
               if ( words.at(i) == _noSampleToken )
               {
                  gene->expressions[i-1] = NAN;
               }

               // else this is a normal floating point expression
               else
               {
                  // read in the floating point value
                  bool ok;
                  Expression value = words.at(i).toDouble(&ok);

                  // make sure reading worked
                  if ( !ok )
                  {
                     E_MAKE_EXCEPTION(e);
                     e.setTitle(tr("Parsing Error"));
                     e.setDetails(tr("Failed reading in expression value \"%1\" for gene %2.")
                                  .arg(words.at(i).toString()).arg(words.at(0).toString()));
                     throw e;
                  }

                  // figure out which transform is being used and do it to expression value adding
                  // it to gene node's list of expressions
                  switch (_transform)
                  {
                  case Transform::None:
                     gene->expressions[i-1] = value;
                     break;
                  case Transform::NLog:
                     gene->expressions[i-1] = log(value);
                     break;
                  case Transform::Log2:
                     gene->expressions[i-1] = log2(value);
                     break;
                  case Transform::Log10:
                     gene->expressions[i-1] = log10(value);
                     break;
                  }
               }
            }
         }
      }

      // figure out percent complete
      int newPercent = _input->pos()*50/_input->size();

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

   // reset gene tail to beginning of linked list and initialize output data
   geneTail = geneRoot.get();
   _output->initialize(geneNames,sampleNames);

   // create gene iterator for expression matrix and iterate through all genes
   EGene gene(_output);
   for (int i = 0; i < _output->getGeneSize() ;++i)
   {
      // make sure next gene tail is valid
      if ( !geneTail )
      {
         E_MAKE_EXCEPTION(e);
         e.setTitle(tr("Internal Error"));
         e.setDetails(tr("Prematurely reached end of linked gene list."));
         throw e;
      }

      // iterate through all gene expressions and copy them to gene iterator
      for (int x = 0; x < _output->getSampleSize() ;++x)
      {
         gene[x] = geneTail->expressions[x];
      }

      // write gene iterator to expression matrix
      gene.write(i);

      // move to next gene node in linked list and calculate percent complete
      geneTail = geneTail->next;
      int newPercent = 50 + i*50/_output->getGeneSize();

      // if precent complete has changed update it and emit progressed
      if ( newPercent != lastPercent )
      {
         lastPercent = newPercent;
         emit progressed(lastPercent);
      }
   }

   // set transform used in expression matrix and return block finished
   _output->setTransform(_transform);
   return false;
}
