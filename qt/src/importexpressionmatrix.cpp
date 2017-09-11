#include "importexpressionmatrix.h"
#include "expressionmatrix.h"
#include "datafactory.h"



const char* ImportExpressionMatrix::NONE {QT_TR_NOOP("none")};
const char* ImportExpressionMatrix::NLOG {QT_TR_NOOP("natural logarithm")};
const char* ImportExpressionMatrix::LOG2 {QT_TR_NOOP("logarithm base 2")};
const char* ImportExpressionMatrix::LOG10 {QT_TR_NOOP("logarithm base 10")};






EAbstractAnalytic::ArgumentType ImportExpressionMatrix::getArgumentData(int argument)
{
   using Type = EAbstractAnalytic::ArgumentType;
   switch (argument)
   {
   case InputFile:
      return Type::FileIn;
   case OutputData:
      return Type::DataOut;
   case SampleSize:
      return Type::Integer;
   case TransformArg:
      return Type::Combo;
   default:
      return Type::Bool;
   }
}






QVariant ImportExpressionMatrix::getArgumentData(int argument, EAbstractAnalytic::Role role)
{
   using Role = EAbstractAnalytic::Role;
   switch (role)
   {
   case Role::CommandLineName:
      switch (argument)
      {
      case InputFile: return QString("input");
      case OutputData: return QString("output");
      case SampleSize: return QString("size");
      case TransformArg: return QString("transform");
      default: return QString();
      }
   case Role::Title:
      switch (argument)
      {
      case InputFile: return tr("Input:");
      case OutputData: return tr("Output:");
      case SampleSize: return tr("Sample Size:");
      case TransformArg: return tr("Transform:");
      default: return QString();
      }
   case Role::WhatsThis:
      switch (argument)
      {
      case InputFile: return tr("Raw input text file containing space/tab divided gene expression"
                                " data.");
      case OutputData: return tr("Output expression matrix that will contain expression data.");
      case SampleSize: return tr("Total number of samples per gene. 0 indicates the text file "
                                 "contains a header of sample names to be read to determine size.");
      case TransformArg: return tr("What type of transformation, if any, should be done to the raw"
                                   " epxression values when imported to the expression matrix.");
      default: return QString();
      }
   case Role::ComboValues:
      switch (argument)
      {
      case TransformArg: return QStringList() << tr(NONE) << tr(NLOG) << tr(LOG2) << tr(LOG10);
      default: return QStringList();
      }
   case Role::Minimum:
      switch (argument)
      {
      case SampleSize: return 0;
      default: return QVariant();
      }
   case Role::FileFilters:
      switch (argument)
      {
      case InputFile: tr("Raw text file %1").arg("(*.txt)");
      default: return QString();
      }
   case Role::DataType:
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
   switch (argument)
   {
   case SampleSize: _sampleSize = value.toInt();
   case TransformArg:
      {
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
   }
}






void ImportExpressionMatrix::setArgument(int argument, QFile* file)
{
   if ( argument == InputFile )
   {
      _input = file;
   }
}






void ImportExpressionMatrix::setArgument(int argument, EAbstractData* data)
{
   if ( argument == OutputData )
   {
      _output = dynamic_cast<ExpressionMatrix*>(data);
   }
}






bool ImportExpressionMatrix::initialize()
{
   if ( !_input || !_output )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Argument Error"));
      e.setDetails(QObject::tr("Did not get valid input and/or output arguments."));
      throw e;
   }
   if ( _sampleSize < 0 )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Argument Error"));
      e.setDetails(QObject::tr("Sample size must be zero or positive."));
      throw e;
   }
   return false;
}






bool ImportExpressionMatrix::runBlock(int block)
{
   struct Gene
   {
      Gene(int size)
      {
         expressions = new float[size];
      }
      ~Gene()
      {
         delete next;
         delete[] expressions;
      }

      float* expressions;
      Gene* next {nullptr};
   };
   Q_UNUSED(block);
   std::unique_ptr<Gene> geneRoot;
   Gene* geneTail {nullptr};
   QStringList geneNames;
   QStringList sampleNames;
   if ( _sampleSize != 0 )
   {
      for (int i = 0; i < _sampleSize ;++i)
      {
         sampleNames.append(QString::number(i));
      }
   }
   QTextStream stream(_input);
   while ( !stream.atEnd() )
   {
      QString line = stream.readLine();
      auto words = line.splitRef(QRegExp("\\s+"),QString::SkipEmptyParts);
      if ( words.size() > 0 )
      {
         if ( _sampleSize == 0 )
         {
            _sampleSize = words.size();
            for (int i = 0; i < words.size() ;++i)
            {
               sampleNames.append(words.at(i).toString());
            }
         }
         else
         {
            if ( words.size() != (_sampleSize+1) )
            {
               E_MAKE_EXCEPTION(e);
               e.setTitle(tr("Parsing Error"));
               e.setDetails(tr("Encountered gene expression line with incorrect amount of fields. "
                               "Read in %1 fields when it should have been %2. Gene name is %3.")
                            .arg(words.size()-1).arg(_sampleSize).arg(words.at(0).toString()));
               throw e;
            }
            Gene* gene {new Gene(_sampleSize)};
            if ( !geneTail )
            {
               geneTail = gene;
               geneRoot.reset(geneTail);
            }
            else
            {
               geneTail->next = gene;
               geneTail = gene;
            }
            geneNames.append(words.at(0).toString());
            for (int i = 1; i < words.size() ;++i)
            {
               if ( words.at(i) == _noSampleToken )
               {
                  gene->expressions[i-1] = NAN;
               }
               else
               {
                  bool ok;
                  gene->expressions[i-1] = words.at(i).toFloat(&ok);
                  if ( !ok )
                  {
                     E_MAKE_EXCEPTION(e);
                     e.setTitle(tr("Parsing Error"));
                     e.setDetails(tr("Failed reading in floating point expression value for gene "
                                     "%1.").arg(words.at(0).toString()));
                     throw e;
                  }
               }
            }
         }
      }
   }
   if ( stream.status() != QTextStream::Ok )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Qt Text Stream encountered an unknown error."));
      throw e;
   }
   geneTail = geneRoot.get();
   _output->initialize(geneNames,sampleNames);
   ExpressionMatrix::Gene egene(_output);
   for (int i = 0; i < _output->getGeneSize() ;++i)
   {
      if ( !geneTail )
      {
         E_MAKE_EXCEPTION(e);
         e.setTitle(tr("Internal Error"));
         e.setDetails(tr("Prematurely reached end of linked gene list."));
         throw e;
      }
      for (int x = 0; x < _output->getSampleSize() ;++x)
      {
         egene[x] = geneTail->expressions[x];
      }
      egene.writeGene(i);
      geneTail = geneTail->next;
   }
   return false;
}
