#include "exportexpressionmatrix.h"
#include "datafactory.h"






EAbstractAnalytic::ArgumentType ExportExpressionMatrix::getArgumentData(int argument)
{
   // use type declaration
   using Type = EAbstractAnalytic::ArgumentType;

   // figure out which argument is being queried and return its type
   switch (argument)
   {
   case InputData: return Type::DataIn;
   case OutputFile: return Type::FileOut;
   case NoSampleToken: return Type::String;
   default: return Type::Bool;
   }
}






QVariant ExportExpressionMatrix::getArgumentData(int argument, Role role)
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
      case OutputFile: return QString("output");
      case NoSampleToken: return QString("nan");
      default: return QString();
      }
   case Role::Title:
      // figure out which argument is being queried and return title
      switch (argument)
      {
      case InputData: return tr("Input:");
      case OutputFile: return tr("Output:");
      case NoSampleToken: return tr("No Sample Token:");
      default: return QString();
      }
   case Role::WhatsThis:
      // figure out which argument is being queried and return "What's This?" text
      switch (argument)
      {
      case InputData: return tr("Input expression matrix containing expression data.");
      case OutputFile: return tr("Raw output text file that will contain space/tab divided gene expression data.");
      case NoSampleToken: return tr("Expected token for expressions that have no value.");
      default: return QString();
      }
   case Role::DataType:
      // if this is input data argument return data type else return nothing
      switch (argument)
      {
      case InputData: return DataFactory::ExpressionMatrixType;
      default: return QString();
      }
   case Role::FileFilters:
      // if this is output file argument return file filter else return nothing
      switch (argument)
      {
      case OutputFile: return tr("Raw text file %1").arg("(*.txt)");
      default: return QString();
      }
   default:
      return QVariant();
   }
}






void ExportExpressionMatrix::setArgument(int argument, QVariant value)
{
   // figure out which argument is being set and set it
   switch (argument)
   {
   case NoSampleToken:
      _noSampleToken = value.toString();
      break;
   }
}






void ExportExpressionMatrix::setArgument(int argument, EAbstractData* data)
{
   // if argument is input data set it
   if ( argument == InputData )
   {
      _input = dynamic_cast<ExpressionMatrix*>(data);
   }
}






void ExportExpressionMatrix::setArgument(int argument, QFile* file)
{
   // if argument is output file set it
   if ( argument == OutputFile )
   {
      _output = file;
   }
}






bool ExportExpressionMatrix::initialize()
{
   // make sure input and output arguments were set
   if ( !_input || !_output )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Argument Error"));
      e.setDetails(QObject::tr("Did not get valid input and/or output arguments."));
      throw e;
   }

   // do not have data pre-allocate
   return false;
}






void ExportExpressionMatrix::runSerial()
{
   // use expression declaration
   using Expression = ExpressionMatrix::Expression;
   using Transform = ExpressionMatrix::Transform;

   // initialize variable used to track percent complete
   int lastPercent {0};

   // get gene names, sample names, and transform
   const QList<EMetadata *> *geneNames = _input->getGeneNames().toArray();
   const QList<EMetadata *> *sampleNames = _input->getSampleNames().toArray();
   Transform transform = _input->getTransform();

   // create text stream to output file and write until end reached
   QTextStream stream(_output);
   stream.setRealNumberPrecision(12);

   for ( int i = 0; i < _input->getSampleSize(); i++ )
   {
      stream << *sampleNames->at(i)->toString() << "\t";
   }
   stream << "\n";

   for ( int i = 0; i < _input->getGeneSize(); i++ )
   {
      // if interruption is requested exit immediately
      if ( isInterruptionRequested() )
      {
         return;
      }

      // write a single row to a line in the output file
      stream << *geneNames->at(i)->toString();

      for ( int j = 0; j < _input->getSampleSize(); j++ )
      {
         Expression value = _input->data(_input->index(i, j), Qt::DisplayRole).toFloat();

         // if value is NAN use the no sample token
         if ( std::isnan(value) )
         {
            stream << "\t" << _noSampleToken;
         }

         // else this is a normal floating point expression
         else
         {
            // figure out which transform is being used and reverse it before writing
            switch (transform)
            {
            case Transform::None:
               break;
            case Transform::NLog:
               value = exp(value);
               break;
            case Transform::Log2:
               value = pow(2, value);
               break;
            case Transform::Log10:
               value = pow(10, value);
               break;
            }

            stream << "\t" << value;
         }
      }

      stream << "\n";

      // figure out percent complete
      int newPercent = 100 * (i + 1) / _input->getGeneSize();

      // if percent complete has changed update it and emit progressed
      if ( newPercent != lastPercent )
      {
         lastPercent = newPercent;
         emit progressed(lastPercent);
      }
   }

   // make sure writing output file worked
   if ( stream.status() != QTextStream::Ok )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Qt Text Stream encountered an unknown error."));
      throw e;
   }
}
