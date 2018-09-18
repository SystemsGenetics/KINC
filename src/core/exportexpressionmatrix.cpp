#include "exportexpressionmatrix.h"
#include "exportexpressionmatrix_input.h"
#include "datafactory.h"
#include "expressionmatrix_gene.h"






int ExportExpressionMatrix::size() const
{
   return 1;
}






void ExportExpressionMatrix::process(const EAbstractAnalytic::Block* result)
{
   Q_UNUSED(result);

   // get gene names, sample names
   EMetaArray geneNames = _input->geneNames().toArray();
   EMetaArray sampleNames = _input->sampleNames().toArray();

   // create text stream to output file
   QTextStream stream(_output);
   stream.setRealNumberPrecision(12);

   // write sample names
   for ( int i = 0; i < _input->sampleSize(); i++ )
   {
      stream << sampleNames.at(i).toString() << "\t";
   }
   stream << "\n";

   // write each gene to a line in output file
   ExpressionMatrix::Gene gene(_input);
   for ( int i = 0; i < _input->geneSize(); i++ )
   {
      // load gene from expression matrix
      gene.read(i);

      // write gene name
      stream << geneNames.at(i).toString();

      // write expression values
      for ( int j = 0; j < _input->sampleSize(); j++ )
      {
         float value {gene.at(j)};

         // if value is NAN use the no sample token
         if ( std::isnan(value) )
         {
            stream << "\t" << _noSampleToken;
         }

         // else this is a normal floating point expression
         else
         {
            stream << "\t" << value;
         }
      }

      stream << "\n";
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






EAbstractAnalytic::Input* ExportExpressionMatrix::makeInput()
{
   return new Input(this);
}






void ExportExpressionMatrix::initialize()
{
   if ( !_input || !_output )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Did not get valid input and/or output arguments."));
      throw e;
   }
}
