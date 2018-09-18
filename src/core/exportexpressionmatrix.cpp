#include "exportexpressionmatrix.h"
#include "exportexpressionmatrix_input.h"
#include "datafactory.h"
#include "expressionmatrix_gene.h"






/*!
 * Return the total number of blocks this analytic must process as steps
 * or blocks of work.
 */
int ExportExpressionMatrix::size() const
{
   return 1;
}






/*!
 * Process the given index with a possible block of results if this analytic
 * produces work blocks. This analytic implementation has no work blocks.
 *
 * @param result
 */
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
            stream << "\t" << _nanToken;
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






/*!
 * Make a new input object and return its pointer.
 */
EAbstractAnalytic::Input* ExportExpressionMatrix::makeInput()
{
   return new Input(this);
}






/*!
 * Initialize this analytic. This implementation checks to make sure the input
 * data object and output file have been set.
 */
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
