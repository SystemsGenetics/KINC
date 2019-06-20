#include "exportexpressionmatrix.h"
#include "exportexpressionmatrix_input.h"
#include "datafactory.h"
#include "expressionmatrix_gene.h"






/*!
 * Return the total number of blocks this analytic must process as steps
 * or blocks of work. This implementation uses a work block for writing the
 * sample names and a work block for writing each gene.
 */
int ExportExpressionMatrix::size() const
{
   EDEBUG_FUNC(this);

   return 1 + _input->geneSize();
}






/*!
 * Process the given index with a possible block of results if this analytic
 * produces work blocks. This implementation uses only the index of the result
 * block to determine which piece of work to do.
 *
 * @param result
 */
void ExportExpressionMatrix::process(const EAbstractAnalyticBlock* result)
{
   EDEBUG_FUNC(this,result);

   // write the sample names in the first step
   if ( result->index() == 0 )
   {
      // get sample names
      EMetaArray sampleNames {_input->sampleNames()};

      // initialize output file stream
      _stream.setDevice(_output);
      _stream.setRealNumberPrecision(_precision);

      // write sample names
      for ( int i = 0; i < _input->sampleSize(); i++ )
      {
         _stream << sampleNames.at(i).toString() << "\t";
      }
      _stream << "\n";
   }

   // write each gene to the output file in a separate step
   else
   {
      // get gene index
      int i = result->index() - 1;

      // get gene name
      QString geneName {_input->geneNames().at(i).toString()};

      // load gene from expression matrix
      ExpressionMatrix::Gene gene(_input);
      gene.read(i);

      // write gene name
      _stream << geneName;

      // write expression values
      for ( int j = 0; j < _input->sampleSize(); j++ )
      {
         float value {gene.at(j)};

         // if value is NAN use the no sample token
         if ( std::isnan(value) )
         {
            _stream << "\t" << _nanToken;
         }

         // else this is a normal floating point expression
         else
         {
            _stream << "\t" << value;
         }
      }

      _stream << "\n";
   }

   // make sure writing output file worked
   if ( _stream.status() != QTextStream::Ok )
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
EAbstractAnalyticInput* ExportExpressionMatrix::makeInput()
{
   EDEBUG_FUNC(this);

   return new Input(this);
}






/*!
 * Initialize this analytic. This implementation checks to make sure the input
 * data object and output file have been set.
 */
void ExportExpressionMatrix::initialize()
{
   EDEBUG_FUNC(this);

   if ( !_input || !_output )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Did not get valid input and/or output arguments."));
      throw e;
   }
}
