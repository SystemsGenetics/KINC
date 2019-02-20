#include "importexpressionmatrix.h"
#include "importexpressionmatrix_input.h"
#include "datafactory.h"
#include "expressionmatrix_gene.h"






/*!
 * Return the total number of blocks this analytic must process as steps
 * or blocks of work. This implementation uses a work block for reading
 * each line of the input file, plus one work block to create the output
 * data object.
 */
int ImportExpressionMatrix::size() const
{
   EDEBUG_FUNC(this);

   return _numLines + 1;
}






/*!
 * Process the given index with a possible block of results if this analytic
 * produces work blocks. This implementation uses only the index of the result
 * block to determine which piece of work to do.
 *
 * @param result
 */
void ImportExpressionMatrix::process(const EAbstractAnalytic::Block* result)
{
   EDEBUG_FUNC(this, result);

   // read or create the sample names in the first step
   if ( result->index() == 0 )
   {
      // seek to the beginning of the input file
      _stream.seek(0);

      // if sample size is not zero then build sample name list
      if ( _sampleSize != 0 )
      {
         for (int i = 0; i < _sampleSize ;++i)
         {
            _sampleNames.append(QString::number(i));
         }
      }

      // otherwise read sample names from first line
      else
      {
         // read a line from the input file
         QString line = _stream.readLine();
         auto words = line.splitRef(QRegExp("\\s+"), QString::SkipEmptyParts);

         // parse the sample names
         _sampleSize = words.size();

         for ( auto& word : words )
         {
            _sampleNames.append(word.toString());
         }

         // make sure reading input file worked
         if ( _stream.status() != QTextStream::Ok )
         {
            E_MAKE_EXCEPTION(e);
            e.setTitle(tr("File IO Error"));
            e.setDetails(tr("Qt Text Stream encountered an unknown error."));
            throw e;
         }
      }
   }

   // read each gene from the input file in a separate step
   else if ( result->index() < _numLines )
   {
      // read a line from the input file
      QString line = _stream.readLine();
      auto words = line.splitRef(QRegExp("\\s+"), QString::SkipEmptyParts);

      // make sure the number of words matches expected sample size
      if ( words.size() == _sampleSize + 1 )
      {
         // read row from text file into gene
         Gene gene(_sampleSize);

         for ( int i = 1; i < words.size(); ++i )
         {
            // if word matches the nan token then set it as such
            if ( words.at(i) == _nanToken )
            {
               gene.expressions[i-1] = NAN;
            }

            // else this is a normal floating point expression
            else
            {
               // read in the floating point value
               bool ok;
               gene.expressions[i-1] = words.at(i).toFloat(&ok);

               // make sure reading worked
               if ( !ok )
               {
                  E_MAKE_EXCEPTION(e);
                  e.setTitle(tr("Parsing Error"));
                  e.setDetails(tr("Failed to read expression value \"%1\" for gene %2.")
                     .arg(words.at(i).toString())
                     .arg(words.at(0).toString()));
                  throw e;
               }
            }
         }

         // append gene data and gene name
         _genes.append(gene);
         _geneNames.append(words.at(0).toString());
      }

      // otherwise throw an error
      else
      {
         E_MAKE_EXCEPTION(e);
         e.setTitle(tr("Parsing Error"));
         e.setDetails(tr("Encountered gene expression line with incorrect amount of fields. "
                         "Read in %1 fields when it should have been %2. Gene name is %3.")
            .arg(words.size()-1)
            .arg(_sampleSize)
            .arg(words.at(0).toString()));
         throw e;
      }

      // make sure reading input file worked
      if ( _stream.status() != QTextStream::Ok )
      {
         E_MAKE_EXCEPTION(e);
         e.setTitle(tr("File IO Error"));
         e.setDetails(tr("Qt Text Stream encountered an unknown error."));
         throw e;
      }
   }

   // create the output data object in the final step
   else if ( result->index() == _numLines )
   {
      // initialize expression matrix
      _output->initialize(_geneNames, _sampleNames);

      // iterate through each gene
      ExpressionMatrix::Gene gene(_output);

      for ( int i = 0; i < _output->geneSize(); ++i )
      {
         // save each gene to expression matrix
         for ( int j = 0; j < _output->sampleSize(); ++j )
         {
            gene[j] = _genes[i].expressions[j];
         }

         gene.write(i);
      }
   }
}






/*!
 * Make a new input object and return its pointer.
 */
EAbstractAnalytic::Input* ImportExpressionMatrix::makeInput()
{
   EDEBUG_FUNC(this);

   return new Input(this);
}






/*!
 * Initialize this analytic. This implementation checks to make sure the input
 * file and output data object have been set.
 */
void ImportExpressionMatrix::initialize()
{
   EDEBUG_FUNC(this);

   // make sure input/output arguments are valid
   if ( !_input || !_output )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Did not get valid input and/or output arguments."));
      throw e;
   }

   // initialize input file stream
   _stream.setDevice(_input);

   // count the number of lines in the input file
   _numLines = 0;

   while ( !_stream.atEnd() )
   {
      _stream.readLine();
      _numLines++;
   }

   // make sure reading input file worked
   if ( _stream.status() != QTextStream::Ok )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Qt Text Stream encountered an unknown error."));
      throw e;
   }
}
