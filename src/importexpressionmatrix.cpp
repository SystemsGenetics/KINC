#include "importexpressionmatrix.h"
#include "importexpressionmatrix_input.h"
#include "datafactory.h"






int ImportExpressionMatrix::size() const
{
   return 1;
}






void ImportExpressionMatrix::process(const EAbstractAnalytic::Block* result)
{
   Q_UNUSED(result);

   // use expression declaration
   using Expression = ExpressionMatrix::Expression;

   // structure for building list of genes
   struct Gene
   {
      Gene(int size)
      {
         expressions = new Expression[size];
      }
      ~Gene()
      {
         delete[] expressions;
      }

      Expression* expressions;
   };

   // initialize gene expression linked list
   QList<Gene *> genes;

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
      // read a line from text file
      QString line = stream.readLine();
      auto words = line.splitRef(QRegExp("\\s+"), QString::SkipEmptyParts);

      // read sample names from first line
      if ( _sampleSize == 0 )
      {
         _sampleSize = words.size();
         for ( auto& word : words )
         {
            sampleNames.append(word.toString());
         }
      }

      // make sure the number of words matches expected sample size
      else if ( words.size() == _sampleSize + 1 )
      {
         // read row from text file into gene
         Gene* gene {new Gene(_sampleSize)};

         for ( int i = 1; i < words.size(); ++i )
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
                  e.setDetails(tr("Failed to read expression value \"%1\" for gene %2.")
                               .arg(words.at(i).toString()).arg(words.at(0).toString()));
                  throw e;
               }

               // apply transform and append expression to gene
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

         // append gene data and gene name
         genes.append(gene);
         geneNames.append(words.at(0).toString());
      }

      // otherwise throw an error
      else
      {
         E_MAKE_EXCEPTION(e);
         e.setTitle(tr("Parsing Error"));
         e.setDetails(tr("Encountered gene expression line with incorrect amount of fields. "
                         "Read in %1 fields when it should have been %2. Gene name is %3.")
                      .arg(words.size()-1).arg(_sampleSize).arg(words.at(0).toString()));
         throw e;
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

   // initialize expression matrix
   _output->initialize(geneNames, sampleNames);

   // iterate through each gene
   ExpressionMatrix::Gene gene(_output);
   for ( int i = 0; i < _output->getGeneSize(); ++i )
   {
      // save each gene to expression matrix
      for ( int j = 0; j < _output->getSampleSize(); ++j )
      {
         gene[j] = genes[i]->expressions[j];
      }

      gene.write(i);
   }

   // set transform used in expression matrix
   _output->setTransform(_transform);
}






EAbstractAnalytic::Input* ImportExpressionMatrix::makeInput()
{
   return new Input(this);
}






void ImportExpressionMatrix::initialize()
{
   if ( !_input || !_output )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Did not get valid input and/or output arguments."));
      throw e;
   }
}
