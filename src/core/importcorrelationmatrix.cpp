#include "importcorrelationmatrix.h"
#include "importcorrelationmatrix_input.h"
#include "datafactory.h"






/*!
 * Return the total number of blocks this analytic must process as steps
 * or blocks of work. This implementation uses a work block for each line
 * of the input file.
 */
int ImportCorrelationMatrix::size() const
{
   EDEBUG_FUNC(this);

   return _numLines;
}







/*!
 * Process the given index with a possible block of results if this analytic
 * produces work blocks. This implementation uses only the index of the result
 * block to determine which piece of work to do.
 *
 * @param result
 */
void ImportCorrelationMatrix::process(const EAbstractAnalyticBlock* result)
{
   EDEBUG_FUNC(this,result);

   // read a line from input file
   QString line = _stream.readLine();
   auto words = line.splitRef(QRegExp("\\s+"), QString::SkipEmptyParts);

   // make sure the line is valid
   if ( words.size() == 7 )
   {
      int geneX = words[0].toInt();
      int geneY = words[1].toInt();
      float correlation = words[5].toFloat();
      QStringRef sampleMask = words[6];

      // make sure sample mask has correct length
      if ( sampleMask.size() != _sampleSize )
      {
         E_MAKE_EXCEPTION(e);
         e.setTitle(tr("Parsing Error"));
         e.setDetails(tr("Encountered sample mask with invalid length %1. Sample size is %2.")
            .arg(sampleMask.size())
            .arg(_sampleSize));
         throw e;
      }

      // save previous pair when new pair is read
      Pairwise::Index nextIndex(geneX, geneY);

      if ( _index != nextIndex )
      {
         // save pairs
         if ( _ccmPair.clusterSize() > 1 )
         {
            _ccmPair.write(_index);
         }

         if ( _cmxPair.clusterSize() > 0 )
         {
            _cmxPair.write(_index);
         }

         // reset pairs
         _ccmPair.clearClusters();
         _cmxPair.clearClusters();

         // update index
         _index = nextIndex;
      }

      // append data to ccm pair and cmx pair
      int cluster = _ccmPair.clusterSize();

      _ccmPair.addCluster();
      _cmxPair.addCluster();

      for ( int i = 0; i < sampleMask.size(); ++i )
      {
         _ccmPair.at(cluster, i) = static_cast<qint8>(sampleMask[i].digitValue());
      }

      _cmxPair.at(cluster) = correlation;
   }

   // otherwise throw an error
   else
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Parsing Error"));
      e.setDetails(tr("Encountered line with incorrect amount of fields. "
                      "Read %1 fields when there should have been %2.")
         .arg(words.size())
         .arg(7));
      throw e;
   }

   // save last pair
   if ( result->index() == _numLines - 1 )
   {
      if ( _ccmPair.clusterSize() > 1 )
      {
         _ccmPair.write(_index);
      }

      if ( _cmxPair.clusterSize() > 0 )
      {
         _cmxPair.write(_index);
      }
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






/*!
 * Make a new input object and return its pointer.
 */
EAbstractAnalyticInput* ImportCorrelationMatrix::makeInput()
{
   EDEBUG_FUNC(this);

   return new Input(this);
}






/*!
 * Initialize this analytic. This implementation checks to make sure the input
 * file and output data objects have been set, and that a correlation name was
 * provided.
 */
void ImportCorrelationMatrix::initialize()
{
   EDEBUG_FUNC(this);

   // make sure input/output arguments are valid
   if ( !_input || !_ccm || !_cmx )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Did not get valid input and/or output arguments."));
      throw e;
   }

   // make sure correlation name is valid
   if ( _correlationName.isEmpty() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Correlation name is required."));
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

   // return stream to beginning of the input file
   _stream.seek(0);

   // make sure reading input file worked
   if ( _stream.status() != QTextStream::Ok )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Qt Text Stream encountered an unknown error."));
      throw e;
   }

   // build gene name metadata
   EMetaArray metaGeneNames;
   for ( int i = 0; i < _geneSize; ++i )
   {
      metaGeneNames.append(QString::number(i));
   }

   // build sample name metadata
   EMetaArray metaSampleNames;
   for ( int i = 0; i < _sampleSize; ++i )
   {
      metaSampleNames.append(QString::number(i));
   }

   // initialize output data
   _ccm->initialize(metaGeneNames, _maxClusterSize, metaSampleNames);
   _cmx->initialize(metaGeneNames, _maxClusterSize, _correlationName);

   // initialize pairwise iterators
   _ccmPair = CCMatrix::Pair(_ccm);
   _cmxPair = CorrelationMatrix::Pair(_cmx);
}
