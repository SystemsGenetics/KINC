#include "importcorrelationmatrix.h"
#include "importcorrelationmatrix_input.h"
#include "ccmatrix_pair.h"
#include "correlationmatrix_pair.h"
#include "datafactory.h"
#include "pairwise_index.h"






int ImportCorrelationMatrix::size() const
{
   return 1;
}







void ImportCorrelationMatrix::process(const EAbstractAnalytic::Block* result)
{
   Q_UNUSED(result);

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

   // build correlation name metadata
   EMetaArray metaCorrelationNames;
   metaCorrelationNames.append(_correlationName);

   // initialize output data
   _ccm->initialize(metaGeneNames, _maxClusterSize, metaSampleNames);
   _cmx->initialize(metaGeneNames, _maxClusterSize, metaCorrelationNames);

   Pairwise::Index index;
   CCMatrix::Pair ccmPair(_ccm);
   CorrelationMatrix::Pair cmxPair(_cmx);

   // create text stream from input file and read until end reached
   QTextStream stream(_input);
   while ( !stream.atEnd() )
   {
      // read a line from text file
      QString line = stream.readLine();
      auto words = line.splitRef(QRegExp("\\s+"), QString::SkipEmptyParts);

      // make sure the line is valid
      if ( words.size() == 11 )
      {
         int geneX = words[0].toInt();
         int geneY = words[1].toInt();
         float correlation = words[9].toFloat();
         QStringRef sampleMask = words[10];

         // make sure sample mask has correct length
         if ( sampleMask.size() != _sampleSize )
         {
            E_MAKE_EXCEPTION(e);
            e.setTitle(tr("Parsing Error"));
            e.setDetails(tr("Encountered sample mask with invalid length %1. "
                            "Sample size is %2.")
                         .arg(sampleMask.size()).arg(_sampleSize));
            throw e;
         }

         // save previous pair when new pair is read
         Pairwise::Index nextIndex(geneX, geneY);

         if ( index != nextIndex )
         {
            // save pairs
            if ( ccmPair.clusterSize() > 1 )
            {
               ccmPair.write(index);
            }

            if ( cmxPair.clusterSize() > 0 )
            {
               cmxPair.write(index);
            }

            // reset pairs
            ccmPair.clearClusters();
            cmxPair.clearClusters();

            // update index
            index = nextIndex;
         }

         // append data to ccm pair and cmx pair
         int cluster = ccmPair.clusterSize();

         ccmPair.addCluster();
         cmxPair.addCluster();

         for ( int i = 0; i < sampleMask.size(); ++i )
         {
            ccmPair.at(cluster, i) = sampleMask[i].digitValue();
         }

         cmxPair.at(cluster, 0) = correlation;
      }

      // save last pair
      if ( ccmPair.clusterSize() > 1 )
      {
         ccmPair.write(index);
      }

      if ( cmxPair.clusterSize() > 0 )
      {
         cmxPair.write(index);
      }

      // skip empty lines and lines with '#' markers
      else if ( words.size() != 1 && words.size() != 0 )
      {
         continue;
      }

      // otherwise throw an error
      else
      {
         E_MAKE_EXCEPTION(e);
         e.setTitle(tr("Parsing Error"));
         e.setDetails(tr("Encountered line with incorrect amount of fields. "
                         "Read %1 fields when there should have been %2.")
                      .arg(words.size()).arg(11));
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
}






EAbstractAnalytic::Input* ImportCorrelationMatrix::makeInput()
{
   return new Input(this);
}






void ImportCorrelationMatrix::initialize()
{
   if ( !_input || !_ccm || !_cmx )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Did not get valid input and/or output arguments."));
      throw e;
   }

   if ( _correlationName.isEmpty() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Correlation name is required."));
      throw e;
   }
}
