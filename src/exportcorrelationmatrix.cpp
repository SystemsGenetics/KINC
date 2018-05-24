#include "exportcorrelationmatrix.h"
#include "exportcorrelationmatrix_input.h"
#include "datafactory.h"



int ExportCorrelationMatrix::size() const
{
   return 1;
}






void ExportCorrelationMatrix::process(const EAbstractAnalytic::Block* /*result*/)
{
   // initialize pair iterators
   CorrelationMatrix::Pair cmxPair(_cmx);
   CCMatrix::Pair ccmPair(_ccm);

   // initialize workspace
   QString sampleMask(_ccm->sampleSize(), '0');

   // create text stream to output file and write until end reached
   QTextStream stream(_output);
   stream.setRealNumberPrecision(6);

   // iterate through all pairs
   while ( cmxPair.hasNext() )
   {
      // read next pair
      cmxPair.readNext();

      if ( cmxPair.clusterSize() > 1 )
      {
         ccmPair.read(cmxPair.index());
      }

      // write pairwise data to output file
      for ( int k = 0; k < cmxPair.clusterSize(); k++ )
      {
         float correlation = cmxPair.at(k, 0);
         int numSamples = 0;
         int numMissing = 0;
         int numPostOutliers = 0;
         int numPreOutliers = 0;
         int numThreshold = 0;

         // if there are multiple clusters then use cluster data
         if ( cmxPair.clusterSize() > 1 )
         {
            // compute summary statistics
            for ( int i = 0; i < _ccm->sampleSize(); i++ )
            {
               switch ( ccmPair.at(k, i) )
               {
               case 1:
                  numSamples++;
                  break;
               case 6:
                  numThreshold++;
                  break;
               case 7:
                  numPreOutliers++;
                  break;
               case 8:
                  numPostOutliers++;
                  break;
               case 9:
                  numMissing++;
                  break;
               }
            }

            // write sample mask to string
            for ( int i = 0; i < _ccm->sampleSize(); i++ )
            {
               sampleMask[i] = '0' + ccmPair.at(k, i);
            }
         }

         // else just initialize empty sample mask
         else
         {
            sampleMask.fill('0');
         }

         // write cluster to output file
         stream
            << cmxPair.index().getX()
            << "\t" << cmxPair.index().getY()
            << "\t" << k
            << "\t" << cmxPair.clusterSize()
            << "\t" << numSamples
            << "\t" << numMissing
            << "\t" << numPostOutliers
            << "\t" << numPreOutliers
            << "\t" << numThreshold
            << "\t" << correlation
            << "\t" << sampleMask
            << "\n";
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






EAbstractAnalytic::Input* ExportCorrelationMatrix::makeInput()
{
   return new Input(this);
}






void ExportCorrelationMatrix::initialize()
{
   if ( !_ccm || !_cmx || !_output )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Did not get valid input and/or output arguments."));
      throw e;
   }
}
