#include "extract.h"
#include "extract_input.h"
#include "datafactory.h"
#include "expressionmatrix_gene.h"



using namespace std;






int Extract::size() const
{
   return 1;
}






void Extract::process(const EAbstractAnalytic::Block* result)
{
   Q_UNUSED(result);

   // initialize pair iterators
   CorrelationMatrix::Pair cmxPair(_cmx);
   CCMatrix::Pair ccmPair(_ccm);

   // get gene names
   EMetaArray geneNames {_cmx->geneNames().toArray()};

   // initialize workspace
   QString sampleMask(_ccm->sampleSize(), '0');

   // create text stream to output file and write until end reached
   QTextStream stream(_output);
   stream.setRealNumberPrecision(6);

   // write header to file
   stream
      << "Source"
      << "\t" << "Target"
      << "\t" << "sc"
      << "\t" << "Interaction"
      << "\t" << "Cluster"
      << "\t" << "Num_Clusters"
      << "\t" << "Cluster_Samples"
      << "\t" << "Missing_Samples"
      << "\t" << "Cluster_Outliers"
      << "\t" << "Pair_Outliers"
      << "\t" << "Too_Low"
      << "\t" << "Samples"
      << "\n";

   // increment through all gene pairs
   while ( cmxPair.hasNext() )
   {
      // read next gene pair
      cmxPair.readNext();

      if ( cmxPair.clusterSize() > 1 )
      {
         ccmPair.read(cmxPair.index());
      }

      // write gene pair data to output file
      for ( int k = 0; k < cmxPair.clusterSize(); k++ )
      {
         auto& source {geneNames.at(cmxPair.index().getX()).toString()};
         auto& target {geneNames.at(cmxPair.index().getY()).toString()};
         float correlation {cmxPair.at(k, 0)};
         QString interaction {"co"};
         int numSamples {0};
         int numMissing {0};
         int numPostOutliers {0};
         int numPreOutliers {0};
         int numThreshold {0};

         // exclude cluster if correlation is not within thresholds
         if ( fabs(correlation) < _minCorrelation || _maxCorrelation < fabs(correlation) )
         {
            continue;
         }

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

         // otherwise use expression data
         else
         {
            // read in gene expressions
            ExpressionMatrix::Gene gene1(_emx);
            ExpressionMatrix::Gene gene2(_emx);

            gene1.read(cmxPair.index().getX());
            gene2.read(cmxPair.index().getY());

            // determine sample mask from expression data
            for ( int i = 0; i < _emx->sampleSize(); ++i )
            {
               if ( isnan(gene1.at(i)) || isnan(gene2.at(i)) )
               {
                  sampleMask[i] = '9';
               }
               else
               {
                  sampleMask[i] = '1';
               }
            }
         }

         // write cluster to output file
         stream
            << source
            << "\t" << target
            << "\t" << correlation
            << "\t" << interaction
            << "\t" << k
            << "\t" << cmxPair.clusterSize()
            << "\t" << numSamples
            << "\t" << numMissing
            << "\t" << numPostOutliers
            << "\t" << numPreOutliers
            << "\t" << numThreshold
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

   // reset gene pair iterator
   cmxPair.reset();

   // create text stream to graphml file and write until end reached
   stream.setDevice(_graphml);

   // write header to file
   stream
      << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
      << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"\n"
      << "    xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
      << "    xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n"
      << "  <graph id=\"G\" edgedefault=\"undirected\">\n";

   // write each node to file
   for ( int i = 0; i < _cmx->geneSize(); i++ )
   {
      auto& id {geneNames.at(i).toString()};

      stream << "    <node id=\"" << id << "\"/>\n";
   }

   // increment through all gene pairs
   while ( cmxPair.hasNext() )
   {
      // read next gene pair
      cmxPair.readNext();

      if ( cmxPair.clusterSize() > 1 )
      {
         ccmPair.read(cmxPair.index());
      }

      // write gene pair edges to file
      for ( int k = 0; k < cmxPair.clusterSize(); k++ )
      {
         auto& source {geneNames.at(cmxPair.index().getX()).toString()};
         auto& target {geneNames.at(cmxPair.index().getY()).toString()};
         float correlation {cmxPair.at(k, 0)};

         // exclude edge if correlation is not within thresholds
         if ( fabs(correlation) < _minCorrelation || _maxCorrelation < fabs(correlation) )
         {
            continue;
         }

         // if there are multiple clusters then use cluster data
         if ( cmxPair.clusterSize() > 1 )
         {
            // write sample mask to string
            for ( int i = 0; i < _ccm->sampleSize(); i++ )
            {
               sampleMask[i] = '0' + ccmPair.at(k, i);
            }
         }

         // otherwise use expression data
         else
         {
            // read in gene expressions
            ExpressionMatrix::Gene gene1(_emx);
            ExpressionMatrix::Gene gene2(_emx);

            gene1.read(cmxPair.index().getX());
            gene2.read(cmxPair.index().getY());

            // determine sample mask from expression data
            for ( int i = 0; i < _emx->sampleSize(); ++i )
            {
               if ( isnan(gene1.at(i)) || isnan(gene2.at(i)) )
               {
                  sampleMask[i] = '9';
               }
               else
               {
                  sampleMask[i] = '1';
               }
            }
         }

         // write edge to file
         stream
            << "    <edge"
            << " source=\"" << source << "\""
            << " target=\"" << target << "\""
            << " samples=\"" << sampleMask << "\""
            << "/>\n";
      }
   }

   // write footer to file
   stream
      << "  </graph>\n"
      << "</graphml>\n";

   // make sure writing graphml file worked
   if ( stream.status() != QTextStream::Ok )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Qt Text Stream encountered an unknown error."));
      throw e;
   }
}






EAbstractAnalytic::Input* Extract::makeInput()
{
   return new Input(this);
}






void Extract::initialize()
{
   if ( !_emx || !_ccm || !_cmx || !_output )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Did not get valid input and/or output arguments."));
      throw e;
   }
}
