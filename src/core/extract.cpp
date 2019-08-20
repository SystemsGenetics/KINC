#include "extract.h"
#include "extract_input.h"
#include "datafactory.h"
#include "expressionmatrix_gene.h"
#include "condition-specificclustersmatrix_pair.h"


using namespace std;






/*!
 * Return the total number of blocks this analytic must process as steps
 * or blocks of work. This implementation uses a work block for writing
 * each pair to the output file.
 */
int Extract::size() const
{
   EDEBUG_FUNC(this);

   return _cmx->size();
}






/*!
 * Process the given index with a possible block of results if this analytic
 * produces work blocks. This implementation uses only the index of the result
 * block to determine which piece of work to do.
 *
 * @param result
 */
void Extract::process(const EAbstractAnalyticBlock* result)
{
   EDEBUG_FUNC(this,result);

   // write pair according to the output format
   switch ( _outputFormat )
   {
   case OutputFormat::Text:
      writeTextFormat(result->index());
      break;
   case OutputFormat::Minimal:
      writeMinimalFormat(result->index());
      break;
   case OutputFormat::GraphML:
      writeGraphMLFormat(result->index());
      break;
   }
}






/*!
 * Write the next pair using the text format.
 *
 * @param index
 */
void Extract::writeTextFormat(int index)
{
   EDEBUG_FUNC(this);

   // get gene names
   EMetaArray geneNames {_cmx->geneNames()};

   // initialize workspace
   QString sampleMask(_ccm->sampleSize(), '0');



   // write header to file
   if ( index == 0 )
   {
      _stream
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
         << "\t" << "Samples";
      if(_amx)
      {
          QVector<QString> headerTitles(formatAnnotations());
          for(int i = 0; i < headerTitles.size(); i++)
          {
              _stream << "\t" << headerTitles.at(i);
          }
      }
      _stream << "\n";
   }

   // read next pair
   _cmxPair.readNext();
   _ccmPair.read(_cmxPair.index());

   // write pairwise data to output file
   for ( int k = 0; k < _cmxPair.clusterSize(); k++ )
   {
      QString source {geneNames.at(_cmxPair.index().getX()).toString()};
      QString target {geneNames.at(_cmxPair.index().getY()).toString()};
      float correlation {_cmxPair.at(k)};
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

      // if cluster data exists then use it
      if ( _ccmPair.clusterSize() > 0 )
      {
         // compute summary statistics
         for ( int i = 0; i < _ccm->sampleSize(); i++ )
         {
            switch ( _ccmPair.at(k, i) )
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
            sampleMask[i] = '0' + _ccmPair.at(k, i);
         }
      }

      // otherwise use expression data
      else
      {
         // read in gene expressions
         ExpressionMatrix::Gene gene1(_emx);
         ExpressionMatrix::Gene gene2(_emx);

         gene1.read(_cmxPair.index().getX());
         gene2.read(_cmxPair.index().getY());

         // determine sample mask, summary statistics from expression data
         for ( int i = 0; i < _emx->sampleSize(); ++i )
         {
            if ( isnan(gene1.at(i)) || isnan(gene2.at(i)) )
            {
               sampleMask[i] = '9';
               numMissing++;
            }
            else
            {
               sampleMask[i] = '1';
               numSamples++;
            }
         }
      }



      // write cluster to output file
      _stream
         << source
         << "\t" << target
         << "\t" << correlation
         << "\t" << interaction
         << "\t" << k
         << "\t" << _cmxPair.clusterSize()
         << "\t" << numSamples
         << "\t" << numMissing
         << "\t" << numPostOutliers
         << "\t" << numPreOutliers
         << "\t" << numThreshold
         << "\t" << sampleMask;

      if(_amx)
      {
          //initialize cscm for data values to tack on
          CSCM::Pair cscm(_cscm);
          cscm.read(_cmxPair.index());
          for(int i = 0; i < _cscm->getTestCount(); i++)
          {
              _stream << "\t" << cscm.at(k, i);
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
 * Write the next pair using the minimal format.
 *
 * @param index
 */
void Extract::writeMinimalFormat(int index)
{
   EDEBUG_FUNC(this);

   // get gene names
   EMetaArray geneNames {_cmx->geneNames()};

   // write header to file
   if ( index == 0 )
   {
      _stream
         << "Source"
         << "\t" << "Target"
         << "\t" << "sc"
         << "\t" << "Cluster"
         << "\t" << "Num_Clusters"
         << "\n";
   }

   // read next pair
   _cmxPair.readNext();

   // write pairwise data to output file
   for ( int k = 0; k < _cmxPair.clusterSize(); k++ )
   {
      QString source {geneNames.at(_cmxPair.index().getX()).toString()};
      QString target {geneNames.at(_cmxPair.index().getY()).toString()};
      float correlation {_cmxPair.at(k)};

      // exclude cluster if correlation is not within thresholds
      if ( fabs(correlation) < _minCorrelation || _maxCorrelation < fabs(correlation) )
      {
         continue;
      }

      // write cluster to output file
      _stream
         << source
         << "\t" << target
         << "\t" << correlation
         << "\t" << k
         << "\t" << _cmxPair.clusterSize()
         << "\n";
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
 * Write the next pair using the GraphML format.
 *
 * @param index
 */
void Extract::writeGraphMLFormat(int index)
{
   EDEBUG_FUNC(this);

   // get gene names
   EMetaArray geneNames {_cmx->geneNames()};

   // initialize workspace
   QString sampleMask(_ccm->sampleSize(), '0');

   if ( index == 0 )
   {
      // write header to file
      _stream
         << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
         << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"\n"
         << "    xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
         << "    xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n"
         << "  <graph id=\"G\" edgedefault=\"undirected\">\n";

      // write node list to file
      for ( int i = 0; i < _cmx->geneSize(); i++ )
      {
         QString id {geneNames.at(i).toString()};

         _stream << "    <node id=\"" << id << "\"/>\n";
      }
   }

   // read next pair
   _cmxPair.readNext();

   if ( _cmxPair.clusterSize() > 1 )
   {
      _ccmPair.read(_cmxPair.index());
   }

   // write pairwise data to net file
   for ( int k = 0; k < _cmxPair.clusterSize(); k++ )
   {
      QString source {geneNames.at(_cmxPair.index().getX()).toString()};
      QString target {geneNames.at(_cmxPair.index().getY()).toString()};
      float correlation {_cmxPair.at(k)};

      // exclude edge if correlation is not within thresholds
      if ( fabs(correlation) < _minCorrelation || _maxCorrelation < fabs(correlation) )
      {
         continue;
      }

      // if there are multiple clusters then use cluster data
      if ( _cmxPair.clusterSize() > 1 )
      {
         // write sample mask to string
         for ( int i = 0; i < _ccm->sampleSize(); i++ )
         {
            sampleMask[i] = '0' + _ccmPair.at(k, i);
         }
      }

      // otherwise use expression data
      else
      {
         // read in gene expressions
         ExpressionMatrix::Gene gene1(_emx);
         ExpressionMatrix::Gene gene2(_emx);

         gene1.read(_cmxPair.index().getX());
         gene2.read(_cmxPair.index().getY());

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
      _stream
         << "    <edge"
         << " source=\"" << source << "\""
         << " target=\"" << target << "\""
         << " samples=\"" << sampleMask << "\""
         << "/>\n";
   }

   // write footer to file
   if ( index == size() - 1 )
   {
      _stream
         << "  </graph>\n"
         << "</graphml>\n";
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
EAbstractAnalyticInput* Extract::makeInput()
{
   EDEBUG_FUNC(this);

   return new Input(this);
}






/*!
 * Initialize this analytic. This implementation checks to make sure the input
 * data objects and output file have been set.
 */
void Extract::initialize()
{
   EDEBUG_FUNC(this);

   // make sure input/output arguments are valid
   if ( !_emx || !_cmx || !_output )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Did not get valid input and/or output arguments."));
      throw e;
   }

   if ( _outputFormat != OutputFormat::Minimal && !_ccm )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("--ccm is required for all output formats except minimal."));
      throw e;
   }

   // initialize pairwise iterators
   _ccmPair = CCMatrix::Pair(_ccm);
   _cmxPair = CorrelationMatrix::Pair(_cmx);

   // initialize output file stream
   _stream.setDevice(_output);
   _stream.setRealNumberPrecision(8);
}




/*!
 * Implements an interface to crate the header names for the tests done in the
 * to create the cscm.
 *
 * @return Vector of header names to append to the network file.
 */
QVector<QString> Extract::formatAnnotations()
{
    //set up variables
    QVector<QString> headerInfo;
    QTextStream file(_amx);
    QVector<QVector<QString>> info;
    QString delimiter = ",";


    auto line = file.readLine();
    //decide on delimiter
    if(line.contains("\t"))
        delimiter = "\t";
    //grab the features
    auto words = line.split(delimiter);
    for(int i = 0; i < words.size(); i++)
    {
        info.append(QVector<QString>());
        info[i].append(words.at(i));
    }
    //grab the rest of the label information
    while(!file.atEnd())
    {
        auto line = file.readLine();
        auto words = line.split("\t");
        for(int i = 0; i < words.size(); i++)
        {
            if(!info.at(i).contains(words.at(i)))
            {
                info[i].append(words.at(i));
            }
        }
    }
    EMetaObject features = _cscm->getFeatures();
    //for each feature
    for(int i = 0; i < info.size(); i++)
    {
        //for each label in the feature
        for(int j = 1; j < info.at(i).size(); j++)
        {
            //for each cluster
            if(features.at(info.at(i).at(0)).toObject().at("Test").toObject().at("Type").toString() == "Catagorical")
            {
                headerInfo.append(info.at(i).at(0) + "__" + info.at(i).at(j));
            }
        }
    }
    return headerInfo;
}
