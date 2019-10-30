#include "extract.h"
#include "extract_input.h"
#include "datafactory.h"
#include "expressionmatrix_gene.h"
#include "conditionspecificclustersmatrix_pair.h"



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
         << "\t" << "Similarity_Score"
         << "\t" << "Interaction"
         << "\t" << "Cluster_Index"
         << "\t" << "Cluster_Size"
         << "\t" << "Samples";

      // If the condition-specific matrix is provided then add some
      // additional headers for each test.
      if ( _csm )
      {
          for ( int i = 0; i < _csm->getTestCount(); i++ )
          {
              QString testName = _csm->getTestName(i);
              QString testType = _csm->getTestType(i);
              _stream << "\t" << testName + + "_pVal";
              // If there is an R-squared value add it as well.
              if (QString::compare(testType, "Quantitative", Qt::CaseInsensitive) == 0 ||
                  QString::compare(testType, "Ordinal", Qt::CaseInsensitive) == 0) {
                  _stream << "\t" << testName + "_RSqr";
              }
          }
      }

      _stream << "\n";
   }

   // read next pair
   _cmxPair.readNext();
   _ccmPair.read(_cmxPair.index());
   if ( _csm )
   {
      _csmPair.read(_cmxPair.index());
   }

   // write pairwise data to output file
   for ( int k = 0; k < _cmxPair.clusterSize(); k++ )
   {
      QString source {geneNames.at(_cmxPair.index().getX()).toString()};
      QString target {geneNames.at(_cmxPair.index().getY()).toString()};
      float correlation {_cmxPair.at(k)};
      QString interaction {"co"};
      int numSamples {0};

      // exclude cluster if correlation is not within thresholds
      if ( fabs(correlation) < _minCorrelation || _maxCorrelation < fabs(correlation) )
      {
         continue;
      }

      // If the condition-specific matrix is provided and there is a p-value filter then
      // we want to exclude values that are higher than the p-value(s).
      if ( _csm )
      {
          // Run some pre-checks on any filters provided.
          pValueFilterCheck();
          int notInclude = 0;
          int include = 0;

          // Iterate through teach test.
          for ( int i = 0; i < _csm->getTestCount(); i++ )
          {
              // Count the number of features that are not signficant for each feature-specific p-value provided.
              if ( _csmPValueFilterFeatureNames.size() != 0 && !PValuefilter(_csm->getTestName(i), _csmPair.at(k, i, "pvalue")))
              {
                  notInclude++;
              }
              // Count the number of features that are signficant if only a single p-value filter was given.
              if ( _csmPValueFilterFeatureNames.size() == 0 && PValuefilter(_csm->getTestName(i), _csmPair.at(k, i, "pvalue")))
              {
                  include++;
              }
          }
          if ( notInclude > 0  || (_csmPValueFilterFeatureNames.size() == 0 && include == 0))
          {
              continue;
          }

          // Run some pre-checks on any filters provided.
          rSquareFilterCheck();
          notInclude = 0;
          include = 0;

          // Iterate through teach test.
          for ( int i = 0; i < _csm->getTestCount(); i++ )
          {
              // Count the number of features that are not signficant for each feature-specific p-value provided.
              if ( _csmRSquareFilterFeatureNames.size() != 0 && !RSquarefilter(_csm->getTestName(i), _csmPair.at(k, i, "r2")))
              {
                  notInclude++;
              }
              // Count the number of features that are signficant if only a single p-value filter was given.
              if ( _csmRSquareFilterFeatureNames.size() == 0 && RSquarefilter(_csm->getTestName(i), _csmPair.at(k, i, "r2")))
              {
                  include++;
              }
          }
          if ( notInclude > 0  || (_csmPValueFilterFeatureNames.size() == 0 && include == 0))
          {
              continue;
          }
      }

      // if cluster data exists then use it
      if ( _ccmPair.clusterSize() > 0 )
      {
         // compute cluster size
         for ( int i = 0; i < _ccm->sampleSize(); i++ )
         {
            if ( _ccmPair.at(k, i) == 1 )
            {
               numSamples++;
            }
         }

         // write sample mask to string
         for ( int i = 0; i < _ccm->sampleSize(); i++ )
         {
            sampleMask[i] = '0' + _ccmPair.at(k, i);
         }
      }

      // otherwise use expression data if provided
      else if ( _emx )
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
            }
            else
            {
               sampleMask[i] = '1';
               numSamples++;
            }
         }
      }

      // otherwise throw an error
      else
      {
         E_MAKE_EXCEPTION(e);
         e.setTitle(tr("Invalid Input"));
         e.setDetails(tr("Expression Matrix was not provided but Cluster Matrix is missing sample data."));
         throw e;
      }

      // write cluster to output file
      _stream
         << source
         << "\t" << target
         << "\t" << correlation
         << "\t" << interaction
         << "\t" << k + 1
         << "\t" << numSamples
         << "\t" << sampleMask;

      if ( _csm )
      {
          for ( int i = 0; i < _csm->getTestCount(); i++ )
          {
              _stream << "\t" << _csmPair.at(k, i, "pvalue");
              QString testType = _csm->getTestType(i);
              // If there is an R-squared value add it as well.
              if (QString::compare(testType, "Quantitative", Qt::CaseInsensitive) == 0 ||
                  QString::compare(testType, "Ordinal", Qt::CaseInsensitive) == 0) {
                  _stream << "\t" << _csmPair.at(k, i, "r2");
              }

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
         << "\t" << k + 1
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

      // otherwise use expression data if provided
      else if ( _emx )
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

      // otherwise throw an error
      else
      {
         E_MAKE_EXCEPTION(e);
         e.setTitle(tr("Invalid Input"));
         e.setDetails(tr("Expression Matrix was not provided but Cluster Matrix is missing sample data."));
         throw e;
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
   if ( !_cmx || !_output )
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
   if ( _csm )
   {
      _csmPair = CSMatrix::Pair(_csm);
      preparePValueFilter();
      prepareRSquareFilter();
   }

   // initialize output file stream
   _stream.setDevice(_output);
   _stream.setRealNumberPrecision(8);
}






/*!
 * Prepares the PValue filter for the csm.
 */
void Extract::preparePValueFilter()
{
    bool ok = false;
    if ( _csmPValueFilter != "" )
    {
        QStringList filters = _csmPValueFilter.split("::");
        for ( int i = 0; i < filters.size(); i++ )
        {
            QStringList data = filters.at(i).split(",");
            data.at(0).toFloat(&ok);
            if (ok)
            {
                _csmPValueFilterThresh.append(data.at(0).toFloat());
                break;
            }
            if(data.size() != 3 && !ok)
            {
                E_MAKE_EXCEPTION(e);
                e.setTitle(tr("Invalid Input"));
                e.setDetails(tr("Invalid filter name given."));
                throw e;
            }
            else
            {
                _csmPValueFilterThresh.append(data.at(2).toFloat());
                _csmPValueFilterFeatureNames.append(data.at(0));
                _csmPValueFilterLabelNames.append(data.at(1));
            }
        }
    }
}

/*!
 * Prepares the RSquared filter for the csm.
 */
void Extract::prepareRSquareFilter()
{
    bool ok = false;
    if ( _csmRSquareFilter != "" )
    {
        QStringList filters = _csmRSquareFilter.split("::");
        for ( int i = 0; i < filters.size(); i++ )
        {
            QStringList data = filters.at(i).split(",");
            data.at(0).toFloat(&ok);
            if (ok)
            {
                _csmRSquareFilterThresh.append(data.at(0).toFloat());
                break;
            }
            if(data.size() != 2 && !ok)
            {
                E_MAKE_EXCEPTION(e);
                e.setTitle(tr("Invalid Input"));
                e.setDetails(tr("Invalid filter name given."));
                throw e;
            }
            else
            {
                _csmRSquareFilterThresh.append(data.at(2).toFloat());
                _csmRSquareFilterFeatureNames.append(data.at(0));
                _csmRSquareFilterLabelNames.append(data.at(1));
            }
        }
    }
}





/*!
 * Used to filter the given data by the name of the label and the pvalue.
 *
 * @param labelName The test name for the label.
 *
 * @param pValue The pValue assosiated with the test.
 *
 * @return True if the test should be included, false otherwise.
 */
bool Extract::PValuefilter(QString labelName, float pValue)
{
    if ( _csmPValueFilter != "" )
    {
        // If there are feature-specific filters then apply the filter to the respective field.
        if(_csmPValueFilterFeatureNames.size() != 0)
        {
            auto names = labelName.split("__");
            for ( int i = 0; i < _csmPValueFilterFeatureNames.size(); i++ )
            {
                if ( names.at(0) == _csmPValueFilterFeatureNames.at(i) && names.at(1) == _csmPValueFilterLabelNames.at(i) )
                {
                    if ( pValue > _csmPValueFilterThresh.at(i) )
                    {
                       return false;
                    }
                    else
                    {
                        return true;
                    }
                }
            }
        }
        // If there are no names for the filter, then check any field.
        else
        {
            if(pValue > _csmPValueFilterThresh.at(0))
            {
                return false;
            }
            else
            {
                return true;
            }
        }
    }
    return true;
}

/*!
 * Filters the given data by the name of the label and the pvalue.
 *
 * @param labelName The test name for the label.
 *
 * @param pValue The pValue assosiated with the test.
 *
 * @return True if the test should be included, false otherwise.
 */
bool Extract::RSquarefilter(QString labelName, float rSquared)
{
    if ( _csmRSquareFilter != "" )
    {
        if(_csmRSquareFilterFeatureNames.size() != 0)
        {
            auto names = labelName.split("__");
            for ( int i = 0; i < _csmRSquareFilterFeatureNames.size(); i++ )
            {
                if ( names.at(0) == _csmRSquareFilterFeatureNames.at(i) && names.at(1) == _csmRSquareFilterLabelNames.at(i) )
                {
                    if ( rSquared > _csmRSquareFilterThresh.at(i) )
                    {
                       return false;
                    }
                    else
                    {
                        return true;
                    }
                }
            }
        }
        else
        {
            if(rSquared < _csmRSquareFilterThresh.at(0))
            {
                return false;
            }
            else
            {
                return true;
            }
        }
    }
    return true;
}




/*!
 * Checks to make sure the filter names are in the tests names. If they
 * are not, it throws an error.
 *
 * @return True if the name apears somewhere in the test names, false
 *         otherwise.
 */
bool Extract::pValueFilterCheck()
{
    //No filter
    if ( _csmPValueFilter == "" )
    {
        return true;
    }

    //default filter all
    if(_csmPValueFilterFeatureNames.size() == 0 && _csmPValueFilterThresh.size() != 0)
    {
        return true;
    }

    //specific filter
    for ( int i = 0; i < _csm->getTestCount(); i++ )
    {
        auto names = _csm->getTestName(i).split("__");
        for ( int j = 0; j < _csmPValueFilterFeatureNames.size(); j++ )
        {
            if ( names.at(0) == _csmPValueFilterFeatureNames.at(j) && names.at(1) == _csmPValueFilterLabelNames.at(j) )
            {
                return true;
            }
        }
    }
    E_MAKE_EXCEPTION(e);
    e.setTitle(tr("Invalid Input"));
    e.setDetails(tr("Invalid filter name given."));
    throw e;
}


/*!
 * Checks to make sure the filter names are in the tests names. If they
 * are not, it throws an error.
 *
 * @return True if the name apears somewhere in the test names, false
 *         otherwise.
 */
bool Extract::rSquareFilterCheck()
{
    //No filter
    if ( _csmRSquareFilter == "" )
    {
        return true;
    }

    //default filter all
    if(_csmRSquareFilterFeatureNames.size() == 0 && _csmRSquareFilterThresh.size() != 0)
    {
        return true;
    }

    //specific filter
    for ( int i = 0; i < _csm->getTestCount(); i++ )
    {
        auto names = _csm->getTestName(i).split("__");
        for ( int j = 0; j < _csmRSquareFilterFeatureNames.size(); j++ )
        {
            if ( names.at(0) == _csmRSquareFilterFeatureNames.at(j) && names.at(1) == _csmRSquareFilterLabelNames.at(j) )
            {
                return true;
            }
        }
    }
    E_MAKE_EXCEPTION(e);
    e.setTitle(tr("Invalid Input"));
    e.setDetails(tr("Invalid filter name given."));
    throw e;
}
