#include "extract.h"
#include "extract_input.h"

using namespace std;



/*!
 * Return the total number of blocks this analytic must process as steps
 * or blocks of work. This implementation uses a work block for writing
 * each pair to the output file.
 */
int Extract::size() const
{
    EDEBUG_FUNC(this);

    return static_cast<int>(_cmx->size());
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

    // Each time this function is called we will read the next cluster pair.
    int clusterSize = _networkWriter->readNext();

    // Write clusters to the output file if they pass filters.
    for ( int k = 0; k < clusterSize; k++ )
    {
        // If the cluster passed all of the tests, then write it to the file.
        QVector<QString> passed = filterEdge(k);
        if (passed.size() > 0)
        {
            _networkWriter->writeEdgeCluster(k, passed);
        }
    }

    // If we're at the last element then finish up the file.
    if ( result->index() == size() - 1 )
    {
        _networkWriter->finish();
    }

    _networkWriter->checkStatus();
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
        E_MAKE_EXCEPTION(e)
        e.setTitle(tr("Invalid Argument"));
        e.setDetails(tr("Did not get valid input and/or output arguments."));
        throw e;
    }

    if ( _outputFormat != OutputFormat::Minimal && !_ccm )
    {
        E_MAKE_EXCEPTION(e)
        e.setTitle(tr("Invalid Argument"));
        e.setDetails(tr("--ccm is required for all output formats except minimal."));
        throw e;
    }

    // Set the proper network output class.
    switch ( _outputFormat )
    {
    case OutputFormat::Text:
        _networkWriter = new FullNetworkWriter(_emx, _cmx, _ccm, _csm, _output);
        break;
    case OutputFormat::Minimal:
        _networkWriter = new MinimalNetworkWriter(_emx, _cmx, _ccm, _csm, _output);
        break;
    case OutputFormat::GraphML:
        _networkWriter = new GMLNetworkWriter(_emx, _cmx, _ccm, _csm, _output);
        break;
    case OutputFormat::Tidy:
        _networkWriter = new TidyNetworkWriter(_emx, _cmx, _ccm, _csm, _output);
        break;
    }

    // Get the network writer object started.
    _networkWriter->initialize();

    // Save the test name columns that will be used by the
    // network writer so we don't have to keep looking them up.
    _testNames = _networkWriter->getTestNames();

    // Set any filters that the user requested. If the filters
    // are incorrectly set this function should throw an error.
    parseFilters(_csmPValueFilter, "pVal");
    parseFilters(_csmRSquareFilter, "rSqr");
}



/*!
 * Parse the filters from a filter string.
 *
 * @param filterString
 * @param type
 */
void Extract::parseFilters(QString filterString, QString type)
{
    // If the user provided no filters then just return;
    if (filterString == "")
    {
        return;
    }

    bool ok = false;
    bool failure = true;

    // If the comparision is unspecified we should set a
    // default based on the type of filter.
    QString defaultComp = (type == "rSqr")
        ? "gt"
        : "lt";

    QStringList filters = filterString.split("::");
    for ( auto& filter : filters )
    {
        QStringList tokens = filter.split(",");

        // Case #1: the user provided a single global threshold (e.g. 1e-3)
        if (tokens.size() == 1)
        {
            tokens.at(0).toFloat(&ok);

            if (ok)
            {
                QPair<QString, float> fpair;
                fpair.first = defaultComp;
                fpair.second = tokens.at(0).toFloat();
                _filters.insert(type, fpair);
                failure = false;
            }
        }

        // Case #2 the user provided two values. There are two cases.
        else if (tokens.size() == 2)
        {
            tokens.at(1).toFloat(&ok);

            if (ok)
            {
                // Case 2a:  global setting. e.g.: gt,1e-3
                if (tokens.at(0) == "gt" || tokens.at(0) == "lt")
                {
                    QPair<QString, float> fpair;
                    fpair.first = tokens.at(0);
                    fpair.second = tokens.at(1).toFloat();
                    QString testName = type;
                    _filters.insert(testName, fpair);
                    failure = false;
                }

                // Case 2b: e.g.:  Subspecies,1e-3
                else
                {
                    QPair<QString, float> fpair;
                    fpair.first = defaultComp;
                    fpair.second = tokens.at(1).toFloat();
                    // If this is a filter on a categorical field then it means
                    // that the user did not specify a category and they want
                    // to apply the filter on any category.  We need to
                    // iterate through the tests and insert a filter for each
                    // one that matches.
                    QVectorIterator<QString> testNamesIter(_testNames);
                    while (testNamesIter.hasNext())
                    {
                        QString testName = testNamesIter.next();
                        if (testName.contains(tokens.at(0) + "__"))
                        {
                            _filters.insert(testName, fpair);
                        }
                    }
                    // If this isn't a categorical field then we can just
                    // insert the filter as is.
                    QString testName = tokens.at(0) + "_" + type;
                    _filters.insert(testName, fpair);
                    failure = false;
                }
            }
        }

        // Case #3: the user provided three values
        else if (tokens.size() == 3)
        {
            tokens.at(2).toFloat(&ok);

            if (ok)
            {
                // Case 3a: Subspecies,gt,1e-3
                if (tokens.at(1) == "gt" || tokens.at(1) == "lt")
                {
                    QPair<QString, float> fpair;
                    fpair.first = tokens.at(1);
                    fpair.second = tokens.at(2).toFloat();
                    // If this is a filter on a categorical field then it means
                    // that the user did not specify a category and they want
                    // to apply the filter on any category.  We need to
                    // iterate through the tests and insert a filter for each
                    // one that matches.
                    QVectorIterator<QString> testNamesIter(_testNames);
                    while (testNamesIter.hasNext())
                    {
                        QString testName = testNamesIter.next();
                        if (testName.contains(tokens.at(0) + "__"))
                        {
                            _filters.insert(testName, fpair);
                        }
                    }
                    // If this isn't a categorical field then we can just
                    // insert the filter as is.
                    QString testName = tokens.at(0) + "_" + type;
                    _filters.insert(testName, fpair);
                    failure = false;
                }

                // Case 3b: Subspecies,Janpoica,1e-3
                else
                {
                    QPair<QString, float> fpair;
                    fpair.first = defaultComp;
                    fpair.second = tokens.at(2).toFloat();
                    QString testName = tokens.at(0) + "__" + tokens.at(1) + "_" + type;
                    _filters.insert(testName, fpair);
                    failure = false;
                }
            }
        }

        // Case #4: the user provided four values (e.g. Subspecies,Japonica,lt,1e-3)
        else if (tokens.size() == 4)
        {
            tokens.at(3).toFloat(&ok);

            if (ok && (tokens.at(2) == "gt" || tokens.at(2) == "lt"))
            {
                QPair<QString, float> fpair;
                fpair.first = tokens.at(2);
                fpair.second = tokens.at(3).toFloat();
                QString testName = tokens.at(0) + "__" + tokens.at(1) + "_" + type;
                _filters.insert(testName, fpair);
                failure = false;
            }
        }

        if (failure)
        {
            E_MAKE_EXCEPTION(e)
            e.setTitle(tr("Invalid P-Value Filter Input"));
            e.setDetails(tr("Invalid P-Value Filter arguments given."));
            throw e;
        }
    }
}



/*!
 * Performs filtering of a cluster of the current edge.
 *
 * @param k
 */
QVector<QString> Extract::filterEdge(int k)
{
    // Stores the tests that passed and failed
    QVector<QString> passed;
    QVector<QString> failed;

    // Exclude cluster if correlation is not within thresholds.
    float correlation = _networkWriter->getEdgeSimilarity(k);
    if ( fabs(correlation) < _minCorrelation || _maxCorrelation < fabs(correlation) )
    {
        return passed;
    }

    // If no filters are present then just return as we've passed
    // the correlation threshold test.
    if (_testNames.size() == 0)
    {
        passed.append("th");
        return passed;
    }

    // Iterate through any filters and check those. If any pass we keep the edge.
    for (auto& testName : _testNames)
    {
        // testName has a _pVal or _rSqr suffix, whereas the
        // official (or real) test name does not.
        QString realTestName = testName;
        realTestName.replace("_pVal", "")
            .replace("_RSqr", "");
        QPair<QString, float> filter;

        // First, check if this is a global pValue setting. If there
        // is a global p-value filter then we'll ignore any other
        // pvalue filter specified.
        if (testName.contains("_pVal") && _filters.contains("pVal"))
        {
            filter = _filters.find("pVal").value();
        }

        // Second, check if this is a global rSqr setting. If there
        // is a global r-squared filter then we'll ignore any other
        // r-squared filter specified.
        else if (testName.contains("_RSqr") && _filters.contains("rSqr"))
        {
            filter = _filters.find("rSqr").value();
        }

        // Third, check if there is a test-specific setting.
        else if (_filters.contains(testName))
        {
            filter = _filters.find(testName).value();
        }

        // If there are no filters then skip to the next one.
        else
        {
            continue;
        }

        // Now perform the filter check.
        float filter_value = filter.second;
        float test_value = static_cast<float>(_networkWriter->getEdgeTestValue(k, testName));

        if (filter.first == "lt" && test_value < filter_value)
        {
            passed.append(realTestName);
        }
        else if (filter.first == "gt" && test_value > filter_value)
        {
            passed.append(realTestName);
        }
        else
        {
            failed.append(realTestName);
        }
    }

    // If the user requested a filter on both the pVal and the RSqr
    // then we need to enforce a logical "and" and we need
    // to check that both tests passed (not just one).
    for (auto& testName : _testNames)
    {
        QString realTestName = testName;
        realTestName.replace("_pVal", "")
            .replace("_RSqr", "");

        if (passed.count(realTestName) == 1 && failed.count(realTestName) == 1)
        {
            passed.removeAll(realTestName);
        }
    }

    // Remove any duplicates from the passed list before returning.
    for (auto& testName : _testNames)
    {
        QString realTestName = testName;
        realTestName.replace("_pVal", "")
            .replace("_RSqr", "");

        if (passed.count(realTestName) > 1)
        {
            passed.removeAll(realTestName);
            passed.append(realTestName);
        }
    }

    // If any of the tests passed then keep this edge.
    if (passed.size() > 0)
    {
        return passed;
    }

    // If we're here then we did not pass any of the filters.
    return passed;
}
