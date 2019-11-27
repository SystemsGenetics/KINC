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
    _cmxPair.readNext();
    _netWriter->setPair(_cmxPair.index());

    // Write clusters to the output file if they pass filters.
    for ( int k = 0; k < _cmxPair.clusterSize(); k++ )
    {
        // If the cluster passed all of the tests, then write it to the file.
        if (filterEdge(k)) {
            _netWriter->writeEdgeCluster(k);
        }
    }

    // If we're at the last element then finish up the file.
    if ( result->index() == size() - 1 ) {
        _netWriter->finish();
    }

    // make sure writing output file worked
    if ( _stream.status() != QTextStream::Ok )
    {
        E_MAKE_EXCEPTION(e)
        e.setTitle(tr("File IO Error"));
        e.setDetails(tr("Qt Text Stream encountered an unknown error."));
        throw e;
    }
}

bool Extract::filterEdge(int k)
{
    float correlation = _netWriter->getEdgeSimilarity(k);

    // Exclude cluster if correlation is not within thresholds.
    if ( fabs(correlation) < _minCorrelation || _maxCorrelation < fabs(correlation) )
    {
        return false;
    }

    // Iterate through any filters and check those. If any pass we keep the edge.
    int num_filters = 0;
    bool keep = false;
    for (int i = 0; i < _testNames.size(); i++)
    {
        QString test_name = _testNames[i];
        QPair<QString, float> filter;

        // First check if there is a global pValue setting.
        if (test_name.contains("_pVal") && _filters.contains("pVal")) {
            filter = _filters.find("pVal").value();
        }
        // Second check if there is a global rSqr setting.
        else if (test_name.contains("_rSqr") && _filters.contains("rSqr")) {
            filter = _filters.find("rSqr").value();
        }
        // Third check if there is a test-specific setting.
        else if (_filters.contains(test_name)) {
            filter = _filters.find(test_name).value();
        }
        // If there are no filters then skip to the next one.
        else {
            continue;
        }
        num_filters++;

        // Now perform the filter check.
        float filter_value = filter.second;
        float test_value = static_cast<float>(_netWriter->getEdgeTestValue(k, i));
        if (filter.first == "lt" && test_value < filter_value)
        {
            keep = true;
        }
        if (filter.first == "gt" && test_value > filter_value)
        {
            keep = true;
        }
    }

    // If any of the filters passed then keep this edge.
    if (keep == true) {
        return true;
    }
    // If no filters passed but we had no filters then keep this edge
    // beause we've already passed the min/max correlation filter above.
    if (keep == false && num_filters == 0) {
        return true;
    }

    // If we're here then we did not pass any of the filters.
    return false;
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

    // initialize pairwise iterators
    _ccmPair = CCMatrix::Pair(_ccm);
    _cmxPair = CorrelationMatrix::Pair(_cmx);

    // initialize output file stream
    _stream.setDevice(_output);
    _stream.setRealNumberPrecision(8);

    // Set the proper network output class.
    switch ( _outputFormat )
    {
    case OutputFormat::Text:
        _netWriter = new FullNetWriter(&_stream, _emx, _cmx, _ccm, _csm);
        break;
    case OutputFormat::Minimal:
        _netWriter = new MinimalNetWriter(&_stream, _emx, _cmx, _ccm, _csm);
        break;
    case OutputFormat::GraphML:
        _netWriter = new GMLNetWriter(&_stream, _emx, _cmx, _ccm, _csm);
        break;
    }

    // Get the network writer object started.
    _netWriter->initialize();

    // Save the test name columns that will be used by the
    // network writer so we don't have to keep looking them up.
    _testNames = _netWriter->getTestNames();

    // Set any filters that the user requested. If the filters
    // are incorrectly set this function should throw an error.
    setFilters(_csmPValueFilter, "pVal");
    setFilters(_csmRSquareFilter, "rSqr");
}

/*!
 * \brief Extract::setFilters
 * \param input_filters
 * \param type
 */
void Extract::setFilters(QString input_filters, QString type) {

    bool ok = false;
    bool failure = true;

    // If the user provided no filters then just return;
    if (input_filters == "")
    {
        return;
    }


    // If the comparision is unspecified we should set a
    // default based on the type of filter.
    QString defaultComp {""};
    if (type == "rSqr")
    {
         defaultComp = "gt";
    }
    else
    {
        defaultComp = "lt";
    }

    QStringList filters = input_filters.split("::");
    for ( int i = 0; i < filters.size(); i++ )
    {
        QStringList data = filters.at(i).split(",");

        // Case #1: the user provided a single global threshold (e.g. 1e-3)
        if (data.size() == 1)
        {
           data.at(0).toFloat(&ok);
           if (ok)
           {
               QPair<QString, float> fpair;
               fpair.first = defaultComp;
               fpair.second = data.at(0).toFloat();
               _filters.insert(type, fpair);
               failure = false;
           }
        }

        // Case #2 the user provided two values. There are two cases.
        else if (data.size() == 2)
        {
           data.at(1).toFloat(&ok);
           if (ok) {

               // Case 2a:  global setting. e.g.: gt,1e-3
               if (data.at(0) == "gt" || data.at(0) == "lt")
               {
                   QPair<QString, float> fpair;
                   fpair.first = data.at(0);
                   fpair.second = data.at(1).toFloat();
                   QString test_name = type;
                   _filters.insert(test_name, fpair);
                   failure = false;
               }
               // Case 2b: e.g.:  Subspecies,1e-3
               else
               {
                   QPair<QString, float> fpair;
                   fpair.first = defaultComp;
                   fpair.second = data.at(1).toFloat();
                   QString test_name = data.at(0) + "_" + type;
                   _filters.insert(test_name, fpair);
                   failure = false;
               }
           }
        }

        // Case #3: the user provided three values
        else if (data.size() == 3) {
           data.at(2).toFloat(&ok);
           if (ok)
           {
               // Case 3a: Subspecies,gt,1e-3
               if (data.at(1) == "gt" || data.at(1) == "lt")
               {
                   QPair<QString, float> fpair;
                   fpair.first = data.at(1);
                   fpair.second = data.at(2).toFloat();
                   QString test_name = data.at(0) + "_" + type;
                   _filters.insert(test_name, fpair);
                   failure = false;
               }
               // Case 3b: Subspecies,Janpoica,1e-3
               else
               {
                   QPair<QString, float> fpair;
                   fpair.first = defaultComp;
                   fpair.second = data.at(2).toFloat();
                   QString test_name = data.at(0) + "__" + data.at(1) + "_" + type;
                   _filters.insert(test_name, fpair);
                   failure = false;
               }
           }
        }

        // Case #4: the user provided four values (e.g. Subspecies,Japonica,lt,1e-3)
        else if (data.size() == 4)
        {
           data.at(3).toFloat(&ok);
           if (ok && (data.at(2) == "gt" || data.at(2) == "lt"))
           {
               QPair<QString, float> fpair;
               fpair.first = data.at(2);
               fpair.second = data.at(3).toFloat();
               QString test_name = data.at(0) + "__" + data.at(1) + "_" + type;
               _filters.insert(test_name, fpair);
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

