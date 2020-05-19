#include "conditionaltest.h"
#include "conditionaltest_input.h"
#include "conditionaltest_resultblock.h"
#include "conditionaltest_workblock.h"
#include "conditionaltest_serial.h"
#include "conditionspecificclustersmatrix.h"
#include "conditionspecificclustersmatrix_pair.h"
#include "correlationmatrix_pair.h"
#include <ace/core/elog.h>
#include <ace/core/ace_qmpi.h>



/*!
 * Supplies ACE with the number of work blocks it is going to create.
 *
 * @return How many pieces you want to break up the task your working on.
 */
int ConditionalTest::size() const
{
    EDEBUG_FUNC(this);

    return (static_cast<qint64>(_cmx->size()) + _workBlockSize - 1) / _workBlockSize;
}



/*!
 * Creates a block of work at the given index.
 *
 * @param index The index at which the work block should be made.
 *
 * @return Pointer to the work block.
 */
std::unique_ptr<EAbstractAnalyticBlock> ConditionalTest::makeWork(int index) const
{
    EDEBUG_FUNC(this,index);

    if ( ELog::isActive() )
    {
        ELog() << tr("Making work index %1 of %2.\n").arg(index).arg(size());
    }

    // compute parameters for work block
    qint64 start {_workBlockStart};
    qint64 size {std::min(_cmx->size() - index * _workBlockSize, static_cast<qint64>(_workBlockSize))};

    // initialize pairwise iterator for cmx file
    CorrelationMatrix::Pair cmxPair(_cmx);

    // iterate to the start index of the next work block
    cmxPair.read(Pairwise::Index(start));

    for ( qint64 i = 0; i < size; i++ )
    {
        cmxPair.readNext();
    }

    // save start index of next work block
    _workBlockStart = cmxPair.index().toRawIndex();

    // construct work block
    return std::unique_ptr<EAbstractAnalyticBlock>(new WorkBlock(index, start, size));
}



/*!
 * create uninitialized work blocks.
 *
 * @return a pointer to an uninitialized work block
 */
std::unique_ptr<EAbstractAnalyticBlock> ConditionalTest::makeWork() const
{
    EDEBUG_FUNC(this);

    return std::unique_ptr<EAbstractAnalyticBlock>(new WorkBlock);
}



/*!
 * create uninitialized result blocks.
 *
 * @return a pointer to an uninitialized result block
 */
std::unique_ptr<EAbstractAnalyticBlock> ConditionalTest::makeResult() const
{
    EDEBUG_FUNC(this);

    return std::unique_ptr<EAbstractAnalyticBlock>(new ResultBlock);
}



/*!
 * An interface that processes result blocks once the serial is done working on
 * them.
 *
 * @param result The processed work block from the serial.
 */
void ConditionalTest::process(const EAbstractAnalyticBlock* result)
{
    EDEBUG_FUNC(this,result);

    if ( ELog::isActive() )
    {
        ELog() << tr("Processing result %1 of %2.\n").arg(result->index()).arg(size());
    }

    // Iterate through the result block pairs.
    const ResultBlock* resultBlock {result->cast<ResultBlock>()};

    for ( auto& pair : resultBlock->pairs() )
    {
        // Create pair objects for the output data file.
        CSMatrix::Pair csmPair(_out);

        // Iterate through the clusters in the pair.
        for ( int k = 0; k < pair.pValues.size(); ++k )
        {
            // add each cluster into the CSM
            csmPair.addCluster(1, _numTests);

            for ( int i = 0; i < pair.pValues.at(k).size(); ++i )
            {
                csmPair.at(k, i, "pvalue") = pair.pValues.at(k).at(i);
                csmPair.at(k, i, "r2") = pair.r2.at(k).at(i);
            }
        }

        // write the info into the CSM
        if ( csmPair.clusterSize() > 0 )
        {
            csmPair.write(pair.index);
        }
    }
}



/*!
 * An interface to create a new input data object.
 *
 * @return Pointer to the new input data object.
 */
EAbstractAnalyticInput* ConditionalTest::makeInput()
{
    EDEBUG_FUNC(this);

    return new Input(this);
}



/*!
 * Create a new serial object.
 *
 * @return Pointer to a serial new object.
 */
EAbstractAnalyticSerial* ConditionalTest::makeSerial()
{
    EDEBUG_FUNC(this);

    return new Serial(this);
}



/*!
 * An interface to initialize the analytic process, it reads in the line number
 * of the annotaion matrix, so we can read in the data properly.
 */
void ConditionalTest::initialize()
{
    EDEBUG_FUNC(this);

    // only the master process needs to validate arguments
    auto& mpi {Ace::QMPI::instance()};
    if ( mpi.isMaster() )
    {

        // make sure input data is valid
        if ( !_ccm )
        {
            E_MAKE_EXCEPTION(e);
            e.setTitle(tr("Invalid Argument"));
            e.setDetails(tr("Did not get valid CCM data argument."));
            throw e;
        }

        if ( !_cmx )
        {
            E_MAKE_EXCEPTION(e);
            e.setTitle(tr("Invalid Argument"));
            e.setDetails(tr("Did not get valid CMX data argument."));
            throw e;
        }

        if ( !_amx )
        {
            E_MAKE_EXCEPTION(e);
            e.setTitle(tr("Invalid Argument"));
            e.setDetails(tr("Did not get valid AMX data argument."));
            throw e;
        }

        if (  !_emx )
        {
            E_MAKE_EXCEPTION(e);
            e.setTitle(tr("Invalid Argument"));
            e.setDetails(tr("Did not get valid EMX data argument."));
            throw e;
        }

        // Initialize work block size.
        if ( _workBlockSize == 0 )
        {
            int numWorkers = std::max(1, mpi.size() - 1);
            _workBlockSize = std::min(32768LL, _cmx->size() / numWorkers);
        }
    }

    // open the stream to the correct file.
    _stream.setDevice(_amx);

    // atain number of lines in the file.
    while ( !_stream.atEnd() )
    {
        _stream.readLine();
        _amxNumLines++;
    }

    // go back to the begginning
    _stream.seek(0);

    // make sure reading input file worked
    if ( _stream.status() != QTextStream::Ok )
    {
        E_MAKE_EXCEPTION(e);
        e.setTitle(tr("File IO Error"));
        e.setDetails(tr("Qt Text Stream encountered an unknown error."));
        throw e;
    }

    // Parse the user specified tests and test types. The test types provided by
    // the user will override the defaults.
    setUserTests();
    setUserTestTypes();

    // Read in the annotation matrix.
    // Set the stream to the right file.
    _stream.setDevice(_amx);
    _stream.seek(0);
    setFeatures();
    setTestTypes();
    setData();
    setNumTests();

    // Order the data read in from the AMX file by the
    // sample order in the EMX file.
    orderLabelsBySample();
}



/*!
 * Validate output arguments for this analytic.
 */
void ConditionalTest::initializeOutputs()
{
    EDEBUG_FUNC(this);

    if ( !_out )
    {
        E_MAKE_EXCEPTION(e);
        e.setTitle(tr("Invalid Argument"));
        e.setDetails(tr("The required output data object was not set."));
        throw e;
    }    

    // Needed by the CSM specific initializer
    EMetaArray features;
    QVector<EMetaArray> featureInfo;
    QVector<EMetaArray> Data;
    Data.resize(_data.size());
    featureInfo.resize(_features.size());

    // Meta Data for Features.
    for ( int i = 0; i < _features.size(); i++ )
    {
        features.append(_features.at(i).at(0));
    }

    // Meta Data for the test types of the catagories.
    for ( int i = 0; i < featureInfo.size(); i++ )
    {
        switch(_testType.at(i))
        {
            case CATEGORICAL : featureInfo[i].append(tr("Categorical"));
                break;
            case ORDINAL : featureInfo[i].append(tr("Ordinal"));
                break;
            case QUANTITATIVE : featureInfo[i].append(tr("Quantitative"));
                break;
            case NONE : featureInfo[i].append(tr("None"));
                break;
            default : featureInfo[i].append(tr("Unknown"));
                break;
        }
    }

    // Meta Data for the sub categories.
    for ( int i = 0; i < _features.size(); i++ )
    {
        for ( int j = 0; j < _features.at(i).size(); j++ )
        {
            featureInfo[i].append(_features.at(i).at(j));
        }
    }

    // Meta Data for all the data in the annotation matrix.
    for ( int i = 0; i < _data.size(); i++ )
    {
        for ( int j = 0; j < _data.at(i).size(); j++ )
        {
            Data[i].append(_data.at(i).at(j).toString());
        }
    }

    // inserts the Meta Data into the CSM Object stored on the disk.
    qint32 maxCluster = 64, subHeaders = 12;
    _out->initialize(features, featureInfo, Data, _numTests, testNames(),
                     _emx->geneNames(), maxCluster, subHeaders);
}



void ConditionalTest::setData()
{
    EDEBUG_FUNC(this);

    QString line;

    _data.resize(_features.size());

    for ( int i = 1; i < _amxNumLines; i++ )
    {
        line = _stream.readLine();
        // splits the file along the delemeter
        auto words2 = line.split(_delimiter, QString::KeepEmptyParts, Qt::CaseInsensitive);

        // add the data to our arrays
        for ( int j = 0; j < words2.size(); j++ )
        {
            _data[j].append(words2[j]);
            if ( _testType.at(j) == CATEGORICAL )
            {
                // this will add treatments types, leaf types, and any other types into the meta data
                if ( !_features.at(j).contains(words2[j]) )
                {
                    _features[j].append(words2[j]);
                }
            }
        }
    }
}



void ConditionalTest::setFeatures()
{
    EDEBUG_FUNC(this);

    // Read the header line from the input file.
    QString line = _stream.readLine();

    // The file can be a csv or a tab delimited AMX. The default is tab delimited.
    if ( _delimiter == "tab" )
    {
        _delimiter = "\t";
    }

    // Splits the line using the delimiter.
    auto headers = line.split(_delimiter, QString::SkipEmptyParts, Qt::CaseInsensitive);


    // Put the column names from the header into the amxdata array.
    for ( int i = 0; i < headers.size(); i++ )
    {
        _features.append(QVector<QString>());
        _features[i].append(headers[i]);
    }
}



void ConditionalTest::setNumTests()
{
    EDEBUG_FUNC(this);
    _numTests = 0;
    // Set the feature information in the meta data
    for ( int i = 0; i < _features.size(); i++ )
    {
        if ( _testType[i] == CATEGORICAL )
        {
            _numTests += _features[i].size() - 1;
        }
        // Ordinal and Quantitative
        else
        {
            if ( _testType[i] != NONE && _testType[i] != UNKNOWN )
            {
                ++_numTests;
            }
        }
    }
}



/*!
 * An interface to decide what the test types are going to be for each feature.
 *
 * @param dataTestType An array storing the test information for each feature.
 */
void ConditionalTest::setTestTypes()
{
    EDEBUG_FUNC(this);

    QString line;
    int num_lines = 0;

    // Counts here for an aproximation for what test type it should be.
    QVector<QVector<qint32>> counts(_features.size());

    _testType.resize(_features.size());
    for (int i = 0; i < _features.size(); i++) {
        _testType[i] = UNKNOWN;
        counts[i].resize(4);
    }

    while(!_stream.atEnd())
    {
        // read in the line
        line = _stream.readLine();
        num_lines++;

        // splits the file along the commas or tabs depending
        if ( _delimiter == "tab" )
        {
            _delimiter = "\t";
        }

        auto words = line.split(_delimiter, QString::SkipEmptyParts, Qt::CaseInsensitive);

        for ( int i = 0; i < words.size(); i++ )
        {
            // Check if the word is an integer or float.
            bool isFloat = 0;
            bool isInt = 0;
            words[i].toFloat(&isFloat);
            words[i].toInt(&isInt);

            // Count the missing values, ints, floats and strings.
            if ( QString::compare(words[i], _missing) == 0 ) {
                counts[i][UNKNOWN]++;
            }
            else if (isInt)
            {
                 counts[i][ORDINAL]++;
            }
            else if (isFloat)
            {
                counts[i][QUANTITATIVE]++;
            }
            else
            {
                 counts[i][CATEGORICAL]++;
            }
        }
    }

    for ( int i = 0; i < counts.size() ; i++ )
    {
        // If all of the values are ordinal (integer) then we'll consdier this
        // column ordinal.
        if ( counts.at(i).at(ORDINAL) + counts[i][UNKNOWN] == num_lines )
        {
            _testType[i] = ORDINAL;
        }
        // If all of the values are quantitative or ordinal (integer) or missing
        // then we'll consider this column quantitative.
        else if ( counts.at(i).at(QUANTITATIVE) + counts.at(i).at(ORDINAL) + counts[i][UNKNOWN] == num_lines )
        {
            _testType[i] = QUANTITATIVE;
        }
        // If all of the values are strings then it's categorical
        else if ( counts.at(i).at(CATEGORICAL) + counts[i][UNKNOWN]== num_lines )
        {
            _testType[i] = CATEGORICAL;
        }
        // If we're here then it means we have some combination of strings and
        // numbers. This column cannot be considered.
        else {
            _testType[i] = NONE;
        }
    }

    // We need to change the test types here if the user has overridden them.
    for ( int i = 0; i < _userTestTypes.size(); i++ )
    {
        for ( int j = 0; j < _testType.size(); j++ )
        {
            if ( _features.at(j).at(0) == _userTestTypes.at(i).at(0) )
            {
                if ( QString::compare(_userTestTypes.at(i).at(1), "CATEGORICAL", Qt::CaseInsensitive) == 0 )
                {
                    _testType[j] = CATEGORICAL;
                }
                else if ( QString::compare(_userTestTypes.at(i).at(1), "ORDINAL", Qt::CaseInsensitive) == 0 )
                {
                    _testType[j] = ORDINAL;
                }
                else if ( QString::compare(_userTestTypes.at(i).at(1), "QUANTITATIVE", Qt::CaseInsensitive) == 0 )
                {
                    _testType[j] = QUANTITATIVE;
                }
                else {
                    E_MAKE_EXCEPTION(e);
                    e.setTitle(tr("Unknown override type in the --feat_type argument."));
                    e.setDetails(tr("Unknown override type: %1.  Valid choices are: %2")
                        .arg(_userTestTypes.at(i).at(1))
                        .arg("categorical, quantitative, ordinal"));
                    throw e;
                }
            }
        }
    }

    // If the user has specified tests that need to be
    // run then set all others to NONE.
    for ( int j = 0; j < _features.size(); j++ )
    {
        int check = 0;
        for ( int k = 0; k < _userTests.size(); k++ )
        {
            if ( _features.at(j).at(0) == _userTests.at(k) )
            {
                check = 1;
            }
        }
        if ( check == 0 )
        {
            _testType[j] = NONE;
        }
    }

    // Reset the file back to the beginning for reparsing.
    _stream.seek(0);
    _stream.readLine();
}



/*!
 * An interface to find the max index of a vector.
 *
 * @param counts Vector to find the max index of.
 *
 * @return The index of the maximum balue in the vector.
 */
int ConditionalTest::max(QVector<qint32> &counts) const
{
    EDEBUG_FUNC(this,&counts);

    int maxnum = 0;

    if ( counts.size() > 0 )
    {
        qint32 max = counts.at(0);

        for ( int i = 1; i < counts.size(); i++ )
        {
            if ( counts.at(i) > max )
            {
                max = counts.at(i);
                maxnum = i;
            }
        }
    }

    return maxnum;
}



/*!
 * An interface to sperate the test out.
 */
void ConditionalTest::setUserTests()
{
    EDEBUG_FUNC(this);

    // make sure input data is valid
    _userTests = _userTestsStr.split(",", QString::SkipEmptyParts, Qt::CaseInsensitive).toVector();

    if ( _userTests.isEmpty() )
    {
        E_MAKE_EXCEPTION(e);
        e.setTitle(tr("Invalid Argument"));
        e.setDetails(tr("Please provide Features to Test"));
        throw e;
    }
}



/*!
 * An interface to override testing types.
 */
void ConditionalTest::setUserTestTypes()
{
    EDEBUG_FUNC(this);

    auto words = _userTestTypesStr.split(",", QString::SkipEmptyParts, Qt::CaseInsensitive).toVector();

    for ( int i = 0; i < words.size(); i++ )
    {
        _userTestTypes.append(QVector<QString>());

        auto word = words.at(i).split(":");

        for ( auto item : word )
        {
            _userTestTypes[i].append(item);
        }
    }
}



/*!
 * An interface to provide the names for the tests, creating an easier way to
 * label the data in the ouptut file.
 */
QString ConditionalTest::testNames()
{
    QString string;

    for ( int i = 0; i < _features.size(); i++ )
    {
        if ( _testType.at(i) == CATEGORICAL )
        {
            for ( int j = 1; j < _features.at(i).size(); j++ )
            {
                string += _features.at(i).at(0);
                string += "__";
                string += _features.at(i).at(j);
                string += ":";
            }
        }
        else if ( _testType.at(i) == QUANTITATIVE || _testType.at(i) == ORDINAL )
        {
            string += _features.at(i).at(0);
            string += ":";
        }
    }

    return string;
}



/*!
 * An interface to move the amx data around so that the emx samples and the
 * annotation matrix samples are in the same order
 */
void ConditionalTest::orderLabelsBySample()
{
    int sampleIndex = 0;
    int startIndex = 0;
    QVariant temp;

    // find the sample index, if none is found with the name "samples" then
    // we default to using the first column.
    for ( int i = 0; i< _features.size(); i++ )
    {
        if ( QString::compare(_features.at(i).at(0), "samples", Qt::CaseInsensitive) ||
             QString::compare(_features.at(i).at(0), "sample", Qt::CaseInsensitive) )
        {
            sampleIndex = i;
            break;
        }
    }

    // if the sample size differes break.
    if ( _data.at(sampleIndex).size() != _emx->sampleSize() )
    {
        E_MAKE_EXCEPTION(e);
        e.setTitle(tr("Sample Size Error"));
        e.setDetails(tr("Sample size in emx does not match annotation matrix."));
        throw e;
    }

    // look through all the samples, if you find one that doesnt match
    // switch it for the right sample.
    for ( int i = 0; i < _emx->sampleSize(); i++ )
    {
        if ( _emx->sampleNames().at(i).toString() != _data.at(sampleIndex).at(i) )
        {
            // find the right sample.
            for ( int j = startIndex; j <  _data.at(sampleIndex).size(); j++ )
            {
                // switch the samples.
                if ( _emx->sampleNames().at(i).toString() == _data.at(sampleIndex).at(j) )
                {
                    for ( int k = 0; k < _data.size(); k++ )
                    {
                        // move wrong sample to temp.
                        temp = _data.at(k).at(i);
                        // move right sample to the right place.
                        _data[k][i] = _data.at(k).at(j);
                        // replace the right sample with the wrong sample.
                        _data[k][j] = temp;
                    }
                    break;
                }

                if ( j == _emx->sampleSize() - 1 )
                {
                    E_MAKE_EXCEPTION(e);
                    e.setTitle(tr("Sample Size Error"));
                    e.setDetails(tr("Sample not in emx."));
                    throw e;
                }
            }

            // increment the start nidex so you dont have to start at 0 every time.
            startIndex++;
        }
    }
}
