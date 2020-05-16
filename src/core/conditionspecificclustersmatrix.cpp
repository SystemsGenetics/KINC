#include "conditionspecificclustersmatrix.h"
#include "conditionspecificclustersmatrix_model.h"



/*!
 * Writes the annotation matrix information into the metadata.
 *
 * @param features The names of all the major feature in the annotation matrix,
 *                 these where the column names.
 *
 * @param featureinfo The information about the feature, their testing type and subcatagory names.
 *
 * @param data All the data in the annotation matrix.
 */
void CSMatrix::initialize(const EMetaArray& features,
                          const QVector<EMetaArray>& featureInfo,
                          const QVector<EMetaArray>& data,
                          int numTests,
                          QString testNames,
                          const EMetaArray& geneNames,
                          int maxClusterSize,
                          int subheader)
{
    EDEBUG_FUNC(this,&features,&featureInfo,&data);

    // make sure all containers have something in them
    if ( features.isEmpty() || featureInfo.isEmpty() || data.isEmpty())
    {
        E_MAKE_EXCEPTION(e);
        e.setTitle(tr("Domain Error"));
        if ( features.isEmpty() )
        {
            e.setDetails(tr("No Features Provided"));
        }
        else if ( featureInfo.isEmpty() )
        {
            e.setDetails(tr("No Feature Info Provided"));
        }
        else if ( data.isEmpty() )
        {
            e.setDetails(tr("No Data Provided"));
        }
        throw e;
    }

    EMetaObject metaObject {meta().toObject()};
    EMetaObject type;
    EMetaArray labels;
    EMetaObject info;
    EMetaObject feature;

    // Set the feature information in the meta data
    for ( int i = 0; i < features.size(); i++ )
    {
        // Set the test type
        type.insert("Type", featureInfo.at(i).at(0));

        // Set the labels
        // Categorical
        if ( featureInfo.at(i).size() > 2 && featureInfo.at(i).at(0).toString() == "Categorical" )
        {
            for ( int j = 2; j < featureInfo.at(i).size(); j++ )
            {
                labels.append(featureInfo.at(i).at(j));
            }
            info.insert("Labels", labels);
        }
        // create the labels and test drop bar aka label information
        info.insert("Test", type);
        // pair label information to feature names
        feature.insert(features.at(i).toString(), info);
        // prepare for the next feature
        type.clear();
        labels.clear();
        info.clear();
    }
    // create the feature drop bar
    metaObject.insert("Features", feature);

    // add the sample names
    EMetaArray sampleNames;
    for ( auto sample : data.at(0) )
    {
        sampleNames.append(sample);
    }
    metaObject.insert("Samples", sampleNames);

    // add all the data for the annotation array, exept the sample names
    EMetaObject Data;
    for ( int i = 0; i < features.size(); i++ )
    {
        Data.insert(features.at(i).toString(), data.at(i));
    }
    metaObject.insert("Data", Data);

    // add num_tests
    _testcount = numTests;
    metaObject.insert("Number of Tests", numTests);

    EMetaArray namesOfTests;
    auto names = testNames.split(":", QString::SkipEmptyParts);
    for ( auto name : names )
    {
        namesOfTests.append(name);
    }
    metaObject.insert("Test Names", namesOfTests);

    // set the root of the meta data to the changed metaobject
    setMeta(metaObject);
    _sampleSize = sampleNames.size();

    Matrix::initialize(geneNames, maxClusterSize, sizeof(double) * 2 * _testcount, subheader);
}



/*!
 * Create a new table model representation of the data.
 *
 * @return Pointer to the table model representation of the data.
 */
QAbstractTableModel* CSMatrix::model()
{
    EDEBUG_FUNC(this);

    if ( !_model )
    {
        _model = new Model(this);
    }
    return _model;
}



/*!
 * Query for the sample size.
 *
 * @return Integer representation of the number of samples in the augmented matrix.
 */
int CSMatrix::sampleSize() const
{
    EDEBUG_FUNC(this);

    return _sampleSize;
}



/*!
 * Write header information.
 */
void CSMatrix::writeHeader()
{
    EDEBUG_FUNC(this);

    stream() << _sampleSize << _testcount;
}



/*!
 * Read header information.
 */
void CSMatrix::readHeader()
{
    EDEBUG_FUNC(this);

    stream() >> _sampleSize >> _testcount;
}



/*!
 * Get the test name information from the meta data.
 *
 * @param index The test name index.
 *
 * @return The test name.
 */
QString CSMatrix::getTestName(int index) const
{
    return meta().toObject().at("Test Names").toArray().at(index).toString();
}



QString CSMatrix::getTestType(int index) const
{
    auto names = getTestName(index).split("__");

    return meta().toObject().at("Features")
                 .toObject().at(names.at(0))
                 .toObject().at("Test")
                 .toObject().at("Type").toString();
}



/*!
 * Get the number of tests that were conducted.
 *
 * @return The number of tests.
 */
qint32 CSMatrix::getTestCount()
{
    return _testcount;
}
