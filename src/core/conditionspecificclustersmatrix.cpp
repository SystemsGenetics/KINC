#include "conditionspecificclustersmatrix.h"
#include "conditionspecificclustersmatrix_model.h"
//





/*!
*  Implements an interface that writes the annotation matrix information into the metadata.
*
* @param features The names of all the major feature in the annotation matrix,
*                 these where the column names.
*
* @param featureinfo The information about the feature, their testing type and subcatagory names.
*
* @param data All the data in the annotation matrix.
*/
void CSMatrix::initialize(const EMetaArray& features, const QVector<EMetaArray>& featureInfo, const QVector<EMetaArray>& data, int& numTests, QString testNames)
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
    int num_tests = 0;

    //set the feature information in the meta data
    for ( int i = 0; i < features.size(); i++ )
    {
        //set the test type
        type.insert("Type", featureInfo.at(i).at(0));

        //set the labels
        if ( featureInfo.at(i).size() > 2 )
        {
            for ( int j = 2; j < featureInfo.at(i).size(); j++ )
            {
                labels.append(featureInfo.at(i).at(j));
                if ( featureInfo.at(i).at(0).toString() != "None" && featureInfo.at(i).at(0).toString() != "Unknown" )
                {
                    ++num_tests;
                }
            }
            info.insert("Labels", labels);
        }
        else
        {
            if ( featureInfo.at(i).at(0).toString() != "None" && featureInfo.at(i).at(0).toString() != "Unknown" )
            {
                ++num_tests;
            }
        }
        //create the labels and test drop bar aka label information
        info.insert("Test", type);
        //pair label information to feature names
        feature.insert(features.at(i).toString(), info);
        //prepare for the next feature
        type.clear();
        labels.clear();
        info.clear();
    }
    //create the feature drop bar
    metaObject.insert("Features", feature);

    //add the sample names
    EMetaArray sampleNames;
    for ( auto sample : data.at(0) )
    {
        sampleNames.append(sample);
    }
    metaObject.insert("Samples", sampleNames);

    //add all the data for the annotation array, exept the sample names
    EMetaObject Data;
    for ( int i = 0; i < features.size(); i++ )
    {
        Data.insert(features.at(i).toString(), data.at(i));
    }
    metaObject.insert("Data", Data);

    //add num_tests
    numTests += num_tests;
    _testcount = numTests;
    metaObject.insert("Number of Tests", numTests);

    EMetaArray namesOfTests;
    auto names = testNames.split(":", QString::SkipEmptyParts);
    for ( auto name : names )
    {
        namesOfTests.append(name);
    }
    metaObject.insert("Test Names", namesOfTests);

    //set the root of the meta data to the changed metaobject
    setMeta(metaObject);
    _sampleSize = sampleNames.size();
}





/*!
*  Implements an interface that writes the pairwise specific metadata.
*
* @param geneNames Names of all the genes in the gene expresion matrix.
*
* @param maxClusterSize Maximum amount of pairs a cell can have in the pairwise matrix,
*                       system max is 64.
*
* @param subheader The size of the subheader, only contains the sample size.
*/
void CSMatrix::initialize(const EMetaArray& geneNames, int maxClusterSize, int subheader)
{
    EDEBUG_FUNC(this,&geneNames,maxClusterSize,subheader);
    Matrix::initialize(geneNames, maxClusterSize, sizeof(double) * 2 * _testcount, subheader);
}





/*!
*  Implements an interface to create a new table model representation of the data.
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
*  Implements an interface to query for the sample size.
*
* @return Integer representation of the number of samples in the augmented matrix.
*/
int CSMatrix::sampleSize() const
{
    EDEBUG_FUNC(this);
    return _sampleSize;
}





/*!
*  Implements an interface to write header information.
*/
void CSMatrix::writeHeader()
{
    EDEBUG_FUNC(this);
    stream() << _sampleSize << _testcount;
}





/*!
*  Implements an interface to read header information.
*/
void CSMatrix::readHeader()
{
    EDEBUG_FUNC(this);
    stream() >> _sampleSize >> _testcount;
}




/*!
*  Implements an interface to set the number of tests that are going to be ran.
*
* @param newData The number of tests you want to set the variable to.
*/
void CSMatrix::setTestCount(qint32 newData)
{
    EDEBUG_FUNC(this, newData);
    _testcount = newData;
}





/*!
*  Implements an interface to get the test name information from the meta data.
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
    QString test_name = getTestName(index);
    return meta().toObject().at("Features")
                 .toObject().at(test_name)
                 .toObject().at("Test")
                 .toObject().at("Type").toString();
}



/*!
*  Implements an interface to get the number of tests that were conducted.
*
* @return The number of tests.
*/
qint32 CSMatrix::getTestCount()
{
    return _testcount;
}
