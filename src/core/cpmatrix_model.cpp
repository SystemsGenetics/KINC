#include "cpmatrix_model.h"
#include "cpmatrix_pair.h"



using namespace std;



/*!
 * Construct a table model for a cluster matrix.
 *
 * @param matrix
 */
CPMatrix::Model::Model(CPMatrix* matrix):
    _matrix(matrix)
{
    EDEBUG_FUNC(this,matrix);

    setParent(matrix);
}



/*!
 * Return a header name for the table model using a given index.
 *
 * @param section
 * @param orientation
 * @param role
 */
QVariant CPMatrix::Model::headerData(int section, Qt::Orientation orientation, int role) const
{
    EDEBUG_FUNC(this,section,orientation,role);

    // orientation is not used
    Q_UNUSED(orientation);

    // if role is not display return nothing
    if ( role != Qt::DisplayRole )
    {
        return QVariant();
    }

    // get gene names
    EMetaArray geneNames {_matrix->geneNames()};

    // make sure section is within limits of gene name array
    if ( section >= 0 && section < geneNames.size() )
    {
        // return gene name
        return geneNames.at(section).toString();
    }

    // no gene found return nothing
    return QVariant();
}



/*!
 * Return the number of rows in the table model.
 *
 * @param index
 */
int CPMatrix::Model::rowCount(const QModelIndex&) const
{
    EDEBUG_FUNC(this);

    return _matrix->geneSize();
}



/*!
 * Return the number of columns in the table model.
 *
 * @param index
 */
int CPMatrix::Model::columnCount(const QModelIndex&) const
{
    EDEBUG_FUNC(this);

    return _matrix->geneSize();
}



/*!
 * Return a data element in the table model using the given index.
 *
 * @param index
 * @param role
 */
QVariant CPMatrix::Model::data(const QModelIndex& index, int role) const
{
    EDEBUG_FUNC(this,&index,role);

    // if role is not display return nothing
    if ( role != Qt::DisplayRole )
    {
        return QVariant();
    }

    // if row and column are equal return empty string
    if ( index.row() == index.column() )
    {
        return "";
    }

    // get constant pair and read in values
    const Pair pair(_matrix);
    int x {index.row()};
    int y {index.column()};
    if ( y > x )
    {
        swap(x,y);
    }
    pair.read({x,y});

    // Return value of pair as a string
    return pair.toString();
}
