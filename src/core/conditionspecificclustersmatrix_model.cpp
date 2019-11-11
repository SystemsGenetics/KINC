#include "conditionspecificclustersmatrix_model.h"
#include "conditionspecificclustersmatrix_pair.h"



/*!
 * Initializes a new CSM model object.
 *
 * @param matrix The parent matrix of the model.
 */
CSMatrix::Model::Model(CSMatrix* matrix) : _matrix(matrix)
{
    EDEBUG_FUNC(this,matrix);

    setParent(matrix);
}



/*!
 * Print the header for the tablular representation of the CSM model to the screen.
 *
 * @param section The column or row index dependig on the orientation type.
 *
 * @param orientation Which header the function is trying to map, vertical or horizontal.
 *
 * @param role The display role the function is querying for, here it only makes sense
 *             with a diaplay role.
 *
 * @return The header cell data at the queried location.
 */
QVariant CSMatrix::Model::headerData(int section, Qt::Orientation orientation, int role) const
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
 * Queries for the number of columns.
 *
 * @param An unused parameter.
 *
 * @return An integer representation of the column count.
 */
int CSMatrix::Model::columnCount(const QModelIndex&) const
{
    EDEBUG_FUNC(this);
    return _matrix->geneSize();
}



/*!
 * Queries for the number of row.
 *
 * @param An unused parameter.
 *
 * @return An integer representation of the row count.
 */
int CSMatrix::Model::rowCount(const QModelIndex&) const
{
    EDEBUG_FUNC(this);
    return _matrix->geneSize();
}



/*!
 * Print the data for the tablular representation of the CSM model to the screen.
 *
 * @param section The column or row index dependig on the orientation type.
 *
 * @param orientation Which header the function is trying to map, vertical or horizontal.
 *
 * @param role The display role the function is querying for, here it only makes sense
 *             with a diaplay role.
 *
 * @return The header cell data at the queried location.
 */
QVariant CSMatrix::Model::data(const QModelIndex& index, int role) const
{
    EDEBUG_FUNC(this, &index, role);

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
        std::swap(x,y);
    }
    pair.read({x,y});
    return pair.toString();
}
