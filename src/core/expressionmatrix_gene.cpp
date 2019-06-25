#include "expressionmatrix_gene.h"
//






/*!
 * Return the expression value at the given index.
 *
 * @param index
 */
float& ExpressionMatrix::Gene::operator[](int index)
{
   EDEBUG_FUNC(this,index);

   // make sure the index is valid
   if ( index < 0 || index >= _matrix->_sampleSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Attempting to access gene expression %1 when maximum is %2.").arg(index)
                   .arg(_matrix->_sampleSize-1));
      throw e;
   }

   // return the specified value
   return _expressions[index];
}






/*!
 * Construct a gene iterator for an expression matrix. Additionally, if the
 * matrix is already initialized, read the first gene into memory.
 *
 * @param matrix
 * @param isInitialized
 */
ExpressionMatrix::Gene::Gene(ExpressionMatrix* matrix, bool isInitialized):
   _matrix(matrix),
   _expressions(new float[matrix->sampleSize()])
{
   EDEBUG_FUNC(this,matrix,isInitialized);

   if ( isInitialized )
   {
      read(_index);
   }
}






/*!
 * Destruct a gene iterator.
 */
ExpressionMatrix::Gene::~Gene()
{
   EDEBUG_FUNC(this);

   delete[] _expressions;
}






/*!
 * Read a row of the expression matrix from the data object file into memory.
 *
 * @param index
 */
void ExpressionMatrix::Gene::read(int index)
{
   EDEBUG_FUNC(this,index);

   // make sure the index is valid
   if ( index < 0 || index >= _matrix->_geneSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Attempting to read gene %1 when maximum is %2.").arg(index)
                   .arg(_matrix->_geneSize-1));
      throw e;
   }

   // seek to the beginning of the specified row in the data object file
   _matrix->seekExpression(index,0);

   // read the entire row into memory
   for ( int i = 0; i < _matrix->sampleSize(); ++i )
   {
      _matrix->stream() >> _expressions[i];
   }

   // set the iterator's current index
   _index = index;
}






/*!
 * Read the next row of the expression matrix into memory.
 */
bool ExpressionMatrix::Gene::readNext()
{
   EDEBUG_FUNC(this);

   // make sure that there is another row in the expression matrix
   if ( (_index + 1) >= _matrix->_geneSize )
   {
      return false;
   }

   // read the next row
   read(_index + 1);

   // return success
   return true;
}






/*!
 * Write the iterator's row data to the data object file corresponding to
 * the given row index.
 *
 * @param index
 */
void ExpressionMatrix::Gene::write(int index)
{
   EDEBUG_FUNC(this,index);

   // make sure the index is valid
   if ( index < 0 || index >= _matrix->_geneSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Attempting to write gene %1 when maximum is %2.").arg(index)
                   .arg(_matrix->_geneSize-1));
      throw e;
   }

   // seek to the beginning of the specified row in the data object file
   _matrix->seekExpression(index,0);

   // write the entire row to the data object
   for ( int i = 0; i < _matrix->sampleSize(); ++i )
   {
      _matrix->stream() << _expressions[i];
   }

   // set the iterator's current index
   _index = index;
}






/*!
 * Write the iterator's row data to the next row in the data object file.
 */
bool ExpressionMatrix::Gene::writeNext()
{
   EDEBUG_FUNC(this);

   // make sure there is another row in the expression matrix
   if ( (_index + 1) >= _matrix->_geneSize )
   {
      return false;
   }

   // write to the next row
   write(_index + 1);

   // return success
   return true;
}






/*!
 * Return the expression value at the given index.
 *
 * @param index
 */
float ExpressionMatrix::Gene::at(int index) const
{
   EDEBUG_FUNC(this,index);

   // make sure the index is valid
   if ( index < 0 || index >= _matrix->_sampleSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Attempting to access gene expression %1 when maximum is %2.").arg(index)
                   .arg(_matrix->_sampleSize-1));
      throw e;
   }

   // return the specified value
   return _expressions[index];
}
