#include "expressionmatrix_gene.h"
//






/*!
 *
 * @param index  
 */
float& ExpressionMatrix::Gene::operator[](int index)
{
   // 
   if ( index < 0 || index >= _matrix->_sampleSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Attempting to access gene expression %1 when maximum is %2.").arg(index)
                   .arg(_matrix->_sampleSize-1));
      throw e;
   }

   // 
   return _expressions[index];
}






/*!
 *
 * @param matrix  
 *
 * @param isInitialized  
 */
ExpressionMatrix::Gene::Gene(ExpressionMatrix* matrix, bool isInitialized):
   _matrix(matrix),
   _expressions(new float[matrix->sampleSize()])
{
   if ( isInitialized )
   {
      read(_index);
   }
}






/*!
 */
ExpressionMatrix::Gene::~Gene()
{
   delete[] _expressions;
}






/*!
 *
 * @param index  
 */
void ExpressionMatrix::Gene::read(int index)
{
   if ( index < 0 || index >= _matrix->_geneSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Attempting to read gene %1 when maximum is %2.").arg(index)
                   .arg(_matrix->_geneSize-1));
      throw e;
   }
   _matrix->seekExpression(index,0);
   for (int i = 0; i < _matrix->sampleSize() ;++i)
   {
      _matrix->stream() >> _expressions[i];
   }
   _index = index;
}






/*!
 */
bool ExpressionMatrix::Gene::readNext()
{
   if ( (_index + 1 ) >= _matrix->_geneSize )
   {
      return false;
   }
   read(_index + 1);
   return true;
}






/*!
 *
 * @param index  
 */
void ExpressionMatrix::Gene::write(int index)
{
   if ( index < 0 || index >= _matrix->_geneSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Attempting to write gene %1 when maximum is %2.").arg(index)
                   .arg(_matrix->_geneSize-1));
      throw e;
   }
   _matrix->seekExpression(index,0);
   for (int i = 0; i < _matrix->sampleSize() ;++i)
   {
      _matrix->stream() << _expressions[i];
   }
   _index = index;
}






/*!
 */
bool ExpressionMatrix::Gene::writeNext()
{
   if ( (_index + 1 ) >= _matrix->_geneSize )
   {
      return false;
   }
   write(_index + 1);
   return true;
}






/*!
 *
 * @param index  
 */
float ExpressionMatrix::Gene::at(int index) const
{
   // 
   if ( index < 0 || index >= _matrix->_sampleSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Attempting to access gene expression %1 when maximum is %2.").arg(index)
                   .arg(_matrix->_sampleSize-1));
      throw e;
   }

   // 
   return _expressions[index];
}
