#include "expressionmatrix_constgene.h"
//






/*!
 *
 * @param matrix  
 *
 * @param isInitialized  
 */
ExpressionMatrix::ConstGene::ConstGene(const ExpressionMatrix* matrix, bool isInitialized):
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
ExpressionMatrix::ConstGene::~ConstGene()
{
   delete[] _expressions;
}






/*!
 *
 * @param index  
 */
void ExpressionMatrix::ConstGene::read(int index)
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
bool ExpressionMatrix::ConstGene::readNext()
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
float ExpressionMatrix::ConstGene::at(int index) const
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
