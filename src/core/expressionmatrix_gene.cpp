#include "expressionmatrix_gene.h"
//






/*!
 *
 * @param matrix  
 *
 * @param read  
 */
ExpressionMatrix::Gene::Gene(ExpressionMatrix* matrix, bool read)
{}






/*!
 */
ExpressionMatrix::Gene::~Gene()
{}






/*!
 *
 * @param index  
 */
float& ExpressionMatrix::Gene::operator[](int index)
{}






/*!
 *
 * @param index  
 */
void ExpressionMatrix::Gene::read(int index)
{}






/*!
 */
void ExpressionMatrix::Gene::readNext()
{}






/*!
 *
 * @param index  
 */
void ExpressionMatrix::Gene::write(int index)
{
   // make sure given gene index is within range
   if ( index < 0 || index >= _matrix->_geneSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Attempting to write gene %1 when maximum is %2.").arg(index)
                   .arg(_matrix->_geneSize-1));
      throw e;
   }

   // write gene expressions to data object
   _matrix->writeGene(index,_expressions);
}






/*!
 */
void ExpressionMatrix::Gene::writeNext()
{}






/*!
 *
 * @param index  
 */
const float& ExpressionMatrix::Gene::at(int index) const
{}






/*! !!! UNKNOWN FUNCTION !!! */
void ExpressionMatrix::Gene::read(int index) const
{
   // make sure given gene index is within range
   if ( index < 0 || index >= _matrix->_geneSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Attempting to read gene %1 when maximum is %2.").arg(index)
                   .arg(_matrix->_geneSize-1));
      throw e;
   }

   // read gene expressions from data object
   _matrix->readGene(index,_expressions);
}






/*! !!! UNKNOWN FUNCTION !!! */
ExpressionMatrix::Expression& ExpressionMatrix::Gene::at(int index)
{
   // make sure given sample index is within range
   if ( index < 0 || index >= _matrix->_sampleSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Attempting to access gene expression %1 when maximum is %2.").arg(index)
                   .arg(_matrix->_sampleSize-1));
      throw e;
   }

   // return gene expression
   return _expressions[index];
}
