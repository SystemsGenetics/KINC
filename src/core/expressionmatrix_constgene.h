#ifndef EXPRESSIONMATRIX_CONSTGENE_H
#define EXPRESSIONMATRIX_CONSTGENE_H
#include "expressionmatrix.h"
//



/*!
 */
class ExpressionMatrix::ConstGene
{
public:
   ConstGene(const ExpressionMatrix* matrix, bool isInitialized = false);
   ~ConstGene();
   void read(int index);
   bool readNext();
   float at(int index) const;
private:
   /*!
    */
   const ExpressionMatrix* _matrix;
   /*!
    */
   int _index {0};
   /*!
    */
   float* _expressions;
};



#endif
