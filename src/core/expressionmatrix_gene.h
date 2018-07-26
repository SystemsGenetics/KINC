#ifndef EXPRESSIONMATRIX_GENE_H
#define EXPRESSIONMATRIX_GENE_H
#include "expressionmatrix.h"
//



/*!
 */
class ExpressionMatrix::Gene
{
public:
   float& operator[](int index);
public:
   Gene(ExpressionMatrix* matrix, bool isInitialized = false);
   ~Gene();
   void read(int index);
   bool readNext();
   void write(int index);
   bool writeNext();
   float at(int index) const;
private:
   /*!
    */
   ExpressionMatrix* _matrix;
   /*!
    */
   int _index {0};
   /*!
    */
   float* _expressions;
};



#endif
