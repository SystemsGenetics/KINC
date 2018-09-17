#ifndef EXPRESSIONMATRIX_GENE_H
#define EXPRESSIONMATRIX_GENE_H
#include "expressionmatrix.h"
//



/*!
 * This class implements the gene iterator for the expression matrix data
 * object. The gene iterator can iterate through each gene in the expression
 * matrix and read or write to each gene. The iterator stores only one row of
 * expression data in memory at a time.
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
    * Pointer to the parent expression matrix.
    */
   ExpressionMatrix* _matrix;
   /*!
    * The iterator's current location in the expression matrix.
    */
   int _index {0};
   /*!
    * Pointer to the expression data of the current gene.
    */
   float* _expressions;
};



#endif
