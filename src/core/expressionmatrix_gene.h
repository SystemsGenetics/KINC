#ifndef EXPRESSIONMATRIX_GENE_H
#define EXPRESSIONMATRIX_GENE_H
#include "expressionmatrix.h"
//



/*!
 * This class implements the gene iterator for the expression matrix data
 * object. The gene iterator can read from or write to any gene (row) in the
 * expression matrix, or simply iterate through each row. The iterator stores
 * only one row of expression data in memory at a time.
 */
class ExpressionMatrix::Gene
{
public:
   float& operator[](int index);
public:
   Gene(ExpressionMatrix* matrix, bool isInitialized = false);
   ~Gene();
public:
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
    * The iterator's current position in the expression matrix.
    */
   int _index {0};
   /*!
    * Pointer to the expression data of the current gene.
    */
   float* _expressions;
};



#endif
