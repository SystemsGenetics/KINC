#ifndef EXPRESSIONMATRIX_GENE_H
#define EXPRESSIONMATRIX_GENE_H
#endif // EXPRESSIONMATRIX_GENE_H
//



/*!
 */
class ExpressionMatrix::Gene
{
public:
   Gene(ExpressionMatrix* matrix, bool read = false);
   ~Gene();
   float& operator[](int index);
   void read(int index);
   void readNext();
   void write(int index);
   void writeNext();
   const float& at(int index) const;
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
