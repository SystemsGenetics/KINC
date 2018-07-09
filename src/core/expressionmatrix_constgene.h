#ifndef EXPRESSIONMATRIX_CONSTGENE_H
#define EXPRESSIONMATRIX_CONSTGENE_H
//



/*!
 */
class ExpressionMatrix::ConstGene
{
public:
   ConstGene(const ExpressionMatrix* matrix, bool read = false);
   ~ConstGene();
   void read(int index);
   void readNext();
   const float& at() const;
private:
   /*!
    */
   const ExpressionMatrix* _matrix;
   /*!
    */
   int _index {0};
   /*!
    */
   int _expressions;
};



#endif
