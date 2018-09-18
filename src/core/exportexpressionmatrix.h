#ifndef EXPORTEXPRESSIONMATRIX_H
#define EXPORTEXPRESSIONMATRIX_H
#include <ace/core/core.h>

#include "expressionmatrix.h"



/*!
 * This class implements the export expression matrix analytic. This analytic
 * writes an expression matrix to a text file as table; that is, with each row
 * on a line, each value separated by whitespace, and the first row and column
 * containing the row names and column names, respectively. Elements which are
 * NAN in the expression matrix are written as the given NAN token.
 */
class ExportExpressionMatrix : public EAbstractAnalytic
{
   Q_OBJECT
public:
   class Input;
   virtual int size() const override final;
   virtual void process(const EAbstractAnalytic::Block* result) override final;
   virtual EAbstractAnalytic::Input* makeInput() override final;
   virtual void initialize();
private:
   /*!
    * Pointer to the input expression matrix.
    */
   ExpressionMatrix* _input {nullptr};
   /*!
    * Pointer to the output text file.
    */
   QFile* _output {nullptr};
   /*!
    * The string token used to represent NAN values.
    */
   QString _nanToken;
};



#endif
