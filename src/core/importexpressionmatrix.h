#ifndef IMPORTEXPRESSIONMATRIX_H
#define IMPORTEXPRESSIONMATRIX_H
#include <ace/core/core.h>

#include "expressionmatrix.h"



/*!
 * This class implements the import expression matrix analytic. This analytic
 * reads in a text file which contains a matrix as a table; that is, with each row
 * on a line, each value separated by whitespace, and the first row and column
 * containing the row names and column names, respectively. Elements which have
 * the given NAN token are read in as NAN. Additionally, the user can specify a
 * fixed number of samples (columns) to import; by default the analytic imports
 * the entire matrix.
 */
class ImportExpressionMatrix : public EAbstractAnalytic
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
    * Pointer to the input text file.
    */
   QFile* _input {nullptr};
   /*!
    * Pointer to the output expression matrix.
    */
   ExpressionMatrix* _output {nullptr};
   /*!
    * The string token used to represent NAN values.
    */
   QString _nanToken {"NA"};
   /*!
    * The number of samples to read.
    */
   qint32 _sampleSize {0};
};



#endif
