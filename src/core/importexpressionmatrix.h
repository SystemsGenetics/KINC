#ifndef IMPORTEXPRESSIONMATRIX_H
#define IMPORTEXPRESSIONMATRIX_H
#include <ace/core/core.h>

#include "expressionmatrix.h"



/*!
 * This class implements the import expression matrix analytic. This analytic
 * reads in a text file which contains a matrix as a table; that is, with each row
 * on a line, each value separated by whitespace, and the first row and column
 * containing the row names and column names, respectively. Elements which have
 * the given NAN token are read in as NAN. If the sample names are not in the
 * input file, the user must provide the number of samples to the analytic, and
 * the samples will be given integer names.
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
   /**
    * Structure used to load gene expression data
    */
   struct Gene
   {
      Gene() = default;
      Gene(int size)
      {
         expressions.resize(size);
      }

      QVector<float> expressions;
   };

   /**
    * Workspace variables to hold gene expression data
    */
   QTextStream _stream;
   int _numLines {0};
   QVector<Gene> _genes;
   QStringList _geneNames;
   QStringList _sampleNames;

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
