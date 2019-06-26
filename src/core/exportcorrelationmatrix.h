#ifndef EXPORTCORRELATIONMATRIX_H
#define EXPORTCORRELATIONMATRIX_H
#include <ace/core/core.h>

#include "ccmatrix_pair.h"
#include "ccmatrix.h"
#include "correlationmatrix_pair.h"
#include "correlationmatrix.h"
#include "expressionmatrix.h"



/*!
 * This class implements the export correlation matrix analytic. This analytic
 * takes two data objects, a correlation matrix and a cluster matrix, and writes
 * a text file of correlations, where each line is a correlation that includes
 * the pairwise index, correlation value, sample mask, and several other summary
 * statistics.
 */
class ExportCorrelationMatrix : public EAbstractAnalytic
{
   Q_OBJECT
public:
   class Input;
   virtual int size() const override final;
   virtual void process(const EAbstractAnalyticBlock* result) override final;
   virtual EAbstractAnalyticInput* makeInput() override final;
   virtual void initialize();
private:
   /**
    * Workspace variables to write to the output file
    */
   QTextStream _stream;
   CCMatrix::Pair _ccmPair;
   CorrelationMatrix::Pair _cmxPair;
   /*!
    * Pointer to the input expression matrix.
    */
   ExpressionMatrix* _emx {nullptr};
   /*!
    * Pointer to the input cluster matrix.
    */
   CCMatrix* _ccm {nullptr};
   /*!
    * Pointer to the input correlation matrix.
    */
   CorrelationMatrix* _cmx {nullptr};
   /*!
    * Pointer to the output text file.
    */
   QFile* _output {nullptr};
};



#endif
