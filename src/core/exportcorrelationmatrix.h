#ifndef EXPORTCORRELATIONMATRIX_H
#define EXPORTCORRELATIONMATRIX_H
#include <ace/core/core.h>

#include "ccmatrix.h"
#include "correlationmatrix.h"
#include "expressionmatrix.h"



/*!
 * This class implements the export correlation matrix analytic. This analytic
 * takes two data objects, a correlation matrix and a cluster matrix, and writes
 * a text file of correlations, where each line is a correlation that includes
 * the pairwise index, correlation value, and sample mask, as well as several
 * other fields which are not used but are required for this format. The analytic
 * attempts to recreate these fields as much as is possible. The expression matrix
 * that was used to produce the correlation matrix must also be provided in order
 * to recreate sample masks for pairs with only one cluster, as these sample masks
 * are not stored in the cluster matrix.
 */
class ExportCorrelationMatrix : public EAbstractAnalytic
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
