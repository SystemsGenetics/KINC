#ifndef EXTRACT_H
#define EXTRACT_H
#include <ace/core/core.h>

#include "ccmatrix.h"
#include "correlationmatrix.h"
#include "expressionmatrix.h"



/*!
 * This class implements the extract analytic. This analytic is very similar to
 * the export correlation matrix analytic, except for a few differences: (1) this
 * analytic uses a slightly different format for the text file, (2) this analytic
 * can apply a correlation threshold, and (3) this analytic can optionally write
 * a GraphML file. The key difference is that this analytic "extracts" a network
 * from the correlation matrix and writes an edge list rather than a correlation
 * list.
 */
class Extract : public EAbstractAnalytic
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
   /*!
    * Pointer to the output GraphML file.
    */
   QFile* _graphml {nullptr};
   /*!
    * The minimum (absolute) correlation threshold.
    */
   float _minCorrelation {0.85};
   /*!
    * The maximum (absolute) correlation threshold.
    */
   float _maxCorrelation {1.00};
};



#endif
