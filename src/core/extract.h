#ifndef EXTRACT_H
#define EXTRACT_H
#include <ace/core/core.h>

#include "ccmatrix_pair.h"
#include "ccmatrix.h"
#include "correlationmatrix_pair.h"
#include "correlationmatrix.h"
#include "expressionmatrix.h"
#include "condition-specificclustersmatrix.h"


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
public:
   virtual int size() const override final;
   virtual void process(const EAbstractAnalyticBlock* result) override final;
   virtual EAbstractAnalyticInput* makeInput() override final;
   virtual void initialize();
   QVector<QString> formatAnnotations();
private:
   /*!
   * Defines the output formats this analytic supports.
   */
   enum class OutputFormat
   {
      /*!
       * Text format
       */
      Text
      /*!
       * Minimal format (does not use CCM)
       */
      ,Minimal
      /*!
       * GraphML format
       */
      ,GraphML
   };
private:
   void writeTextFormat(int index);
   void writeMinimalFormat(int index);
   void writeGraphMLFormat(int index);
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
    * Pointer to the input condition specific cluster matrix.
    */
   CSCM* _cscm {nullptr};
   /*!
    * Pointer to the annotation file.
    */
   QFile* _amx {nullptr};
   /*!
    * The output format to use.
    */
   OutputFormat _outputFormat {OutputFormat::Text};
   /*!
    * Pointer to the output text file.
    */
   QFile* _output {nullptr};
   /*!
    * The minimum (absolute) correlation threshold.
    */
   float _minCorrelation {0.85f};
   /*!
    * The maximum (absolute) correlation threshold.
    */
   float _maxCorrelation {1.00f};
};



#endif
