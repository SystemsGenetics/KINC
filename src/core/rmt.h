#ifndef RMT_H
#define RMT_H
#include <ace/core/core.h>



class CorrelationMatrix;



class RMT : public EAbstractAnalytic
{
   Q_OBJECT
public:
   class Input;
   virtual int size() const override final;
   virtual void process(const EAbstractAnalytic::Block* result) override final;
   virtual EAbstractAnalytic::Input* makeInput() override final;
   virtual void initialize();
private:
   QVector<float> computeMaximums(const QVector<float>& matrix);
   QVector<float> computePruneMatrix(const QVector<float>& matrix, const QVector<float>& maximums, float threshold, int* size);
   QVector<float> computeEigenvalues(QVector<float>* pruneMatrix, int size);
   float computeChiSquare(const QVector<float>& eigens);
   float computePaceChiSquare(const QVector<float>& eigens, int pace);
   QVector<float> degenerate(const QVector<float>& eigens);
   QVector<float> unfold(const QVector<float>& eigens, int pace);

   CorrelationMatrix* _input {nullptr};
   QFile* _logfile {nullptr};
   float _thresholdStart {0.99};
   float _thresholdStep {0.001};
   float _thresholdStop {0.5};
   float _chiSquareThreshold1 {99.607};
   float _chiSquareThreshold2 {200};
   int _minEigenvalueSize {50};
   int _minUnfoldingPace {10};
   int _maxUnfoldingPace {40};
   int _histogramBinSize {60};
};



#endif
