#ifndef POWERLAW_H
#define POWERLAW_H
#include <ace/core/core.h>



class CorrelationMatrix;



class PowerLaw : public EAbstractAnalytic
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
   QVector<bool> computeAdjacencyMatrix(const QVector<float>& matrix, const QVector<float>& maximums, float threshold, int* size);
   QVector<int> computeDegreeDistribution(const QVector<bool>& matrix, int size);
   float computeCorrelation(const QVector<int>& histogram);

   CorrelationMatrix* _input {nullptr};
   QFile* _logfile {nullptr};
   float _thresholdStart {0.99};
   float _thresholdStep {0.01};
   float _thresholdStop {0.5};
};



#endif
