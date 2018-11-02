#ifndef SIMILARITY_INPUT_H
#define SIMILARITY_INPUT_H
#include "similarity.h"



/*!
 * This class implements the abstract input of the similarity analytic.
 */
class Similarity::Input : public EAbstractAnalytic::Input
{
   Q_OBJECT
public:
   /*!
    * Defines all arguments for its parent analytic.
    */
   enum Argument
   {
      InputData = 0
      ,ClusterData
      ,CorrelationData
      ,ClusteringType
      ,CorrelationType
      ,MinExpression
      ,MinSamples
      ,MinClusters
      ,MaxClusters
      ,CriterionType
      ,RemovePreOutliers
      ,RemovePostOutliers
      ,MinCorrelation
      ,MaxCorrelation
      ,WorkBlockSize
      ,GlobalWorkSize
      ,LocalWorkSize
      ,Total
   };
   explicit Input(Similarity* parent);
   virtual int size() const override final;
   virtual EAbstractAnalytic::Input::Type type(int index) const override final;
   virtual QVariant data(int index, Role role) const override final;
   virtual void set(int index, const QVariant& value) override final;
   virtual void set(int index, QFile* file) override final;
   virtual void set(int index, EAbstractData* data) override final;
private:
   static const QStringList CLUSTERING_NAMES;
   static const QStringList CORRELATION_NAMES;
   static const QStringList CRITERION_NAMES;
   /*!
    * Pointer to the base analytic for this object.
    */
   Similarity* _base;
};



#endif
