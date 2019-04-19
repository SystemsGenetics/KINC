#ifndef SIMILARITY_SERIAL_H
#define SIMILARITY_SERIAL_H
#include "similarity.h"
#include "pairwise_clusteringmodel.h"
#include "pairwise_correlationmodel.h"



/*!
 * This class implements the serial working class of the similarity analytic.
 */
class Similarity::Serial : public EAbstractAnalyticSerial
{
   Q_OBJECT
public:
   explicit Serial(Similarity* parent);
   virtual std::unique_ptr<EAbstractAnalyticBlock> execute(const EAbstractAnalyticBlock* block) override final;
private:
   int fetchPair(const Pairwise::Index& index, QVector<Pairwise::Vector2>& data, QVector<qint8>& labels);
   int removeOutliersCluster(const QVector<Pairwise::Vector2>& data, QVector<qint8>& labels, qint8 cluster, qint8 marker);
   int removeOutliers(const QVector<Pairwise::Vector2>& data, int numSamples, QVector<qint8>& labels, qint8 clusterSize, qint8 marker);
   /*!
    * Pointer to the base analytic for this object.
    */
   Similarity* _base;
   /*!
    * Pointer to the clustering model to use.
    */
   Pairwise::ClusteringModel* _clusModel {nullptr};
   /*!
    * Pointer to the correlation model to use.
    */
   Pairwise::CorrelationModel* _corrModel {nullptr};
};



#endif
