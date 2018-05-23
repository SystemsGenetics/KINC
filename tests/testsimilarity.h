#ifndef TESTSIMILARITY_H
#define TESTSIMILARITY_H
#include <QtTest/QtTest>

#include "pairwise_index.h"



class TestSimilarity : public QObject
{
	Q_OBJECT

private:
	struct Pair
	{
		Pairwise::Index index;
		QVector<QVector<qint8>> sampleMasks;
		QVector<float> correlations;
	};

private slots:
	void test();
};



#endif
