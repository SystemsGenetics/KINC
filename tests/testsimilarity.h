#ifndef TESTSIMILARITY_H
#define TESTSIMILARITY_H
#include <QtTest/QtTest>

#include "genepair_index.h"



class TestSimilarity : public QObject
{
	Q_OBJECT

private:
	struct Pair
	{
		Pairwise::Index index;
		QVector<QVector<qint8>> clusters;
		QVector<float> correlations;
	};

private slots:
	void test();
};



#endif
