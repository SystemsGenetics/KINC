#ifndef TESTRMT_H
#define TESTRMT_H
#include <QtTest/QtTest>

#include "genepair_index.h"



class TestRMT : public QObject
{
	Q_OBJECT

private:
	struct Pair
	{
		GenePair::Index index;
		QVector<QVector<qint8>> clusters;
		QVector<float> correlations;
	};

private slots:
	void test();
};



#endif
