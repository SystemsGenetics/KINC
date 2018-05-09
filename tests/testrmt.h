#ifndef TESTRMT_H
#define TESTRMT_H
#include <QtTest/QtTest>

#include "pairwise_index.h"



class TestRMT : public QObject
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
