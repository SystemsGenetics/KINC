#ifndef TESTCORRELATIONMATRIX_H
#define TESTCORRELATIONMATRIX_H
#include <QtTest/QtTest>

#include "pairwise_index.h"



class TestCorrelationMatrix : public QObject
{
	Q_OBJECT

private:
	struct Pair
	{
		Pairwise::Index index;
		QVector<float> correlations;
	};

private slots:
	void test();
};



#endif
