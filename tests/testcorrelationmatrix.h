#ifndef TESTCORRELATIONMATRIX_H
#define TESTCORRELATIONMATRIX_H
#include <QtTest/QtTest>

#include "genepair_index.h"



class TestCorrelationMatrix : public QObject
{
	Q_OBJECT

private:
	struct Pair
	{
		GenePair::Index index;
		QVector<float> clusters;
	};

private slots:
	void test();
};



#endif
