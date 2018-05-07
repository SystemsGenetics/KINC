#ifndef TESTEXPORTCORRELATIONMATRIX_H
#define TESTEXPORTCORRELATIONMATRIX_H
#include <QtTest/QtTest>

#include "genepair_index.h"



class TestExportCorrelationMatrix : public QObject
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
