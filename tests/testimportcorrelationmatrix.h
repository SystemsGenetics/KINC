#ifndef TESTIMPORTCORRELATIONMATRIX_H
#define TESTIMPORTCORRELATIONMATRIX_H
#include <QtTest/QtTest>

#include "genepair_index.h"



class TestImportCorrelationMatrix : public QObject
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
