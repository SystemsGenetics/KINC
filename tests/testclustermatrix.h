#ifndef TESTClUSTERMATRIX_H
#define TESTCLUSTERMATRIX_H
#include <QtTest/QtTest>

#include "genepair_index.h"



class TestClusterMatrix : public QObject
{
	Q_OBJECT

private:
	struct Pair
	{
		Pairwise::Index index;
		QVector<QVector<qint8>> clusters;
	};

private slots:
	void test();
};



#endif
