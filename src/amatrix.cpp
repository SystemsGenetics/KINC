#include "amatrix.h"



AMatrix::AMatrix() {}
void AMatrix::init() {}
void AMatrix::load(Ace::GetOpts &ops, Ace::Terminal &tm) {}
void AMatrix::dump(Ace::GetOpts &ops, Ace::Terminal &tm) {}
void AMatrix::query(Ace::GetOpts &ops, Ace::Terminal &tm) {}
bool AMatrix::empty() {}
void AMatrix::initialize(std::vector<std::string>&& geneNames) {}
AMatrix::Iterator AMatrix::begin() {}
AMatrix::Iterator AMatrix::end() {}
AMatrix::Iterator& AMatrix::at(int x, int y) {}
AMatrix::Iterator& AMatrix::ref(int x, int y) {}
void AMatrix::flip_endian() {}
