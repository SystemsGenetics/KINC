#ifndef AMATRIX_H
#define AMATRIX_H
#include <ace.h>



namespace Ace = AccelCompEng;



class AMatrix : public Ace::Data, private Ace::NVMemory::Node
{
public:
   class Iterator;
   AMatrix();
   void init() override final;
   void load(Ace::GetOpts &ops, Ace::Terminal &tm) override final;
   void dump(Ace::GetOpts &ops, Ace::Terminal &tm) override final;
   void query(Ace::GetOpts &ops, Ace::Terminal &tm) override final;
   bool empty() override final;
   void initialize(std::vector<std::string>&& geneNames);
   Iterator begin();
   Iterator end();
   Iterator& at(int x, int y);
   Iterator& ref(int x, int y);
private:
   struct __attribute__ ((__packed__)) Header
   {
      int32_t _geneSize {0};
      int64_t _netData {fnullptr};
   };
   void null_data() override final {}
   void flip_endian() override final;
   bool _isNew {true};
   std::unique_ptr<Iterator> _iterator {nullptr};
};



class AMatrix::Iterator : private Ace::NVMemory::Node
{
public:
   void read();
   void write();
   bool operator*();
   bool operator!=(const Iterator&);
   void operator++();
private:
   Iterator(AMatrix* p, int i);
   void null_data() override final {}
   void flip_endian() override final {}
   AMatrix* _p;
   int _i;
};



#endif
