#ifndef AMATRIX_H
#define AMATRIX_H
#include <ace.h>



namespace Ace = AccelCompEng;



class AMatrix : public Ace::Data, private Ace::NVMemory::Node
{
public:
   struct AlreadySet : public Ace::Exception { using Ace::Exception::Exception; };
   struct NullMatrix : public Ace::Exception { using Ace::Exception::Exception; };
   struct OutOfRange : public Ace::Exception { using Ace::Exception::Exception; };
   struct InvalidSize : public Ace::Exception { using Ace::Exception::Exception; };
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
      int64_t _genePtr {fnullptr};
      int64_t _netData {fnullptr};
   };
   Header& data() { return get<Header>(); }
   const Header& data() const { return get<Header>(); }
   void null_data() override final {}
   void flip_endian() override final;
   bool _isNew {true};
   std::vector<std::string> _geneNames;
   std::unique_ptr<Iterator> _iterator {nullptr};
   constexpr static int _strSize {1024};
};



class AMatrix::Iterator : private Ace::NVMemory::Node
{
public:
   friend class AMatrix;
   void read();
   void write();
   int x() const;
   int y() const;
   uint8_t& operator*();
   bool operator!=(const Iterator&);
   void operator++();
private:
   Iterator(AMatrix* p, int x, int y);
   void set(int x, int y);
   static int64_t initialize(Ace::NVMemory& mem, int geneSize);
   void null_data() override final {}
   void flip_endian() override final {}
   AMatrix* _p;
   int _x;
   int _y;
};



#endif
