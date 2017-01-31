#ifndef AMATRIX_H
#define AMATRIX_H
#include <ace.h>



namespace Ace = AccelCompEng;



/// Data plugin class that holds adjacency data for gene statistics.
class AMatrix : public Ace::Data, private Ace::NVMemory::Node
{
public:
   struct AlreadySet : public Ace::Exception { using Ace::Exception::Exception; };
   struct NullMatrix : public Ace::Exception { using Ace::Exception::Exception; };
   struct OutOfRange : public Ace::Exception { using Ace::Exception::Exception; };
   struct InvalidSize : public Ace::Exception { using Ace::Exception::Exception; };
   class Iterator;
   AMatrix();
   /// ACE initialize call.
   void init() override final;
   /// ACE load call.
   void load(Ace::GetOpts &ops, Ace::Terminal &tm) override final;
   /// ACE dump call.
   void dump(Ace::GetOpts &ops, Ace::Terminal &tm) override final;
   /// ACE query call.
   void query(Ace::GetOpts &ops, Ace::Terminal &tm) override final;
   /// ACE empty call.
   bool empty() override final;
   /// Initialize a new adjacency matrix for this object.
   /// @param geneNames List of all gene names.
   void initialize(std::vector<std::string>&& geneNames);
   /// First gene edge iterator.
   Iterator begin();
   /// End of list gene edge iterator.
   Iterator end();
   /// Get gene edge.
   /// @param x Index of first gene.
   /// @param y Index of second gene.
   /// @return Gene edge iterator.
   Iterator& at(int x, int y);
   /// Get reference to gene edge.
   /// @param x Index of first gene.
   /// @param y Index of second gene.
   /// @return Gene edge reference.
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
