#ifndef HISTORY_H
#define HISTORY_H
#include <string>
#include "filemem.h"
#include "exception.h"



namespace HistoryData
{
   struct Header;
   struct Node;
   struct String;
   constexpr static int _idLen = 9;
   constexpr static int _nodeLen = 60;
   constexpr static int _hdrLen = _idLen+16;
}

struct HistoryData::Header : FileMem::Static<_hdrLen>
{
   using Static<_hdrLen>::Static;
   FileMem::Ptr& head() { get<FileMem::Ptr>(0); }
   FileMem::Ptr& dataHead() { get<FileMem::Ptr>(8); }
   char* ident() { &get<char>(16); }
};

struct HistoryData::Node : FileMem::Static<_nodeLen>
{
   using Static<_nodeLen>::Static;
   int64_t& timeStamp() { get<int64_t>(0); }
   FileMem::Ptr& fileNamePtr() { get<FileMem::Ptr>(8); }
   int32_t& fileNameSize() { get<int32_t>(16); }
   FileMem::Ptr& objectPtr() { get<FileMem::Ptr>(20); }
   int32_t& objectSize() { get<int32_t>(28); }
   FileMem::Ptr& commandPtr() { get<FileMem::Ptr>(32); }
   int32_t& commandSize() { get<int32_t>(40); }
   FileMem::Ptr& childHead() { get<FileMem::Ptr>(44); }
   FileMem::Ptr& next() { get<FileMem::Ptr>(52); }
};

struct HistoryData::String : FileMem::Object
{
   using Object::Object;
   String(FileMem::SizeT size): Object(size) {}
   char* c_str() { &get<char>(0); }
};



class History
{
public:
   // *
   // * EXCEPTIONS
   // *
   struct Exception;
   struct InvalidFile;
   struct AlreadySet;
   struct NotSet;
   // *
   // * DECLERATIONS
   // *
   class Iterator;
   // *
   // * BASIC METHODS
   // *
   History(const History&) = delete;
   History(History&&) = delete;
   History& operator=(const History&) = delete;
   History& operator=(History&&) = delete;
   History(const std::string&);
   // *
   // * FUNCTIONS
   // *
   void set_timestamp(int64_t);
   void set_filename(const std::string&);
   void set_object(const std::string&);
   void set_command(const std::string&);
   void add_child(History&);
   Iterator begin();
   Iterator end();
private:
   // *
   // * DECLERATIONS
   // *
   using Header = HistoryData::Header;
   using Node = HistoryData::Node;
   using String = HistoryData::String;
   // *
   // * FUNCTIONS
   // *
   FileMem::Ptr write_string(const std::string&);
   FileMem::Ptr rec_add_child(Iterator);
   void copy_child(Node&,Iterator&);
   // *
   // * CONSTANTS
   // *
   constexpr static int _idLen = HistoryData::_idLen;
   constexpr static int _nodeLen = HistoryData::_nodeLen;
   constexpr static int _hdrLen = HistoryData::_hdrLen;
   constexpr static int _emptyLen = _hdrLen + _nodeLen;
   constexpr const static char* _identString = "\113\111\116\103\102\106\111"
                                               "\114\105";
   // *
   // * VARIABLES
   // *
   FileMem _mem;
   Header _header;
   Node _info;
};



class History::Iterator
{
public:
   // *
   // * DECLERATIONS
   // *
   friend class History;
   // *
   // * FUNCTIONS
   // *
   int64_t timestamp();
   std::string filename();
   std::string object();
   std::string command();
   Iterator childHead();
   // *
   // * OPERATORS
   // *
   Iterator& operator++();
   bool operator!=(const Iterator&); ///TODO
private:
   // *
   // * BASIC METHODS
   // *
   Iterator(FileMem::Ptr,FileMem&);
   // *
   // * OPERATORS
   // *
   Node& operator*();
   // *
   // * VARIABLES
   // *
   FileMem& _mem;
   FileMem::Ptr _ptr;
   Node _info;
};



struct History::Exception : public ::Exception
{
   using ::Exception::Exception;
};

struct History::InvalidFile : public History::Exception
{
   InvalidFile(const char* file, int line):
      Exception(file,line,"History::InvalidFile")
   {}
};

struct History::AlreadySet : public History::Exception
{
   AlreadySet(const char* file, int line):
      Exception(file,line,"History::AlreadySet")
   {}
};

struct History::NotSet : public History::Exception
{
   NotSet(const char* file, int line):
      Exception(file,line,"History::NotSet")
   {}
};

#endif
