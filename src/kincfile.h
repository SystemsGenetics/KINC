#ifndef KINCFILE_H
#define KINCFILE_H
#include <string>
#include <memory>
#include "filemem.h"
#include "history.h"
#include "fstring.h"



namespace KincFileData
{
   struct Header;
   constexpr FileMem::SizeT idSz = 4;
   constexpr FileMem::SizeT hdrSz = idSz+24;
   constexpr auto idString = "\113\111\116\103";
}

struct KincFileData::Header : FileMem::Static<hdrSz>
{
   using FPtr = FileMem::Ptr;
   using FSizeT = FileMem::SizeT;
   using Static<hdrSz>::Static;
   char* idString() { &get<char>(0); }
   FPtr& histHead() { get<FPtr>(idSz); }
   FPtr& dataHead() { get<FPtr>(idSz+8); }
   const FPtr& dataHead() const { get<FPtr>(idSz+8); }
   FPtr& ident() { get<FPtr>(idSz+16); }
};



/// @ingroup dataplugin
/// @brief Base file utility class for data plugin.
///
/// Opens and manages a file memory object for a data plugin object. Provides
/// functions for the data plugin to interface with the file memory object. Also
/// provides a history object that is stored within the same file.
class KincFile
{
public:
   // *
   // * EXCEPTIONS
   // *
   struct Exception;
   struct InvalidFile;
   struct AlreadySet;
   // *
   // * DECLERATIONS
   // *
   using string = std::string;
   using FPtr = FileMem::Ptr;
   // *
   // * BASIC METHODS
   // *
   KincFile(const string&);
   // *
   // * COPY METHODS
   // *
   KincFile(const KincFile&) = delete;
   KincFile& operator=(const KincFile&) = delete;
   // *
   // * MOVE METHODS
   // *
   KincFile(KincFile&&) = delete;
   KincFile& operator=(KincFile&&) = delete;
   // *
   // * FUNCTIONS
   // *
   void clear();
   bool is_new();
   History& history();
protected:
   // *
   // * FUNCTIONS
   // *
   string ident() const;
   void ident(const string&);
   FPtr head() const;
   void head(FPtr);
   FileMem* mem();//NOT TESTED.
private:
   // *
   // * DECLERATIONS
   // *
   using Header = KincFileData::Header;
   using hptr = std::unique_ptr<History>;
   // *
   // * FUNCTIONS
   // *
   void create();
   // *
   // * CONSTANTS
   // *
   constexpr static auto _idString  = KincFileData::idString;
   constexpr static auto _idSz = KincFileData::idSz;
   constexpr static auto _hdrSz = KincFileData::hdrSz;
   // *
   // * VARIABLES
   // *
   /// File memory object.
   FileMem _mem;
   /// Remembers if the file opened is new or not.
   bool _new {true};
   /// Pointer to history object of file.
   hptr _hist {nullptr};
   /// Header for KINC file.
   Header _hdr;
   /// Custom data plugin ident value.
   FString _ident;
};



/// Tests to see if this object created a new KINC file when constructed.
///
/// @return True if this is a new file, else false.
inline bool KincFile::is_new()
{
   return _new;
}



/// Get reference of history object for this object.
///
/// @return History object.
inline History& KincFile::history()
{
   return *_hist;
}



/// Get data plugin ident value for this object.
///
/// @return data plugin ident.
inline KincFile::string KincFile::ident() const
{
   return *_ident;
}



/// @brief Get data plugin memory head.
///
/// Get file memory location for the beginning or head of data plugin memory for
/// this object.
///
/// @return Location to beginning of data plugin memory.
inline FileMem::Ptr KincFile::head() const
{
   return _hdr.dataHead();
}



/// Get pointer of file memory object for this object.
///
/// @return File memory object.
inline FileMem* KincFile::mem()
{
   return &_mem;
}



/// Generic base exception class for all exceptions thrown in KincFile class.
struct KincFile::Exception : public ::Exception
{
   using ::Exception::Exception;
};

/// The file being opened is not a valid KINC file.
struct KincFile::InvalidFile : public KincFile::Exception
{
   InvalidFile(const char* file, int line):
      Exception(file,line,"KincFile::InvalidFile")
   {}
};

/// The data plugin ident value has already been set.
struct KincFile::AlreadySet : public KincFile::Exception
{
   AlreadySet(const char* file, int line):
      Exception(file,line,"KincFile::AlreadySet")
   {}
};



#endif
