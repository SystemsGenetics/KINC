#ifndef FILEMEM_H
#define FILEMEM_H
#include <string>
#include <cstdint>
#include <type_traits>
#include <utility>
#include "exception.h"



/// @brief File memory interface class.
///
/// Defines the behavior of a file memory class. An object implementing this
/// interface takes a file and maps it as a memory object; providing pointers
/// to portions of the file that the user can read and write from using the
/// FileMem::VPtr class defined elsewhere.
///
/// @author Josh Burns
/// @date 2 Feb 2016
class FileMem
{
public:
   // *
   // * EXCEPTIONS
   // *
   struct Exception;
   struct SystemError;
   struct InvalidFile;
   struct FileSegFault;
   struct OutOfMemory;
   // *
   // * DECLERATIONS
   // *
   using VPtr = uint64_t;
   template<class T> class Ptr;
   // *
   // * BASIC METHODS
   // *
   virtual ~FileMem() = default;
   // *
   // * VIRTUAL FUNCTIONS
   // *
   /// Must clear any allocated memory on the object's file.
   virtual void clear() = 0;
   /// @brief Add additional space to file.
   ///
   /// Must reserve additional unallocated space on the file for memory
   /// allocation.
   ///
   /// @param newBytes The number of new bytes to reserve from object's file.
   /// @return True if the reservation of additional space was successful.
   ///
   /// @post The allocation of additional space from the reserved amount will
   /// not fail.
   virtual bool reserve(int64_t newBytes) = 0;
   /// Must return the total size in bytes of the file memory object.
   virtual uint64_t size() const = 0;
   /// @brief Return available space.
   ///
   /// Must return the total amount of bytes available for memory
   /// allocation.
   virtual uint64_t available() const = 0;
   // *
   // * CONSTANTS
   // *
   constexpr static VPtr nullPtr = 0xffffffffffffffffll;
   /// Must be at the beginning of any file used as a memory object of this
   /// interface.
   constexpr const static char* _identString = "KINCBINDAT";
protected:
   // *
   // * VIRTUAL FUNCTIONS
   // *
   /// Must write a block of binary data to the object's file.
   ///
   /// @param data Binary data in memory that will be written to the file.
   /// @param ptr Specifies location in file to write.
   /// @param size Total number of bytes that will be written.
   virtual void write(const void* data, VPtr ptr, uint64_t size) = 0;
   /// Must read a block of binary data from the object's file.
   ///
   /// @param data Binary data in memory that will be written to from the file.
   /// @param ptr Specifies location in file to read.
   /// @param size Total number of bytes that will be read.
   virtual void read(void*,VPtr,uint64_t) const = 0;
   /// @brief Allocate new file memory space.
   ///
   /// Must allocate new space from the object file and return a new file memory
   /// pointer that points to the beginning of that newly allocated space.
   ///
   /// @param size Total size in bytes of space to be allocated.
   /// @return Points to newly allocated space on file object.
   virtual VPtr allocate(uint64_t size) = 0;
   /// @brief Provide head pointer.
   ///
   /// Must return a file memory pointer that points to the very beginning of
   /// file memory space.
   ///
   /// @return Points to beginning of file memory space.
   virtual VPtr head() const = 0;
};



template<class T> class FileMem::Ptr
{
public:
   static_assert(std::is_class<T>::value,
                 "FileMem::Ptr<T>, T is of an unsupported type.");
   Ptr(const Ptr<T>&) = delete;
   Ptr<T>& operator=(const Ptr<T>&) = delete;
   Ptr(FileMem&,T*);
   Ptr(FileMem&,const T&);
   Ptr(FileMem&,T&&);
   Ptr(FileMem&,VPtr);
   Ptr(Ptr<T>&&);
   ~Ptr();
   inline Ptr<T>& operator=(Ptr<T>&&);
   inline Ptr<T>& operator=(VPtr);
   inline T& operator*();
   inline const T& operator*() const;
   inline T* operator->();
   inline const T* operator->() const;
   inline VPtr addr() const;
   inline void save();
private:
   T* _val;
   FileMem& _base;
   VPtr _ptr;
};



template<class T> FileMem::Ptr<T>::Ptr(FileMem& mem, T* val):
   _val(val),
   _base(mem),
   _ptr(nullPtr)
{
   _ptr = _base.allocate(val.bin_size());
}



template<class T> FileMem::Ptr<T>::Ptr(FileMem& mem, const T& val):
   _val(nullptr),
   _base(mem),
   _ptr(nullPtr)
{
   _ptr = _base.allocate(val.bin_size());
   _val = new T {val};
}



template<class T> FileMem::Ptr<T>::Ptr(FileMem& mem, T&& val):
   _val(nullptr),
   _base(mem),
   _ptr(nullPtr)
{
   _val = new T {std::move(val)};
   _ptr = _base.allocate(val.bin_size());
}



template<class T> FileMem::Ptr<T>::Ptr(FileMem& mem, VPtr loc):
   _val(nullptr),
   _base(mem),
   _ptr(loc)
{
   _val = new T;
   _base.read(_val->bin_data(),_ptr,_val->bin_size());
}



template<class T> FileMem::Ptr<T>::Ptr(FileMem::Ptr<T>&& tmp):
   _val(tmp._val),
   _base(tmp._base),
   _ptr(tmp._ptr)
{
   tmp._val = nullptr;
   tmp._ptr = nullPtr;
}



template<class T> FileMem::Ptr<T>::~Ptr()
{
   if (_val)
   {
      delete _val;
   }
}



template<class T> inline FileMem::Ptr<T>&
FileMem::Ptr<T>::operator=(Ptr<T>&& tmp)
{
   _val = tmp._val;
   _base = tmp._base;
   _ptr = tmp._ptr;
   tmp._val = nullptr;
   tmp._ptr = nullPtr;
}



template<class T> inline FileMem::Ptr<T>& FileMem::Ptr<T>::operator=(VPtr ptr)
{
   _ptr = ptr;
   _base.read(_val->bin_data(),_ptr,_val->bin_size());
}



template<class T> inline T& FileMem::Ptr<T>::operator*()
{
   return *_val;
}



template<class T> inline const T& FileMem::Ptr<T>::operator*() const
{
   return *_val;
}



template<class T> inline T* FileMem::Ptr<T>::operator->()
{
   return _val;
}



template<class T> inline const T* FileMem::Ptr<T>::operator->() const
{
   return _val;
}



template<class T> inline FileMem::VPtr FileMem::Ptr<T>::addr() const
{
   return _ptr;
}



template<class T> inline void FileMem::Ptr<T>::save()
{
   _base.write(_val->bin_data(),_ptr,_val->bin_size());
}



/// Generic base exception class for all exceptions thrown in FileMem class.
struct FileMem::Exception : public ::Exception
{
   using ::Exception::Exception;
};

/// Exception that is thrown when a system error occurs.
struct FileMem::SystemError : public ::SystemError
{
   using ::SystemError::SystemError;
};

/// Exception that is thrown when a file that is not a valid file memory object
/// is opened.
struct FileMem::InvalidFile : public FileMem::Exception
{
   InvalidFile(const char* file, int line):
      Exception(file,line,"FileMem::InvalidFile")
   {}
};

/// Exception that is thrown when a read or write operation is attempted that
/// goes beyond the limits of the total file object's size.
struct FileMem::FileSegFault : public FileMem::Exception
{
   FileSegFault(const char* file, int line):
      Exception(file,line,"FileMem::FileSegFault")
   {}
};

/// Exception that is thrown when an allocation of new file memory is attempted
/// that is greater than the available amount of bytes that can be allocated.
struct FileMem::OutOfMemory : public FileMem::Exception
{
   OutOfMemory(const char* file, int line):
      Exception(file,line,"FileMem::OutOfMemory")
   {}
};



#endif
