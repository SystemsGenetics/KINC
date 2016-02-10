#ifndef FILEMEM_H
#define FILEMEM_H
#include "linuxfile.h"



template<class T> class RawFileMem : public T
{
public:
   using T::T;
   using T::Exception;
   using T::SystemError;
   using T::InvalidFile;
   using T::FileSegFault;
   using T::OutOfMemory;
   enum class Sync { read,write };
   using Ptr = typename T::Ptr;
   struct VPtr
   {
      Ptr ptr;
      RawFileMem<T>* mem;
   };
   using SizeT = typename T::SizeT;
   template<SizeT S> struct Object;
   struct DObject;
   template<class V> class Map;
   constexpr static Ptr nullPtr {T::nullPtr};
   void clear() { T::clear(); }
   bool reserve(size_t size) { return T::reserve(size); }
   size_t size() { return T::size(); }
   size_t capacity() { return T::capacity(); }
   VPtr head() { return {T::head(),this}; }
   VPtr allocate(size_t size) { return {T::allocate(size),this}; }
   void read(void* data,Ptr ptr,SizeT size) const
   { T::read(data,ptr,size); }
   void write(const void* data,Ptr ptr,SizeT size) { T::write(data,ptr,size); }
};
using FileMem = RawFileMem<LinuxFile>;



template<> template<FileMem::SizeT S> struct FileMem::Object
{
   template<class T> inline T& get(int n);
   char bytes[S];
   constexpr static SizeT size = S;
};



template<> struct FileMem::DObject
{
   inline DObject(const DObject&);
   inline DObject(DObject&&);
   inline DObject& operator=(const DObject&);
   inline DObject& operator=(DObject&&);
   inline DObject(SizeT);
   template<class T> inline T& get(int);
   inline ~DObject();
   char* bytes;
   SizeT size;
};



template<> template<class V> class FileMem::Map
{
public:
   inline void operator=(Ptr);//
   inline Map(VPtr);//
   inline Ptr addr() const;//
   inline void sync(Sync,uint64_t = 0);//
   inline V& operator*();
   inline V* operator->();
private:
   V _val;
   FileMem* _mem;
   Ptr _ptr;
};



template<> template<FileMem::SizeT S> template<class T>
inline T& FileMem::Object<S>::get(int n)
{
   return *reinterpret_cast<T*>(&bytes[n]);
}



inline FileMem::DObject::DObject(const DObject& obj):
   size(obj.size)
{
   bytes = new char[size];
   memcpy(bytes,obj.bytes,size);
}



inline FileMem::DObject::DObject(DObject&& tmp):
   size(tmp.size),
   bytes(tmp.bytes)
{
   tmp.bytes = nullptr;
   tmp.size = 0;
}



inline FileMem::DObject& FileMem::DObject::operator=(const DObject& obj)
{
   size = obj.size;
   if (bytes)
   {
      delete[] bytes;
   }
   bytes = new char[size];
   memcpy(bytes,obj.bytes,size);
}



inline FileMem::DObject& FileMem::DObject::operator=(DObject&& tmp)
{
   size = tmp.size;
   bytes = tmp.bytes;
   tmp.bytes = nullptr;
   tmp.size = 0;
}



inline FileMem::DObject::DObject(uint64_t dSize):
   size(dSize)
{
   bytes = new char[size];
}



inline FileMem::DObject::~DObject()
{
   if (bytes)
   {
      delete[] bytes;
   }
}



template<class T> inline T& FileMem::DObject::get(int n)
{
   return *reinterpret_cast<T*>(&bytes[n]);
}



template<> template<class V>
inline void FileMem::Map<V>::operator=(Ptr ptr)
{
   _ptr = ptr;
}



template<> template<class V>
inline FileMem::Map<V>::Map(VPtr vptr):
   _ptr(vptr.ptr),
   _mem(vptr.mem)
{}



template<> template<class V>
inline FileMem::Ptr FileMem::Map<V>::addr() const
{
   return _ptr;
}



template<> template<class V>
inline void FileMem::Map<V>::sync(Sync which, uint64_t inc)
{
   uint64_t offset = _ptr + inc*_val.size;
   switch (which)
   {
   case Sync::read:
      _mem->read(_val.bytes,offset,_val.size);
      break;
   case Sync::write:
      _mem->write(_val.bytes,offset,_val.size);
      break;
   }
}



template<> template<class V> inline V& FileMem::Map<V>::operator*()
{
   return _val;
}



template<> template<class V> inline V* FileMem::Map<V>::operator->()
{
   return &_val;
}



#endif
