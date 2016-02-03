#ifndef LINUXFILE_H
#define LINUXFILE_H
#include <string>
#include <cstdint>
#include <vector>
#include "filemem.h"



/// @brief Linux OS FileMem implementation.
///
/// Implements the file memory object interface within the Linux OS environment.
/// Uses basic file IO calls in linux such as open/close/read/write and
/// posix_fallocate.
class LinuxFile : public FileMem
{
public:
   // *
   // * BASIC METHODS
   // *
   LinuxFile(const LinuxFile&) = delete;
   LinuxFile(LinuxFile&&) = delete;
   LinuxFile& operator=(const LinuxFile&) = delete;
   LinuxFile& operator=(LinuxFile&&) = delete;
   LinuxFile(const std::string&);
   ~LinuxFile();
   // *
   // * FUNCTIONS
   // *
   void clear() override final;
   bool reserve(int64_t) override final;
   uint64_t size() override final;
   uint64_t available() override final;
protected:
   // *
   // * FUNCTIONS
   // *
   void write(const void*,VPtr,uint64_t) override final;
   void read(void*,VPtr,uint64_t) override final;
   VPtr allocate(uint64_t) override final;
   VPtr head() override final;
private:
   // *
   // * VARIABLES
   // *
   /// File descriptor for file memory object.
   int _fd;
   /// Length of identification string in bytes.
   int _idLen;
   /// Total size of the file memory object in bytes.
   uint64_t _size;
   /// Space available for allocation in bytes.
   uint64_t _available;
   /// The next location in file memory that will be returned when new space is
   /// allocated.
   VPtr _next;
};



#endif
