#ifndef LINUXFILE_H
#define LINUXFILE_H
#include <string>
#include <cstdint>
#include <vector>
#include "filemem.h"



// posix_fallocate(int fd, off_t offset, off_t len);
class LinuxFile : public FileMem
{
public:
   LinuxFile(const LinuxFile&) = delete;
   LinuxFile(LinuxFile&&) = delete;
   LinuxFile& operator=(const LinuxFile&) = delete;
   LinuxFile& operator=(LinuxFile&&) = delete;
   LinuxFile(const std::string&);
   ~LinuxFile();
   void clear() override final;
   bool reserve(int64_t) override final;
   int64_t size() override final;
   int64_t available() override final;
protected:
   void write(void*,VPtr,int) override final;
   void read(void*,VPtr,int) override final;
   VPtr allocate(int) override final;
private:
   int _fd;
   int _idLen;
   int64_t _size;
   int64_t _available;
   VPtr _next;
};



#endif
