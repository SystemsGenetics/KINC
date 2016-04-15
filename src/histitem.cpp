#include "histitem.h"



/// @brief Initializes history item.
///
/// Initializes history item object as being empty. If file location given is
/// not nullptr then loads item from file memory location.
///
/// @param mem File memory object.
/// @param ptr Location in file memory of item, nullptr if none.
HistItem::HistItem(FileMem* mem, FileMem::Ptr ptr):
   _mem(mem),
   _item(ptr),
   _fileName(mem),
   _object(mem),
   _command(mem)
{
   if (_item.addr()!=FileMem::nullPtr)
   {
      load_item();
   }
}



/// Moves given history object.
///
/// @param tmp Object to take data from.
HistItem::HistItem(HistItem&& tmp):
   _mem(tmp._mem),
   _item(tmp._item),
   _fileName(std::move(tmp._fileName)),
   _object(std::move(tmp._object)),
   _command(std::move(tmp._command))
{
   tmp._item = FileMem::nullPtr;
}



/// Moves given history object, overwriting any data this object holds.
///
/// @param tmp Object to take data from.
HistItem& HistItem::operator=(HistItem&& tmp)
{
   _mem = tmp._mem;
   _item = tmp._item;
   _fileName = std::move(tmp._fileName);
   _object = std::move(tmp._object);
   _command = std::move(tmp._command);
   tmp._item = FileMem::nullPtr;
}



/// Allocate a new history item to file memory.
///
/// @exception IsAllocated This history object has already been allocated.
void HistItem::allocate()
{
   bool cond = _item.addr()==FileMem::nullPtr;
   assert<IsAllocated>(cond,__FILE__,__LINE__);
   _mem->allot(_item);
   _item.timeStamp() = 0;
   _item.fileNamePtr() = FileMem::nullPtr;
   _item.objectPtr() = FileMem::nullPtr;
   _item.commandPtr() = FileMem::nullPtr;
   _item.childHead() = FileMem::nullPtr;
   _item.next() = FileMem::nullPtr;
}



/// @brief Make recursive copy of history.
///
/// Make a copy of the given history item, including new copies of all of given
/// history item's children.
///
/// @param hist History item to recursively copy.
///
/// @exception IsAllocated This history object has already been allocated.
/// @exception IsNullPtr The history item given for copying is not set(its
/// address is nullptr).
void HistItem::copy_from(const HistItem& hist)
{
   bool cond = _item.addr()==FileMem::nullPtr;
   assert<IsAllocated>(cond,__FILE__,__LINE__);
   cond = hist.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   _mem->allot(_item);
   _fileName = hist.fileName();
   _object = hist.object();
   _command = hist.command();
   _item.timeStamp() = hist.timeStamp();
   _item.fileNamePtr() = _fileName.addr();
   _item.objectPtr() = _object.addr();
   _item.commandPtr() = _command.addr();
   _item.childHead() = rec_add_item(hist._mem,hist.childHead());
   _item.next() = rec_add_item(hist._mem,hist.next());
   _mem->sync(_item,FileSync::write);
}



/// Write any changes made to history item to file memory.
///
/// @exception IsNullPtr This object has not been set(address is nullptr).
void HistItem::sync()
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   _mem->sync(_item,FileSync::write);
}



/// @brief Set time stamp.
///
/// Set this object's time stamp to new value given. This does not effect what
/// is written in file memory.
///
/// @param n New time stamp.
///
/// @exception IsNullPtr This object has not been set(address is nullptr).
void HistItem::timeStamp(int64_t n)
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   _item.timeStamp() = n;
}



/// Get time stamp value of this object.
///
/// @return Time stamp value.
///
/// @exception IsNullPtr This object has not been set(address is nullptr).
int64_t HistItem::timeStamp() const
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   return _item.timeStamp();
}



/// Set this object's file name to value given.
///
/// @param str File name value.
///
/// @exception IsNullPtr This object has not been set(address is nullptr).
/// @exception AlreadySet This object's file name has already been set.
void HistItem::fileName(const std::string& str)
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   cond = _item.fileNamePtr()==FileMem::nullPtr;
   assert<AlreadySet>(cond,__FILE__,__LINE__);
   _fileName = str;
   _item.fileNamePtr() = _fileName.addr();
}



/// Get file name value of this object.
///
/// @return File name value.
///
/// @exception IsNullPtr This object has not been set(address is nullptr).
const HistItem::string& HistItem::fileName() const
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   return *_fileName;
}



/// Set this object's object to value given.
///
/// @param str Object value.
///
/// @exception IsNullPtr This object has not been set(address is nullptr).
/// @exception AlreadySet This object's object has already been set.
void HistItem::object(const std::string& str)
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   cond = _item.objectPtr()==FileMem::nullPtr;
   assert<AlreadySet>(cond,__FILE__,__LINE__);
   _object = str;
   _item.objectPtr() = _object.addr();
}



/// Get object value of this object.
///
/// @return Object value.
///
/// @exception IsNullPtr This object has not been set(address is nullptr).
const HistItem::string& HistItem::object() const
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   return *_object;
}



/// Set this object's command to value given.
///
/// @param str Command value.
///
/// @exception IsNullPtr This object has not been set(address is nullptr).
/// @exception AlreadySet This object's command has already been set.
void HistItem::command(const std::string& str)
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   cond = _item.commandPtr()==FileMem::nullPtr;
   assert<AlreadySet>(cond,__FILE__,__LINE__);
   _command = str;
   _item.commandPtr() = _command.addr();
}



/// Get command value of this object.
///
/// @return Command value.
///
/// @exception IsNullPtr This object has not been set(address is nullptr).
const HistItem::string& HistItem::command() const
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   return *_command;
}



/// Set this object's next pointer to value given.
///
/// @param str Next pointer value.
///
/// @exception IsNullPtr This object has not been set(address is nullptr).
/// @exception AlreadySet This object's next pointer has already been set.
void HistItem::next(FileMem::Ptr ptr)
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   cond = _item.next()==FileMem::nullPtr;
   assert<AlreadySet>(cond,__FILE__,__LINE__);
   _item.next() = ptr;
}



/// Get next pointer of this object.
///
/// @return Next pointer.
///
/// @exception IsNullPtr This object has not been set(address is nullptr).
FileMem::Ptr HistItem::next() const
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   return _item.next();
}



/// Set this object's child pointer to value given.
///
/// @param str Child pointer value.
///
/// @exception IsNullPtr This object has not been set(address is nullptr).
/// @exception AlreadySet This object's child pointer has already been set.
void HistItem::childHead(FileMem::Ptr ptr)
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   cond = _item.childHead()==FileMem::nullPtr;
   assert<AlreadySet>(cond,__FILE__,__LINE__);
   _item.childHead() = ptr;
}



/// Get child pointer of this object.
///
/// @return Child pointer.
///
/// @exception IsNullPtr This object has not been set(address is nullptr).
FileMem::Ptr HistItem::childHead() const
{
   bool cond = _item.addr()!=FileMem::nullPtr;
   assert<IsNullPtr>(cond,__FILE__,__LINE__);
   return _item.childHead();
}



/// @brief Load new address.
///
/// Set file memory location to new location given. If the new location is not
/// nullptr then load history item data from the new location for this object,
/// else leave this object empty and unset.
///
/// @param ptr New file memory location, or nullptr to unset.
void HistItem::operator=(FileMem::Ptr ptr)
{
   _item = ptr;
   _fileName.addr(FileMem::nullPtr);
   _object.addr(FileMem::nullPtr);
   _command.addr(FileMem::nullPtr);
   if (_item.addr()!=FileMem::nullPtr)
   {
      load_item();
   }
}



/// @brief Load item from file memory.
///
/// Loads History item from file memory and populates this object, using this
/// object's file memory pointer as location to load from.
inline void HistItem::load_item()
{
   _mem->sync(_item,FileSync::read);
   if (_item.fileNamePtr()!=FileMem::nullPtr)
   {
      _fileName.addr(_item.fileNamePtr());
   }
   if (_item.objectPtr()!=FileMem::nullPtr)
   {
      _object.addr(_item.objectPtr());
   }
   if (_item.commandPtr()!=FileMem::nullPtr)
   {
      _command.addr(_item.commandPtr());
   }
}



/// @brief Recursively copy memory file location.
///
/// If the file memory location is not nullptr, recursively copies given history
/// item from given file memory object into this object's file memory object.
/// Recursive in this context is defined as also copying all children and any
/// items in the forward list pointed to by the next pointer.
///
/// @param mem File memory object to copy items from.
/// @param ptr File memory location in given object to copy from.
FileMem::Ptr HistItem::rec_add_item(FileMem* mem, FileMem::Ptr ptr)
{
   FileMem::Ptr ret = ptr;
   if (ptr!=FileMem::nullPtr)
   {
      HistItem from(mem,ptr);
      HistItem to(_mem);
      to.copy_from(from);
      ret = to.addr();
   }
   return ret;
}
