/*
 * See GPL.txt for copyright information.
 *
 * Author: Joshua Burns
 *
 */
#ifndef CLDEVLIST_H
#define CLDEVLIST_H
#include <vector>



class CLDevice;



/// @brief Container list of all OpenCL devices.
///
/// List of all available OpenCL devices. The list is constructed upon
/// instantiation of the object and can be recompiled by called the refresh()
/// function.
class CLDevList
{
   class Iterator;
public:
   // ****************************** Basic Methods **************************
   CLDevList();
   // ****************************** Functions ******************************
   Iterator begin();
   Iterator end();
   void refresh();
   bool exist(int,int);
   CLDevice& at(int,int);
private:
   void build();
   std::vector<std::vector<CLDevice>> _list;
};



/// Iterate through OpenCL device list.
///
/// @pre CLDevList object this iterate references cannot change during lifetime
/// of an iterator.
class CLDevList::Iterator
{
   friend class CLDevList;
public:
   // ****************************** Operators ******************************
   CLDevice& operator*();
   void operator++();
   void operator--();
   bool operator!=(const Iterator&);
private:
   // ****************************** Basic Methods **************************
   Iterator(int,int,CLDevList*);
   CLDevList* _devList;
   int _pi;
   int _di;
};



#endif
