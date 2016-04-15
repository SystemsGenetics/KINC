#ifndef GETOPTS_H
#define GETOPTS_H
#include <string>
#include <sstream>
#include <list>
#include <initializer_list>
#include "exception.h"



/// @ingroup dataplugin
/// @brief Decomposes line of user input.
///
/// Takes a single line from the user and seperates it into commands and
/// options. A command is any word seperated by either spaces or tabs that does
/// not begin with a '-' character. An option is any word that begins with a '-'
/// character, having a key and optional value. They key and value of an option
/// are seperated using the '=' character. For example, 'one two --option=val'
/// has two commands(one and two) and one option(with a key of option and a
/// value of val).
class GetOpts
{
public:
   // *
   // * EXCEPTIONS
   // *
   struct Exception;
   struct InvalidType;
   struct InvalidSyntax;
   // *
   // * DECLERATIONS
   // *
   class Iterator;
   using string = std::string;
   using initlist = std::initializer_list<string>;
   // *
   // * BASIC METHODS
   // *
   GetOpts(const string&);
   // *
   // * FUNCTIONS
   // *
   const string& orig() const;
   int com_size() const;
   bool com_empty() const;
   int com_get(initlist);
   string& com_front();
   void com_pop();
   int size() const;
   bool empty() const;
   bool has_opt(const string&,bool = false);
   Iterator begin();
   Iterator end();
   Iterator erase(Iterator);
private:
   // *
   // * DECLERATIONS
   // *
   using clist = std::list<string>;
   using oplist = std::list<std::pair<string,string>>;
   // *
   // * VARIABLES
   // *
   /// Origional line of user input.
   string _orig;
   /// List of commands.
   clist _comms;
   /// List of options.
   oplist _opts;
};



/// Iterates through list of all options decomposed from user input.
///
/// @author Josh Burns
/// @date 24 March 2016
class GetOpts::Iterator
{
public:
   // *
   // * DECLERATIONS
   // *
   friend class GetOpts;
   using string = GetOpts::string;
   using oplist = GetOpts::oplist;
   // *
   // * FUNCTIONS
   // *
   const string& key() const;
   template<class T> T value() const;
   bool is_key(const string&) const;
   bool val_empty() const;
   // *
   // * OPERATORS
   // *
   void operator++();
   bool operator!=(const Iterator&);
   bool operator==(const Iterator&);
private:
   // *
   // * BASIC METHODS
   // *
   Iterator(const oplist::iterator&);
   // *
   // * VARIABLES
   // *
   /// Internal iterator that points to an option in GetOpts list.
   oplist::iterator _i;
};



/// Get origional user input line for this object.
///
/// @return Origional user input.
inline const GetOpts::string& GetOpts::orig() const
{
   return _orig;
}



/// Get total number of commands left for this object.
///
/// @return Number of commands.
inline int GetOpts::com_size() const
{
   return _comms.size();
}



/// Tests if the commands list for this object is empty.
///
/// @return True if commands list is empty, else false.
inline bool GetOpts::com_empty() const
{
   return _comms.empty();
}



/// Returns reference to string of command in front of list for this object.
///
/// @return Front of list command.
inline GetOpts::string& GetOpts::com_front()
{
   return _comms.front();
}



/// Get total number of options left for this object.
///
/// @return Number of options.
inline int GetOpts::size() const
{
   return _opts.size();
}



/// Tests if the options list for this object is empty.
///
/// @return True if options list is empty, else false.
inline bool GetOpts::empty() const
{
   return _opts.empty();
}



/// Gets iterator for options list that points at the beginning of list.
///
/// @return Start of list options iterator.
inline GetOpts::Iterator GetOpts::begin()
{
   return {_opts.begin()};
}



/// Gets iterator for options list that points to one past end of list.
///
/// @return One past end of list iterator.
inline GetOpts::Iterator GetOpts::end()
{
   return {_opts.end()};
}



/// Removes option from list by iterator given.
///
/// @param i Iterator that points to option that will be removed.
/// @return Iterator that points to next option in list or one past end of list
/// if the removed option was at back of list.
inline GetOpts::Iterator GetOpts::erase(Iterator i)
{
   return {_opts.erase(i._i)};
}



/// Get key value of this option.
///
/// @return Key value.
inline const GetOpts::string& GetOpts::Iterator::key() const
{
   return _i->first;
}



/// Get value of this option, if any.
///
/// @param T The variable type for option's value.
/// @return The option's value.
///
/// @exception InvalidType The given type does not match the option's value or
/// there is no value set for this option.
template<class T> T GetOpts::Iterator::value() const
{
   T ret;
   using istring = std::istringstream;
   istring buffer(_i->second);
   bool cond = buffer >> ret;
   assert<InvalidType>(cond,__FILE__,__LINE__);
   return ret;
}



/// Tests if the given key value matches the key value of this option.
///
/// @return True if they match, else fail.
inline bool GetOpts::Iterator::is_key(const string& cmp) const
{
   return _i->first==cmp;
}



/// Tests if the value for this option is empty or not.
///
/// @return True if there is no value set, else false.
inline bool GetOpts::Iterator::val_empty() const
{
   return _i->second.empty();
}



/// Iterates to next option in list of options.
inline void GetOpts::Iterator::operator++()
{
   ++_i;
}



/// Tests if this option iterator is not equal to the given option iterator.
///
/// @return True if they are not equal, else fail.
inline bool GetOpts::Iterator::operator!=(const Iterator& cmp)
{
   return _i!=cmp._i;
}



/// Tests if this option iterator is equal to the given option iterator.
///
/// @return True if they are equal, else fail.
inline bool GetOpts::Iterator::operator==(const Iterator& cmp)
{
   return _i==cmp._i;
}



/// Sets internal list iterator to iterator given.
///
/// @param i Internal iterator used to store pointer to option in list of
/// iterators.
inline GetOpts::Iterator::Iterator(const oplist::iterator& i):
   _i(i)
{}



/// Generic base exception class for all exceptions thrown in GetOpts class.
struct GetOpts::Exception : public ::Exception
{
   using ::Exception::Exception;
};

/// The type given for fetching the value of an option does not match the value
/// of the option or the option's value is not set.
struct GetOpts::InvalidType : public GetOpts::Exception
{
   InvalidType(const char* file, int line):
      Exception(file,line,"GetOpts::InvalidType")
   {}
};

/// A syntax error was encountered while decomposing the user's input line.
struct GetOpts::InvalidSyntax : public GetOpts::Exception
{
   InvalidSyntax(const char* file, int line):
      Exception(file,line,"GetOpts::InvalidSyntax")
   {}
};



#endif
