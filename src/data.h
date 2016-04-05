#ifndef DATA_H
#define DATA_H
#include "terminal.h"
#include "getopts.h"



/// @defgroup dataplugin Data Plugin Framework
/// @brief Basic framework supplied to all data plug-in implementations.
///
/// This is a the basic framework provided for data plug-in implementations. The
/// virtual Data class provides the virtual functions that a data plug-in class
/// is required to implement. The DataPlugin class is the base class which all
/// data plug-in classes are required to have as their public parent. The
/// KincFile class provides the interface to the file memory object which all
/// data plug-in classes use to interact with the binary file associated with
/// their object. Lastly, the FileMem class and its subclasses give the basic
/// file input/output paradigm used to directly interact and read from or write
/// to the binary file attached to data plug-in objects. FString is also
/// provided as a helper class for storing character strings if FileMem binary
/// files.
///
/// To create a data plug-in implementation, make a folder in plugins/data
/// directory. There must be a Makefile inside the root of your new directory,
/// an the directory name must match the name of the class(case sensitive
/// match). The Makefile must compile all source into object files and then
/// store all static objects into a static library to the location provided by
/// the variable TARGET.
///
/// The data plug-in class itself must provide overloaded functions for all
/// virtual functions given in the Data class, along with providing a single
/// constructor that takes the same arguments as the constructor of DataPlugin.
/// You can optionally just pull in the constructor of DataPlugin if desired.
///
/// All exception handling in data plug-in classes should be handled through the
/// base exception class DataException. Throwing any other exception class that
/// does not inherit from this base class will result in undefined behavior.
///
/// Arguments and commands from the user are decoded in the GetOpts class and a
/// reference is given to all data plug-in virtual functions that implement
/// user commands. For the same functions, a reference is also given to the
/// Terminal interface which operates as the output to the user screen the data
/// plug-in can use to display information to the user.
/// @{ @}
///
/// @ingroup dataplugin
/// Defines the interface between the main program and any data object class.
/// Any data object is instantiated when the file is opened and destroyed when
/// the file is closed.
class Data
{
public:
   virtual ~Data() = default;
   /// Executes data load command.
   ///
   /// Must execute this as a user load command with arguments given, printing
   /// any output to terminal interface given.
   ///
   /// @param ops User command line object.
   /// @param tm Terminal interface for output.
   virtual void load(GetOpts& ops, Terminal& tm) = 0;
   /// Executes data dump command.
   ///
   /// Must execute this as a user dump command with arguments given, printing
   /// any output to terminal interface given.
   ///
   /// @param ops User command line object.
   /// @param tm Terminal interface for output.
   virtual void dump(GetOpts& ops, Terminal& tm) = 0;
   /// Executes data query command.
   ///
   /// Must execute this as a user query command with arguments given, printing
   /// any output to terminal interface given.
   ///
   /// @param ops User command line object.
   /// @param tm Terminal interface for output.
   virtual void query(GetOpts& ops, Terminal& tm) = 0;
   /// Must return if the data object is empty or not.
   ///
   /// @return True if data object is empty else false.
   virtual bool empty() = 0;
};



#endif
