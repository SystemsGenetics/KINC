#ifndef DATA_H
#define DATA_H
#include "terminal.h"
#include "getopts.h"



/// Defines the interface between the main program and any data object class.
/// Any data object is instantiated when the file is opened and destroyed when
/// the file is closed.
///
/// @author Josh Burns
/// @date 18 March 2016
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
