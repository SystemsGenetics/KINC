#ifndef TERMINAL_H
#define TERMINAL_H
#include <string>
#include "exception.h"



/// @brief Terminal interface class.
///
/// Defines the behavior of the terminal interface this program uses as basic
/// console based input and output. Some functions are not virtual, either be
/// direct or indirect wraparounds to virtual functions.
///
/// @author Josh Burns
/// @date 2 February 2016
class Terminal
{
public:
   // *
   // * EXCEPTIONS
   // *
   struct Exception;
   struct SystemError;
   struct InvalidUse;
   // *
   // * BASIC METHODS
   // *
   virtual ~Terminal() = default;
   // *
   // * STATIC FUNCTIONS
   // *
   static Terminal& newline(Terminal&);
   static Terminal& endl(Terminal&);
   static Terminal& flush(Terminal&);
   static Terminal& general(Terminal&);
   static Terminal& warning(Terminal&);
   static Terminal& error(Terminal&);
   // *
   // * FUNCTIONS
   // *
   Terminal& print(short);
   Terminal& print(unsigned short);
   Terminal& print(int);
   Terminal& print(unsigned int);
   Terminal& print(long);
   Terminal& print(unsigned long);
   Terminal& print(float);
   Terminal& print(double);
   Terminal& print(const char*);
   Terminal& print(const std::string&);
   // *
   // * OPERATORS
   // *
   Terminal& operator<<(Terminal& (*pf)(Terminal&));
   // *
   // * VIRTUAL FUNCTIONS
   // *
   /// Must set the header of the terminal interface to the string provided.
   ///
   /// @param newHeader New header
   virtual void header(const std::string& newHeader) = 0;
   /// @brief Set precision for floating point numbers.
   ///
   /// Must set the precision for all proceeding floating point outputs to the
   /// precision given.
   ///
   /// @param prec The precision defined by the number of digits displayed past
   /// the point.
   virtual void precision(int prec) = 0;
   // Output operators, self-explanatory.
   virtual Terminal& operator<<(short n) = 0;
   virtual Terminal& operator<<(unsigned short n) = 0;
   virtual Terminal& operator<<(int n) = 0;
   virtual Terminal& operator<<(unsigned int n) = 0;
   virtual Terminal& operator<<(long n) = 0;
   virtual Terminal& operator<<(unsigned long n) = 0;
   virtual Terminal& operator<<(float n) = 0;
   virtual Terminal& operator<<(double n) = 0;
   virtual Terminal& operator<<(const char* str) = 0;
   virtual Terminal& operator<<(const std::string& str) = 0;
   /// @brief Reads line from user and returns result.
   ///
   /// Must read a single line of input from the user terminal. Must also be
   /// blocking until a new line of input has been read.
   ///
   /// @param input Will be overwritten with new line of input from user.
   virtual void operator>>(std::string& input) = 0;
protected:
   // *
   // * ENUMERATIONS
   // *
   /// Defines all possible operators for terminal output.
   enum class Ops
   {
      newline, ///< Make new line to output.
      flush, ///< Flush any buffered output.
      general, ///< Set proceeding output as general message.
      warning, ///< Set proceeding output as warning message.
      error ///< Set proceeding output as error message.
   };
   // *
   // * VIRTUAL FUNCTIONS
   // *
   /// @brief Outputs operator type.
   ///
   /// Must take any operator from the enumerated type Ops and process it to
   /// output of the terminal. Some operators change state information which
   /// must be saved after call of this function.
   ///
   /// @param oper Specific operator being given to output of terminal.
   virtual void set_ops(Ops oper) = 0;
};



/// Generic base exception class for all exceptions thrown in Terminal class.
struct Terminal::Exception : public ::Exception
{
   using ::Exception::Exception;
};

/// Exception that is thrown when a system error occurs.
struct Terminal::SystemError : public ::SystemError
{
   using ::SystemError::SystemError;
};

/// Exception that is thrown when an implementation of Terminal is used in an
/// invalid way.
struct Terminal::InvalidUse : public Terminal::Exception
{
   InvalidUse(const char* file, int line):
      Exception(file,line,"Terminal::InvalidUse")
   {}
};



#endif
