#ifndef ANALYTIC_H
#define ANALYTIC_H
#include <CL/cl.hpp>
#include <string>
#include <memory>
#include "dataplugin.h"
#include "terminal.h"
#include "cldevice.h"
#include "clprogram.h"


/// @brief Analytic interface class
///
/// Defines the interface between the program and any analytic plugin class.
/// This virtual class is instantiated for running a specific analytic command
/// and destroyed once the command has executed.
///
/// @author Josh Burns
/// @date 17 March 2016
class Analytic
{
public:
   virtual ~Analytic() = default;
   /// @brief Input data object.
   ///
   /// Must take a single input data object that will be used in execution. All
   /// inputs will be given through this function before execution is called.
   /// The pointer will not be deleted by the analytic object.
   ///
   /// @param in Pointer to the data object to be used as input.
   virtual void input(DataPlugin* in) = 0;
   /// @brief Output data object.
   ///
   /// Must take a single output data object that will be used in execution. All
   /// outputs will be given through this function before execution is called.
   /// The pointer will not be deleted by the analytic object.
   ///
   /// @param out Pointer to the data object to be used as output.
   virtual void output(DataPlugin* out) = 0;
protected:
   virtual void execute_cl(GetOpts& ops, Terminal& tm) = 0;
   virtual void execute_pn(GetOpts& ops, Terminal& tm) = 0;
};



#endif
