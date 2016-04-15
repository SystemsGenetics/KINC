#ifndef ANALYTIC_H
#define ANALYTIC_H
#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>
#include <string>
#include <memory>
#include "dataplugin.h"
#include "terminal.h"


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
   /// @brief Executes analytic on data objects.
   ///
   /// Must execute the analytic method of this class on the given input and
   /// output data objects given to the object before this function is called.
   /// This function must not return control until execution of the analytic on
   /// all data objects are finished and written to.
   ///
   /// @param ops Additional arguments and options passed to the analytic from
   /// the user.
   /// @param tm The terminal interface that can be used by the analytic to
   /// output information about processing.
   /// @param dev An optional OpenCL device that can be used for accelerated
   /// computation of the analytic method. This can be nullptr if the user has
   /// not selected an OpenCL device in the console.
   virtual void execute(GetOpts& ops, Terminal& tm, cl::Device* dev) = 0;
};



#endif
