#ifndef CHEMUTIL_HPP
#define CHEMUTIL_HPP

#include <string>
#include <sstream>
#include <iostream>

// these are required for TChem
// #include "tchem/TChem_SourceTermToyProblem.hpp"
#include "tchem/TChem_KineticModelData.hpp"
#include "tchem/TChem_Util.hpp"
#include "tchem/TChem_Impl_RateOfProgress.hpp"

namespace chemUtil {

// some aliases
using ordinal_type = TChem::ordinal_type;
using real_type = TChem::real_type;
using real_type_1d_view = TChem::real_type_1d_view;
using real_type_2d_view = TChem::real_type_2d_view;
using policy_type = typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;

// these are the file paths for the inputs and output Tchem requires
class chemFiles{
  public:
    std::string prefixPath;
    std::string chemFile;
    std::string thermFile;
    std::string outputFile;
    chemFiles(std::string chemDir);
};

// sets up everything necessary for TChem Toy Problem and gets results (tendencies)
class chemSolver{
  private:
    // path to input files
    std::string chemDir;
    // chemFiles object containing file paths for inputs/output
    chemFiles cfiles;
    // verbose output flag (required, but doesn't seem to actually do much)
    bool verbose;
    // number of chemical batches/samples that will have the given composition
    int nBatch;
    // timer variables
    real_type t_deepcopy, t_device_batch;
    // kokkos team policy
    policy_type policy;
    // timer object for kokkos
    Kokkos::Impl::Timer timer;
    // the below are required for the call to get_results()
    // *-----------------------------------------------------------------------*
    // lat/lon views
    real_type_1d_view theta, lambda;
    // 2D views containing chemical state, and omega for tendency results
    real_type_2d_view state, omega;
    // TChem kinetic model data
    TChem::KineticModelData kmd;
    // a const version of the object that contains the data describing the
      // kinetic model
    TChem::KineticModelConstData<TChem::exec_space> kmcd;
    // *-----------------------------------------------------------------------*
    // prints a summary of the TChem run
    void print_summary(const chemFiles& cfiles);
  public:
    // constructor
    chemSolver(std::string chemDir, bool detail, int inBatch, bool iverbose,
               real_type itheta, real_type ilambda,
               real_type initX, real_type initX2);
    // runs chemical model and saves the results (tendencies) to the output file
    real_type_2d_view get_results();
};

} // namespace chemUtil

// NOTE: everything in this namespace is copied directly from TChem
namespace tchem_stuff {

  struct SourceTermToyProblem
  {

    template<typename KineticModelConstDataType>
    static inline ordinal_type getWorkSpaceSize(
      const KineticModelConstDataType& kmcd)
    {
      return SourceTermToyProblem::getWorkSpaceSize(kmcd);
    }

    //
    static void runDeviceBatch( /// input
      typename UseThisTeamPolicy<exec_space>::type& policy,
      const real_type_1d_view& theta,
      const real_type_1d_view& lambda,
      const real_type_2d_view& state,
      /// output
      const real_type_2d_view& SourceTermToyProblem,
      /// const data from kinetic model
      const KineticModelConstDataDevice& kmcd);

  }; // end STTP struct

  template<typename KineticModelConstDataType>
  KOKKOS_INLINE_FUNCTION static ordinal_type getWorkSpaceSize(
    const KineticModelConstDataType& kmcd);
} // end tchem_stuff namespace
#endif
