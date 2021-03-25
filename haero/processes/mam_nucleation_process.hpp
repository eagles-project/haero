#ifndef HAERO_MAM_NUCLEATION_PROCESS_HPP
#define HAERO_MAM_NUCLEATION_PROCESS_HPP

#include "haero/process.hpp"


namespace haero {


/// @class MAMNucleationProcess
/// This type is an interface (base class) to an aerosol process quantified by
/// a parametrization that computes a set of tendencies for prognostic variables
/// within an aerosol system. Each subclass of this type implements a particular
/// **implementation** of a specific **parametrization** for a particular **process**.
///
/// To make these ideas more complete, consider the following examples of
/// important **physical processes** in the ProcessType above.
///
/// Each of these processes has one or more **parametrizations**--mathematical
/// models that quantify the outcomes of these processes in specific
/// circumstances. For example, **surface emissions** may be parametrized by:
/// * time series input from a file
/// * a low-order polynomial in the time variable
/// * estimated from, e.g., an agent-based representation of human activity.
///
/// Finally, each of these parametrizations can have one or more
/// **implementations**. For example, every parametrization can have a Fortran
/// implementation that runs only on CPUs, as well as a C++ implementation that
/// can run on CPUs and GPUs.
///
/// The PrognosticProcess class provides an interface for all
/// implementations of all parametrizations for all physical processes that
/// compute tendencies for aerosol systems.
class MAMNucleationProcess : public PrognosticProcess {
  public:

  /// Constructor, called by all PrognosticProcess subclasses.
  /// @param [in] type The type of aerosol process modeled by the subclass.
  /// @param [in] name A descriptive name that captures the aerosol process,
  ///                  its underlying parametrization, and its implementation.
  MAMNucleationProcess();

  /// Destructor.
  KOKKOS_INLINE_FUNCTION
  virtual ~MAMNucleationProcess() {}

  /// Default copy constructor. For use in moving host instance to device.
  KOKKOS_INLINE_FUNCTION
  MAMNucleationProcess(const MAMNucleationProcess& pp) :
    PrognosticProcess(pp) {}

  /// MAMNucleationProcess objects are not assignable.
  PrognosticProcess& operator=(const MAMNucleationProcess&) = delete;

  //------------------------------------------------------------------------
  //                                Accessors
  //------------------------------------------------------------------------


  //------------------------------------------------------------------------
  //                Methods to be overridden by subclasses
  //------------------------------------------------------------------------

  /// Override this method if your aerosol process needs to be initialized
  /// with information about the model. The default implementation does nothing.
  /// @param [in] modal_aerosol_config The aerosol configuration describing the
  ///                                  aerosol system to which this process
  ///                                  belongs.
  virtual void init(const ModalAerosolConfig& modal_aerosol_config);

  /// Override this method to implement the aerosol process using the specific
  /// parameterization for the subclass.
  /// @param [in] modal_aerosol_config The aerosol configuration describing the
  ///                                  aerosol system to which this process
  ///                                  belongs.
  /// @param [in] t The simulation time at which this process is being invoked
  ///               (in seconds).
  /// @param [in] dt The simulation time interval ("timestep size") over which
  ///                this process occurs.
  /// @param [in] prognostics The prognostic variables used by and affected by
  ///                         this process.
  /// @param [in] atmosphere The atmosphere state variables used by this process.
  /// @param [in] diagnostics The prognostic variables used by and affected by
  ///                         this process.
  /// @param [out] tendencies A container that stores time derivatives for
  ///                         prognostic variables evolved by this process.
  KOKKOS_FUNCTION
  virtual void run(const ModalAerosolConfig& modal_aerosol_config,
                   Real t, Real dt,
                   const Prognostics& prognostics,
                   const Atmosphere& atmosphere,
                   const Diagnostics& diagnostics,
                   Tendencies& tendencies) const;


  KOKKOS_INLINE_FUNCTION 
  static void ternary_nuc_merik2007(const double t, 
                                    const double rh, 
                                    const double c2,
                                    const double c3, 
                                    double &j_log, 
                                    double &ntot, 
                                    double &nacid, 
                                    double &namm, 
                                    double &r)
 {
    using namespace std;


    const double t_onset=143.6002929064716 + 1.0178856665693992*rh + \
      10.196398812974294*log(c2) - \
      0.1849879416839113*(log(c2)*log(c2)) - 17.161783213150173*log(c3) + \
      (109.92469248546053*log(c3))/log(c2) + \
      0.7734119613144357*log(c2)*log(c3) - 0.15576469879527022*(log(c3)*log(c3));

    if (t_onset > t) {
      j_log = -12.861848898625231 + 4.905527742256349*c3 - 358.2337705052991*rh - \
        0.05463019231872484*c3*t + 4.8630382337426985*rh*t + \
        0.00020258394697064567*c3*(t*t) - 0.02175548069741675*rh*(t*t) - \
        2.502406532869512e-7*c3*(t*t*t) + 0.00003212869941055865*rh*(t*t*t) - \
        4.39129415725234e6/(log(c2)*log(c2)) + (56383.93843154586*t)/(log(c2)*log(c2)) - \
        (239.835990963361*(t*t))/(log(c2)*log(c2)) + \
        (0.33765136625580167*(t*t*t))/(log(c2)*log(c2)) - \
        (629.7882041830943*rh)/((c3*c3*c3)*log(c2)) + \
        (7.772806552631709*rh*t)/((c3*c3*c3)*log(c2)) - \
        (0.031974053936299256*rh*(t*t))/((c3*c3*c3)*log(c2)) + \
        (0.00004383764128775082*rh*(t*t*t))/((c3*c3*c3)*log(c2)) + \
        1200.472096232311*log(c2) - 17.37107890065621*t*log(c2) + \
        0.08170681335921742*(t*t)*log(c2) - 0.00012534476159729881*(t*t*t)*log(c2) - \
        14.833042158178936*(log(c2)*log(c2)) + 0.2932631303555295*t*(log(c2)*log(c2)) - \
        0.0016497524241142845*(t*t)*(log(c2)*log(c2)) + \
        2.844074805239367e-6*(t*t*t)*(log(c2)*log(c2)) - 231375.56676032578*log(c3) - \
        100.21645273730675*rh*log(c3) + 2919.2852552424706*t*log(c3) + \
        0.977886555834732*rh*t*log(c3) - 12.286497122264588*(t*t)*log(c3) - \
        0.0030511783284506377*rh*(t*t)*log(c3) + \
        0.017249301826661612*(t*t*t)*log(c3) + 2.967320346100855e-6*rh*(t*t*t)*log(c3) + \
        (2.360931724951942e6*log(c3))/log(c2) - \
        (29752.130254319443*t*log(c3))/log(c2) + \
        (125.04965118142027*(t*t)*log(c3))/log(c2) - \
        (0.1752996881934318*(t*t*t)*log(c3))/log(c2) + \
        5599.912337254629*log(c2)*log(c3) - 70.70896612937771*t*log(c2)*log(c3) + \
        0.2978801613269466*(t*t)*log(c2)*log(c3) - \
        0.00041866525019504*(t*t*t)*log(c2)*log(c3) + 75061.15281456841*(log(c3)*log(c3)) - \
        931.8802278173565*t*(log(c3)*log(c3)) + 3.863266220840964*(t*t)*(log(c3)*log(c3)) - \
        0.005349472062284983*(t*t*t)*(log(c3)*log(c3)) - \
        (732006.8180571689*(log(c3)*log(c3)))/log(c2) + \
        (9100.06398573816*t*(log(c3)*log(c3)))/log(c2) - \
        (37.771091915932004*(t*t)*(log(c3)*log(c3)))/log(c2) + \
        (0.05235455395566905*(t*t*t)*(log(c3)*log(c3)))/log(c2) - \
        1911.0303773001353*log(c2)*(log(c3)*log(c3)) + \
        23.6903969622286*t*log(c2)*(log(c3)*log(c3)) - \
        0.09807872005428583*(t*t)*log(c2)*(log(c3)*log(c3)) + \
        0.00013564560238552576*(t*t*t)*log(c2)*(log(c3)*log(c3)) - \
        3180.5610833308*(log(c3)*log(c3)*log(c3)) + 39.08268568672095*t*(log(c3)*log(c3)*log(c3)) - \
        0.16048521066690752*(t*t)*(log(c3)*log(c3)*log(c3)) + \
        0.00022031380023793877*(t*t*t)*(log(c3)*log(c3)*log(c3)) + \
        (40751.075322248245*(log(c3)*log(c3)*log(c3)))/log(c2) - \
        (501.66977622013934*t*(log(c3)*log(c3)*log(c3)))/log(c2) + \
        (2.063469732254135*(t*t)*(log(c3)*log(c3)*log(c3)))/log(c2) - \
        (0.002836873785758324*(t*t*t)*(log(c3)*log(c3)*log(c3)))/log(c2) + \
        2.792313345723013*(log(c2)*log(c2))*(log(c3)*log(c3)*log(c3)) - \
        0.03422552111802899*t*(log(c2)*log(c2))*(log(c3)*log(c3)*log(c3)) + \
        0.00014019195277521142*(t*t)*(log(c2)*log(c2))*(log(c3)*log(c3)*log(c3)) - \
        1.9201227328396297e-7*(t*t*t)*(log(c2)*log(c2))*(log(c3)*log(c3)*log(c3)) - \
        980.923146020468*log(rh) + 10.054155220444462*t*log(rh) - \
        0.03306644502023841*(t*t)*log(rh) + 0.000034274041225891804*(t*t*t)*log(rh) + \
        (16597.75554295064*log(rh))/log(c2) - \
        (175.2365504237746*t*log(rh))/log(c2) + \
        (0.6033215603167458*(t*t)*log(rh))/log(c2) - \
        (0.0006731787599587544*(t*t*t)*log(rh))/log(c2) - \
        89.38961120336789*log(c3)*log(rh) + 1.153344219304926*t*log(c3)*log(rh) - \
        0.004954549700267233*(t*t)*log(c3)*log(rh) + \
        7.096309866238719e-6*(t*t*t)*log(c3)*log(rh) + \
        3.1712136610383244*(log(c3)*log(c3)*log(c3))*log(rh) - \
        0.037822330602328806*t*(log(c3)*log(c3)*log(c3))*log(rh) + \
        0.0001500555743561457*(t*t)*(log(c3)*log(c3)*log(c3))*log(rh) - \
        1.9828365865570703e-7*(t*t*t)*(log(c3)*log(c3)*log(c3))*log(rh);

      const double j = exp(j_log);

      ntot = 57.40091052369212 - 0.2996341884645408*t + \
        0.0007395477768531926*(t*t) - \
        5.090604835032423*log(c2) + 0.011016634044531128*t*log(c2) + \
        0.06750032251225707*(log(c2)*log(c2)) - 0.8102831333223962*log(c3) + \
        0.015905081275952426*t*log(c3) - 0.2044174683159531*log(c2)*log(c3) + \
        0.08918159167625832*(log(c3)*log(c3)) - 0.0004969033586666147*t*(log(c3)*log(c3)) + \
        0.005704394549007816*(log(c3)*log(c3)*log(c3)) + 3.4098703903474368*log(j) - \
        0.014916956508210809*t*log(j) + 0.08459090011666293*log(c3)*log(j) - \
        0.00014800625143907616*t*log(c3)*log(j) + 0.00503804694656905*(log(j)*log(j));

      r = 3.2888553966535506e-10 - 3.374171768439839e-12*t + \
        1.8347359507774313e-14*(t*t) + 2.5419844298881856e-12*log(c2) - \
        9.498107643050827e-14*t*log(c2) + 7.446266520834559e-13*(log(c2)*log(c2)) + \
        2.4303397746137294e-11*log(c3) + 1.589324325956633e-14*t*log(c3) - \
        2.034596219775266e-12*log(c2)*log(c3) - 5.59303954457172e-13*(log(c3)*log(c3)) - \
        4.889507104645867e-16*t*(log(c3)*log(c3)) + 1.3847024107506764e-13*(log(c3)*log(c3)*log(c3)) + \
        4.141077193427042e-15*log(j) - 2.6813110884009767e-14*t*log(j) + \
        1.2879071621313094e-12*log(c3)*log(j) - \
        3.80352446061867e-15*t*log(c3)*log(j) - 1.8790172502456827e-14*(log(j)*log(j));

      nacid = -4.7154180661803595 + 0.13436423483953885*t - \
        0.00047184686478816176*(t*t) - \
        2.564010713640308*log(c2) + 0.011353312899114723*t*log(c2) + \
        0.0010801941974317014*(log(c2)*log(c2)) + 0.5171368624197119*log(c3) - \
        0.0027882479896204665*t*log(c3) + 0.8066971907026886*(log(c3)*log(c3)) - \
        0.0031849094214409335*t*(log(c3)*log(c3)) - 0.09951184152927882*(log(c3)*log(c3)*log(c3)) + \
        0.00040072788891745513*t*(log(c3)*log(c3)*log(c3)) + 1.3276469271073974*log(j) - \
        0.006167654171986281*t*log(j) - 0.11061390967822708*log(c3)*log(j) + \
        0.0004367575329273496*t*log(c3)*log(j) + 0.000916366357266258*(log(j)*log(j));

      namm = 71.20073903979772 - 0.8409600103431923*t + \
        0.0024803006590334922*(t*t) + \
        2.7798606841602607*log(c2) - 0.01475023348171676*t*log(c2) + \
        0.012264508212031405*(log(c2)*log(c2)) - 2.009926050440182*log(c3) + \
        0.008689123511431527*t*log(c3) - 0.009141180198955415*log(c2)*log(c3) + \
        0.1374122553905617*(log(c3)*log(c3)) - 0.0006253227821679215*t*(log(c3)*log(c3)) + \
        0.00009377332742098946*(log(c3)*log(c3)*log(c3)) + 0.5202974341687757*log(j) - \
        0.002419872323052805*t*log(j) + 0.07916392322884074*log(c3)*log(j) - \
        0.0003021586030317366*t*log(c3)*log(j) + 0.0046977006608603395*(log(j)*log(j));

    } else {
      // nucleation rate less that 5e-6, setting j_log arbitrary small
      j_log = -300.0 ;
    }
  }
};


}

#endif
