#ifndef HAERO_MAM_NUCLEATION_PROCESS_HPP
#define HAERO_MAM_NUCLEATION_PROCESS_HPP

#include "haero/process.hpp"
#include "haero/physical_constants.hpp"


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
  double adjust_factor_pbl_ratenucl = 0;

  public:

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

  virtual void init(const ModalAerosolConfig& modal_aerosol_config) override;

  KOKKOS_FUNCTION
  virtual void run(const ModalAerosolConfig& modal_aerosol_config,
                   Real t, Real dt,
                   const Prognostics& prognostics,
                   const Atmosphere& atmosphere,
                   const Diagnostics& diagnostics,
                   Tendencies& tendencies) const override;


  /// Set the Adjustment factor for nucleation rate corrected for the planetary boundary
  /// layer.  This is used in calculating the boundary nucleation rate in
  /// pbl_nuc_wang2008.
  void set_adjust_factor_pbl_ratenucl(const double v) 
  {
     adjust_factor_pbl_ratenucl = v;
  }

/// pbl_nuc_wang2008 calculates boundary nucleation rate
/// using the first or second-order parameterization in
///     wang, m., and j.e. penner, 2008,
///        aerosol indirect forcing in a global model with particle nucleation,
///        atmos. chem. phys. discuss., 8, 13943-13998

/// @param [in]   so4vol            ! concentration of h2so4 (molecules cm-3)
/// @param [in]   newnuc_method_flagaa [11,12] value selects [first,second]-order parameterization

/// @param [inout]  newnuc_method_flagaa2
/// @param [inout]  ratenucl         ! binary nucleation rate, j (# cm-3 s-1)
/// @param [inout]  rateloge         ! log( ratenucl )
/// @param [inout]  cnum_tot         ! total number of molecules

/// in the critical nucleus
/// @param [inout]  cnum_h2so4       ! number of h2so4 molecules
/// @param [inout]  cnum_nh3         ! number of nh3 molecules
/// @param [inout]  radius_cluster   ! the radius of cluster (nm)


KOKKOS_INLINE_FUNCTION
void pbl_nuc_wang2008(const double so4vol,
                      const int    newnuc_method_flagaa,
                      int    & newnuc_method_flagaa2,
                      double & ratenucl,
                      double & rateloge,
                      double & cnum_tot,
                      double & cnum_h2so4,
                      double & cnum_nh3,
                      double & radius_cluster ) const
{
  using namespace std;

  double tmp_ratenucl = 0;
  // nucleation rate
  if (newnuc_method_flagaa == 11)
    tmp_ratenucl = 1.0e-6 * so4vol;
  else if (newnuc_method_flagaa == 12)
    tmp_ratenucl = 1.0e-12 * so4vol*so4vol;
  else
    return;

  tmp_ratenucl = tmp_ratenucl * adjust_factor_pbl_ratenucl;
  const double tmp_rateloge = log( max( 1.0e-38, tmp_ratenucl ) );

  //! exit if pbl nuc rate is lower than (incoming) ternary/binary rate
  if (tmp_rateloge <= rateloge) return;

  rateloge = tmp_rateloge;
  ratenucl = tmp_ratenucl;
  newnuc_method_flagaa2 = newnuc_method_flagaa;

  // following wang 2002, assume fresh nuclei are 1 nm diameter
  //    subsequent code will "grow" them to aitken mode size
  radius_cluster = 0.5;

  // assume fresh nuclei are pure h2so4
  //    since aitken size >> initial size, the initial composition
  //    has very little impact on the results
  const double tmp_diam = radius_cluster * 2.0e-7;                  // diameter in cm
  const double tmp_volu = (tmp_diam*tmp_diam*tmp_diam) * (constants::pi/6.0);  // volume in cm^3
  const double tmp_mass = tmp_volu * 1.8;                           // mass in g
  cnum_h2so4 = (tmp_mass / 98.0) * 6.023e23;                        // no. of h2so4 molec assuming pure h2so4
  cnum_tot = cnum_h2so4;
  cnum_nh3 = 0.0;
}


/// binary_nuc_vehk2002 calculates binary nucleation rate and critical cluster size
/// using the parameterization in
///     vehkamäki, h., m. kulmala, i. napari, k.e.j. lehtinen,
///        c. timmreck, m. noppel and a. laaksonen, 2002,
///        an improved parameterization for sulfuric acid-water nucleation
///        rates for tropospheric and stratospheric conditions,
///        j. geophys. res., 107, 4622, doi:10.1029/2002jd002184

///  @param [in]  temp          temperature (k)
///  @param [in]  rh            relative humidity (0-1)
///  @param [in]  so4vol        concentration of h2so4 (molecules cm-3)
///
///  @param [out] ratenucl      binary nucleation rate, j (# cm-3 s-1)
///  @param [out] rateloge      log( ratenucl )

///  @param [out] cnum_h2so4    number of h2so4 molecules in the critical nucleus
///  @param [out] cnum_tot      total number of molecules in the critical nucleus

KOKKOS_INLINE_FUNCTION
static void binary_nuc_vehk2002(const double temp,
                         const double rh,
                         const double so4vol,
                         double &ratenucl,
                         double &rateloge,
                         double &cnum_h2so4,
                         double &cnum_tot,
                         double &radius_cluster)

{
  using namespace std;

  //calc sulfuric acid mole fraction in critical cluster
  const double crit_x = 0.740997 - 0.00266379 * temp
    - 0.00349998 * log (so4vol)
    + 0.0000504022 * temp * log (so4vol)
    + 0.00201048 * log (rh)
    - 0.000183289 * temp * log(rh)
    + 0.00157407 * log(rh) * log(rh)
    - 0.0000179059 * temp * log(rh) * log(rh)
    + 0.000184403 * log(rh) * log(rh) * log(rh)
    - 1.50345e-6 * temp * log(rh) * log(rh) * log(rh);

  // calc nucleation rate
  double acoe = 0.14309+2.21956*temp
    - 0.0273911 * (temp*temp)
    + 0.0000722811 * (temp*temp*temp) + 5.91822/crit_x;;

  double bcoe = 0.117489 + 0.462532 *temp
    - 0.0118059 * (temp*temp)
    + 0.0000404196 * (temp*temp*temp) + 15.7963/crit_x;

  double ccoe = -0.215554-0.0810269 * temp
    + 0.00143581 * (temp*temp)
    - 4.7758e-6 * (temp*temp*temp)
    - 2.91297/crit_x;

  double dcoe = -3.58856+0.049508 * temp
    - 0.00021382 * (temp*temp)
    + 3.10801e-7 * (temp*temp*temp)
    - 0.0293333/crit_x;

  double ecoe = 1.14598 - 0.600796 * temp
    + 0.00864245 * (temp*temp)
    - 0.0000228947 * (temp*temp*temp)
    - 8.44985/crit_x;

  double fcoe = 2.15855 + 0.0808121 * temp
    - 0.000407382 * (temp*temp)
    - 4.01957e-7 * (temp*temp*temp)
    + 0.721326/crit_x;

  double gcoe = 1.6241 - 0.0160106 * temp
    + 0.0000377124 * (temp*temp)
    + 3.21794e-8 * (temp*temp*temp)
    - 0.0113255/crit_x;

  double hcoe = 9.71682 - 0.115048 * temp
    + 0.000157098 * (temp*temp)
    + 4.00914e-7 * (temp*temp*temp)
    + 0.71186/crit_x;

  double icoe = -1.05611 + 0.00903378 * temp
    - 0.0000198417 * (temp*temp)
    + 2.46048e-8  * (temp*temp*temp)
    - 0.0579087/crit_x;

  double jcoe = -0.148712 + 0.00283508 * temp
    - 9.24619e-6  * (temp*temp)
    + 5.00427e-9 * (temp*temp*temp)
    - 0.0127081/crit_x;

  double tmpa = (
    acoe
    + bcoe * log (rh)
    + ccoe * log (rh) * log(rh)
    + dcoe * log (rh) * log(rh) * log(rh)
    + ecoe * log (so4vol)
    + fcoe * log (rh) * log (so4vol)
    + gcoe * log (rh) * log(rh)
    * (log (so4vol))
    + hcoe * log (so4vol) * log (so4vol)
    + icoe * log (rh)
    * log (so4vol) * log (so4vol)
    + jcoe * log (so4vol) * log (so4vol) * log (so4vol)
  );
  rateloge = tmpa;
  tmpa = min( tmpa, log(1.0e38) );
  ratenucl = exp ( tmpa );

  // calc number of molecules in critical cluster
  acoe = -0.00295413 - 0.0976834*temp
    + 0.00102485 * (temp*temp)
    - 2.18646e-6 * (temp*temp*temp) - 0.101717/crit_x;

  bcoe = -0.00205064 - 0.00758504*temp
    + 0.000192654 * (temp*temp)
    - 6.7043e-7 * (temp*temp*temp) - 0.255774/crit_x;

  ccoe = +0.00322308 + 0.000852637 * temp
    - 0.0000154757 * (temp*temp)
    + 5.66661e-8 * (temp*temp*temp)
    + 0.0338444/crit_x;

  dcoe = +0.0474323 - 0.000625104 * temp
    + 2.65066e-6 * (temp*temp)
    - 3.67471e-9 * (temp*temp*temp)
    - 0.000267251/crit_x;

  ecoe = -0.0125211 + 0.00580655 * temp
    - 0.000101674 * (temp*temp)
    + 2.88195e-7 * (temp*temp*temp)
    + 0.0942243/crit_x;

  fcoe = -0.038546 - 0.000672316 * temp
    + 2.60288e-6 * (temp*temp)
    + 1.19416e-8 * (temp*temp*temp)
    - 0.00851515/crit_x;

  gcoe = -0.0183749 + 0.000172072 * temp
    - 3.71766e-7 * (temp*temp)
    - 5.14875e-10 * (temp*temp*temp)
    + 0.00026866/crit_x;

  hcoe = -0.0619974 + 0.000906958 * temp
    - 9.11728e-7 * (temp*temp)
    - 5.36796e-9 * (temp*temp*temp)
    - 0.00774234/crit_x;

  icoe = +0.0121827 - 0.00010665 * temp
    + 2.5346e-7 * (temp*temp)
    - 3.63519e-10 * (temp*temp*temp)
    + 0.000610065/crit_x;

  jcoe = +0.000320184 - 0.0000174762 * temp
    + 6.06504e-8 * (temp*temp)
    - 1.4177e-11 * (temp*temp*temp)
    + 0.000135751/crit_x;

  cnum_tot = exp (
    acoe
    + bcoe * log (rh)
    + ccoe * log (rh) * log (rh)
    + dcoe * log (rh) * log (rh) * log (rh)
    + ecoe * log (so4vol)
    + fcoe * log (rh) * log (so4vol)
    + gcoe * log (rh) * log (rh)
    * (log (so4vol))
    + hcoe * log (so4vol) * log (so4vol)
    + icoe * log (rh)
    * log (so4vol) * log (so4vol)
    + jcoe * log (so4vol) * log (so4vol) * log (so4vol)
  );

  cnum_h2so4 = cnum_tot * crit_x;

  // calc radius (nm) of critical cluster
  radius_cluster = exp( -1.6524245 + 0.42316402*crit_x
    + 0.3346648*log(cnum_tot) );
}

/// ternary_nuc_merik2007 calculates the parameterized composition
/// and nucleation rate of critical clusters in h2o-h2so4-nh3 vapor
///
/// warning: the fit should not be used outside its limits of validity
/// (limits indicated below)
///
/// @param [in]  t:     temperature (k), limits 235-295 k
/// @param [in] rh:    relative humidity as fraction (eg. 0.5=50%) limits 0.05-0.95
/// @param [in] c2:    sulfuric acid concentration (molecules/cm3) limits 5x10^4 - 10^9 molecules/cm3
/// @param [in] c3:    ammonia mixing ratio (ppt) limits 0.1 - 1000 ppt
///
/// @param [out] j_log: logarithm of nucleation rate (1/(s cm3))
/// @param [out] ntot:  total number of molecules in the critical cluster
/// @param [out] nacid: number of sulfuric acid molecules in the critical cluster
/// @param [out] namm:  number of ammonia molecules in the critical cluster
/// @param [out] r:     radius of the critical cluster (nm)
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
      // nucleation rate less than 5e-6, setting j_log arbitrary small
      j_log = -300.0 ;
    }
  }
};


}

#endif
