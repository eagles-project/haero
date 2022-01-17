module merikanto2007

  use haero_precision, only: wp

  implicit none

  private

  public :: log_nucleation_rate, &
            onset_temperature, &
            critical_radius, &
            num_critical_molecules, &
            num_h2so4_molecules, &
            num_nh3_molecules

contains

  !> Computes the logarithm of the ternary nucleation rate [cm-3 s-1] as
  !> parameterized by Merikanto et al (2007), eq 8.
  !> @param [in] temp The atmospherіc temperature [K]
  !> @param [in] rel_hum The relative humidity [-]
  !> @param [in] c_h2so4 The number concentration of H2SO4 gas [cm-3]
  !> @param [in] xi_nh3 The molar mixing ratio of NH3 [ppt]
  function log_nucleation_rate(temp, rel_hum, c_h2so4, xi_nh3) result(rate)
    real(wp), intent(in) :: temp, rel_hum, c_h2so4, xi_nh3

    real(wp) :: rate

    rate = -12.861848898625231_wp &
           + 4.905527742256349_wp*xi_nh3 &
           - 358.2337705052991_wp*rel_hum &
           - 0.05463019231872484_wp*xi_nh3*temp &
           + 4.8630382337426985_wp*rel_hum*temp &
           + 0.00020258394697064567_wp*xi_nh3*temp**2 &
           - 0.02175548069741675_wp*rel_hum*temp**2 &
           - 2.502406532869512e-7_wp*xi_nh3*temp**3 &
           + 0.00003212869941055865_wp*rel_hum*temp**3 &
           - 4.39129415725234e6_wp/log(c_h2so4)**2 &
           + (56383.93843154586_wp*temp)/log(c_h2so4)**2 &
           - (239.835990963361_wp*temp**2)/log(c_h2so4)**2 &
           + (0.33765136625580167_wp*temp**3)/log(c_h2so4)**2 &
           - (629.7882041830943_wp*rel_hum)/(xi_nh3**3*log(c_h2so4)) &
           + (7.772806552631709_wp*rel_hum*temp)/(xi_nh3**3*log(c_h2so4)) &
           - (0.031974053936299256_wp*rel_hum*temp**2)/(xi_nh3**3*log(c_h2so4)) &
           + (0.00004383764128775082_wp*rel_hum*temp**3)/(xi_nh3**3*log(c_h2so4)) &
           + 1200.472096232311_wp*log(c_h2so4) &
           - 17.37107890065621_wp*temp*log(c_h2so4) &
           + 0.08170681335921742_wp*temp**2*log(c_h2so4) &
           - 0.00012534476159729881_wp*temp**3*log(c_h2so4) &
           - 14.833042158178936_wp*log(c_h2so4)**2 &
           + 0.2932631303555295_wp*temp*log(c_h2so4)**2 &
           - 0.0016497524241142845_wp*temp**2*log(c_h2so4)**2 &
           + 2.844074805239367e-6_wp*temp**3*log(c_h2so4)**2 &
           - 231375.56676032578_wp*log(xi_nh3) &
           - 100.21645273730675_wp*rel_hum*log(xi_nh3) &
           + 2919.2852552424706_wp*temp*log(xi_nh3) &
           + 0.977886555834732_wp*rel_hum*temp*log(xi_nh3) &
           - 12.286497122264588_wp*temp**2*log(xi_nh3) &
           - 0.0030511783284506377_wp*rel_hum*temp**2*log(xi_nh3) &
           + 0.017249301826661612_wp*temp**3*log(xi_nh3) &
           + 2.967320346100855e-6_wp*rel_hum*temp**3*log(xi_nh3) &
           + (2.360931724951942e6_wp*log(xi_nh3))/log(c_h2so4) &
           - (29752.130254319443_wp*temp*log(xi_nh3))/log(c_h2so4) &
           + (125.04965118142027_wp*temp**2*log(xi_nh3))/log(c_h2so4) &
           - (0.1752996881934318_wp*temp**3*log(xi_nh3))/log(c_h2so4) &
           + 5599.912337254629_wp*log(c_h2so4)*log(xi_nh3) &
           - 70.70896612937771_wp*temp*log(c_h2so4)*log(xi_nh3) &
           + 0.2978801613269466_wp*temp**2*log(c_h2so4)*log(xi_nh3) &
           - 0.00041866525019504_wp*temp**3*log(c_h2so4)*log(xi_nh3) &
           + 75061.15281456841_wp*log(xi_nh3)**2 &
           - 931.8802278173565_wp*temp*log(xi_nh3)**2 &
           + 3.863266220840964_wp*temp**2*log(xi_nh3)**2 &
           - 0.005349472062284983_wp*temp**3*log(xi_nh3)**2 &
           - (732006.8180571689_wp*log(xi_nh3)**2)/log(c_h2so4) &
           + (9100.06398573816_wp*temp*log(xi_nh3)**2)/log(c_h2so4) &
           - (37.771091915932004_wp*temp**2*log(xi_nh3)**2)/log(c_h2so4) &
           + (0.05235455395566905_wp*temp**3*log(xi_nh3)**2)/log(c_h2so4) &
           - 1911.0303773001353_wp*log(c_h2so4)*log(xi_nh3)**2 &
           + 23.6903969622286_wp*temp*log(c_h2so4)*log(xi_nh3)**2 &
           - 0.09807872005428583_wp*temp**2*log(c_h2so4)*log(xi_nh3)**2 &
           + 0.00013564560238552576_wp*temp**3*log(c_h2so4)*log(xi_nh3)**2 &
           - 3180.5610833308_wp*log(xi_nh3)**3 &
           + 39.08268568672095_wp*temp*log(xi_nh3)**3 &
           - 0.16048521066690752_wp*temp**2*log(xi_nh3)**3 &
           + 0.00022031380023793877_wp*temp**3*log(xi_nh3)**3 &
           + (40751.075322248245_wp*log(xi_nh3)**3)/log(c_h2so4) &
           - (501.66977622013934_wp*temp*log(xi_nh3)**3)/log(c_h2so4) &
           + (2.063469732254135_wp*temp**2*log(xi_nh3)**3)/log(c_h2so4) &
           - (0.002836873785758324_wp*temp**3*log(xi_nh3)**3)/log(c_h2so4) &
           + 2.792313345723013_wp*log(c_h2so4)**2*log(xi_nh3)**3 &
           - 0.03422552111802899_wp*temp*log(c_h2so4)**2*log(xi_nh3)**3 &
           + 0.00014019195277521142_wp*temp**2*log(c_h2so4)**2*log(xi_nh3)**3 &
           - 1.9201227328396297e-7_wp*temp**3*log(c_h2so4)**2*log(xi_nh3)**3 &
           - 980.923146020468_wp*log(rel_hum) &
           + 10.054155220444462_wp*temp*log(rel_hum) &
           - 0.03306644502023841_wp*temp**2*log(rel_hum) &
           + 0.000034274041225891804_wp*temp**3*log(rel_hum) &
           + (16597.75554295064_wp*log(rel_hum))/log(c_h2so4) &
           - (175.2365504237746_wp*temp*log(rel_hum))/log(c_h2so4) &
           + (0.6033215603167458_wp*temp**2*log(rel_hum))/log(c_h2so4) &
           - (0.0006731787599587544_wp*temp**3*log(rel_hum))/log(c_h2so4) &
           - 89.38961120336789_wp*log(xi_nh3)*log(rel_hum) &
           + 1.153344219304926_wp*temp*log(xi_nh3)*log(rel_hum) &
           - 0.004954549700267233_wp*temp**2*log(xi_nh3)*log(rel_hum) &
           + 7.096309866238719e-6_wp*temp**3*log(xi_nh3)*log(rel_hum) &
           + 3.1712136610383244_wp*log(xi_nh3)**3*log(rel_hum) &
           - 0.037822330602328806_wp*temp*log(xi_nh3)**3*log(rel_hum) &
           + 0.0001500555743561457_wp*temp**2*log(xi_nh3)**3*log(rel_hum) &
           - 1.9828365865570703e-7_wp*temp**3*log(xi_nh3)**3*log(rel_hum)
  end function

  !> Computes the "onset temperature" [K] (eq 10) above which Merikanto's
  !> parameterization for the nucleation rate (eq 8) cannot be used (in which
  !> case the authors suggest setting the nucleation rate to zero).
  !> @param [in] rel_hum The relative humidity [-]
  !> @param [in] c_h2so4 The number concentration of H2SO4 gas [cm-3]
  !> @param [in] xi_nh3 The molar mixing ratio of NH3 [ppt]
  function onset_temperature(rel_hum, c_h2so4, xi_nh3) result(t_onset)

    real(wp), intent(in) :: rel_hum, c_h2so4, xi_nh3

    real(wp) :: t_onset

    t_onset = 143.6002929064716_wp &
              + 1.0178856665693992_wp*rel_hum &
              + 10.196398812974294_wp*log(c_h2so4) &
              - 0.1849879416839113_wp*log(c_h2so4)**2 &
              - 17.161783213150173_wp*log(xi_nh3) &
              + (109.92469248546053_wp*log(xi_nh3))/log(c_h2so4) &
              + 0.7734119613144357_wp*log(c_h2so4)*log(xi_nh3) &
              - 0.15576469879527022_wp*log(xi_nh3)**2
  end function

  !> Computes the radius of a critical cluster [nm] as parameterized in Merikanto
  !> et al (2007), eq 11.
  !> @param [in] log_j The logarithm of the nucleation rate ["log (cm-3 s-1)"]
  !> @param [in] temp The atmospherіc temperature [K]
  !> @param [in] c_h2so4 The number concentration of H2SO4 gas [cm-3]
  !> @param [in] xi_nh3 The molar mixing ratio of NH3 [ppt]
  function critical_radius(log_j, temp, c_h2so4, xi_nh3) result(r_crit)
    real(wp), intent(in) :: log_j, temp, c_h2so4, xi_nh3

    real(wp) :: r_crit

    r_crit = 3.2888553966535506e-10_wp &
             - 3.374171768439839e-12_wp*temp &
             + 1.8347359507774313e-14_wp*temp**2 &
             + 2.5419844298881856e-12_wp*log(c_h2so4) &
             - 9.498107643050827e-14_wp*temp*log(c_h2so4) &
             + 7.446266520834559e-13_wp*log(c_h2so4)**2 &
             + 2.4303397746137294e-11_wp*log(xi_nh3) &
             + 1.589324325956633e-14_wp*temp*log(xi_nh3) &
             - 2.034596219775266e-12_wp*log(c_h2so4)*log(xi_nh3) &
             - 5.59303954457172e-13_wp*log(xi_nh3)**2 &
             - 4.889507104645867e-16_wp*temp*log(xi_nh3)**2 &
             + 1.3847024107506764e-13_wp*log(xi_nh3)**3 &
             + 4.141077193427042e-15_wp*log_j &
             - 2.6813110884009767e-14_wp*temp*log_j &
             + 1.2879071621313094e-12_wp*log(xi_nh3)*log_j &
             - 3.80352446061867e-15_wp*temp*log(xi_nh3)*log_j &
             - 1.8790172502456827e-14_wp*log_j**2
  end function

  !> Computes the total number of molecules in a critical cluster as
  !> parameterized in Merikanto et al (2007), eq 12.
  !> @param [in] log_j The logarithm of the nucleation rate ["log (cm-3 s-1)"]
  !> @param [in] temp The atmospherіc temperature [K]
  !> @param [in] c_h2so4 The number concentration of H2SO4 gas [cm-3]
  !> @param [in] xi_nh3 The molar mixing ratio of NH3 [ppt]
  function num_critical_molecules(log_j, temp, c_h2so4, xi_nh3) result(n_crit)

    real(wp), intent(in) :: log_j, temp, c_h2so4, xi_nh3

    real(wp) :: n_crit

    n_crit = 57.40091052369212_wp &
      - 0.2996341884645408_wp*temp &
      + 0.0007395477768531926_wp*temp**2 &
      - 5.090604835032423_wp*log(c_h2so4) &
      + 0.011016634044531128_wp*temp*log(c_h2so4) &
      + 0.06750032251225707_wp*log(c_h2so4)**2 &
      - 0.8102831333223962_wp*log(xi_nh3) &
      + 0.015905081275952426_wp*temp*log(xi_nh3) &
      - 0.2044174683159531_wp*log(c_h2so4)*log(xi_nh3) &
      + 0.08918159167625832_wp*log(xi_nh3)**2 &
      - 0.0004969033586666147_wp*temp*log(xi_nh3)**2 &
      + 0.005704394549007816_wp*log(xi_nh3)**3 &
      + 3.4098703903474368_wp*log_j &
      - 0.014916956508210809_wp*temp*log_j &
      + 0.08459090011666293_wp*log(xi_nh3)*log_j &
      - 0.00014800625143907616_wp*temp*log(xi_nh3)*log_j &
      + 0.00503804694656905_wp*log_j**2
  end function

  !> Computes the total number of H2SO4 molecules in a critical cluster as
  !> parameterized in Merikanto et al (2007), eq 13.
  !> @param [in] log_j The logarithm of the nucleation rate ["log (cm-3 s-1)"]
  !> @param [in] temp The atmospherіc temperature [K]
  !> @param [in] c_h2so4 The number concentration of H2SO4 gas [cm-3]
  !> @param [in] xi_nh3 The molar mixing ratio of NH3 [ppt]
  function num_h2so4_molecules(log_j, temp, c_h2so4, xi_nh3) result(n_h2so4)
    real(wp), intent(in) :: log_j, temp, c_h2so4, xi_nh3

    real(wp) :: n_h2so4

    n_h2so4 = -4.7154180661803595_wp &
      + 0.13436423483953885_wp*temp &
      - 0.00047184686478816176_wp*temp**2 &
      - 2.564010713640308_wp*log(c_h2so4) &
      + 0.011353312899114723_wp*temp*log(c_h2so4) &
      + 0.0010801941974317014_wp*log(c_h2so4)**2 &
      + 0.5171368624197119_wp*log(xi_nh3) &
      - 0.0027882479896204665_wp*temp*log(xi_nh3) &
      + 0.8066971907026886_wp*log(xi_nh3)**2 &
      - 0.0031849094214409335_wp*temp*log(xi_nh3)**2 &
      - 0.09951184152927882_wp*log(xi_nh3)**3 &
      + 0.00040072788891745513_wp*temp*log(xi_nh3)**3 &
      + 1.3276469271073974_wp*log_j &
      - 0.006167654171986281_wp*temp*log_j &
      - 0.11061390967822708_wp*log(xi_nh3)*log_j &
      + 0.0004367575329273496_wp*temp*log(xi_nh3)*log_j &
      + 0.000916366357266258_wp*log_j**2
  end function

  !> Computes the total number of NH3 molecules in a critical cluster as
  !> parameterized in Merikanto et al (2007), eq 14.
  !> @param [in] log_j The logarithm of the nucleation rate ["log (cm-3 s-1)"]
  !> @param [in] temp The atmospherіc temperature [K]
  !> @param [in] c_h2so4 The number concentration of H2SO4 gas [cm-3]
  !> @param [in] xi_nh3 The molar mixing ratio of NH3 [ppt]
  function num_nh3_molecules(log_j, temp, c_h2so4, xi_nh3) result(n_nh3)
    real(wp), intent(in) :: log_j, temp, c_h2so4, xi_nh3

    real(wp) :: n_nh3

    n_nh3 = 71.20073903979772_wp &
      - 0.8409600103431923_wp*temp &
      + 0.0024803006590334922_wp*temp**2 &
      + 2.7798606841602607_wp*log(c_h2so4) &
      - 0.01475023348171676_wp*temp*log(c_h2so4) &
      + 0.012264508212031405_wp*log(c_h2so4)**2 &
      - 2.009926050440182_wp*log(xi_nh3) &
      + 0.008689123511431527_wp*temp*log(xi_nh3) &
      - 0.009141180198955415_wp*log(c_h2so4)*log(xi_nh3) &
      + 0.1374122553905617_wp*log(xi_nh3)**2 &
      - 0.0006253227821679215_wp*temp*log(xi_nh3)**2 &
      + 0.00009377332742098946_wp*log(xi_nh3)**3 &
      + 0.5202974341687757_wp*log_j &
      - 0.002419872323052805_wp*temp*log_j &
      + 0.07916392322884074_wp*log(xi_nh3)*log_j &
      - 0.0003021586030317366_wp*temp*log(xi_nh3)*log_j &
      + 0.0046977006608603395_wp*log_j**2
  end function

end module
