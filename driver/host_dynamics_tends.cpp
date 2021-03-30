#include "host_dynamics_tends.hpp"

namespace haero {
namespace driver {

void DynamicsTendencies::compute(const Real t,
  ConstColumnView phi, ConstColumnView rho, const AtmosphericConditions& conds) {

  const int nlp1 = phi.extent(0);
  const int nl = rho.extent(0);

  Kokkos::parallel_for("DynamicsTendencies::ComputeInterfaces", nlp1,
    KOKKOS_LAMBDA (const int pack_idx) {
      for (int vi=0; vi<PackInfo::vec_end(nlp1, pack_idx); ++vi) {
        w_tend(pack_idx)[vi] = wtend(t, phi(pack_idx)[vi], conds);
        phi_tend(pack_idx)[vi] = phitend(t, phi(pack_idx)[vi], conds);
      }
    });

  Kokkos::parallel_for("DynamicsTendencies::ComputeMidpoints", nl,
    KOKKOS_LAMBDA (const int pack_idx) {
      for (int vi=0; vi<PackInfo::vec_end(nl, pack_idx); ++vi) {
        const int k = PackInfo::array_idx(pack_idx,vi);
        const int kphalf_pack = PackInfo::pack_idx(k+1);
        const int kphalf_vec = PackInfo::vec_idx(k+1);
        const Real phimid = 0.5*(phi(pack_idx)[vi] + phi(kphalf_pack)[kphalf_vec]);

        rho_tend(pack_idx)[vi] = rhotend(t, phimid, rho(pack_idx)[vi], conds);
        thetav_tend(pack_idx)[vi] = thetavtend();
        qv_tend(pack_idx)[vi] = qvtend();
      }
    });
}

} // namespace driver
} // namespace haero
