#include "haero/processes/mam_calcsize_process.hpp"

#include <cmath>

namespace haero {

MAMCalcsizeProcess::MAMCalcsizeProcess()
    : DeviceAerosolProcess<MAMCalcsizeProcess>(CalcsizeProcess, "MAMCalcsizeProcess") {}
}  // namespace haero
