#ifndef SOLARBUNDLERCERES_H
#define SOLARBUNDLERCERES_H


#include "api/solver/map/IBundler.h"
#include "xpcf/component/ComponentBase.h"
#include <vector>
#include "SolARCeresAPI.h"





#include <string>

namespace SolAR {
    using namespace datastructure;
    namespace MODULES {
        namespace CERES {
            class SOLARCERES_EXPORT_API SolARBundlerCeres : public org::bcom::xpcf::ComponentBase,
                public api::solver::map::IBundler {
            public:
                SolARBundlerCeres();
                ~SolARBundlerCeres() = default;
                bool adjustBundle(const std::string&path_bundle,
                                  std::vector<SRef<CloudPoint>>& cloud_before,
                                  std::vector<SRef<CloudPoint>>& cloud_after) override final;
                void unloadComponent () override final;

            };
        }
    }
}



#endif // SOLARBUNDLERCERES_H
