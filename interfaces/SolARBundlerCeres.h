#ifndef SOLARBUNDLERCERES_H
#define SOLARBUNDLERCERES_H


#include "api/solver/map/IBundler.h"
#include "api/storage/ICovisibilityGraph.h"
#include "api/storage/IKeyframesManager.h"
#include "api/storage/IPointCloudManager.h"

#include "xpcf/component/ConfigurableBase.h"

#include "SolARCeresAPI.h"


#include "ceres/ceres.h"
#include "ceres/rotation.h"

#define POINT_DIM 3
#define OBSERV_DIM 2
#define EXT_DIM 6
#define INT_DIM 9



namespace xpcf = org::bcom::xpcf;

namespace SolAR {
	using namespace datastructure;
    using namespace api::storage;
	namespace MODULES {
		namespace CERES {

			/**
			 * @class SolARBundlerCeres
			 * @brief <B>Applies a bundle adjustment to optimize a 3D map and keyframes.</B>
			 * <TT>UUID: 4897fc13-682c-4e95-8aba-abd9f7a17193</TT>
			 *
             * This Bundler component has to inject 3 storage components:
             *         -IKeyframesManager
             *         -IPointCloudManager
             *         -ICovisibilityGraph
             * The definition of the injection of this three storage components will have to be added in the xml configuration file of your application.
			 */

			class SOLARCERES_EXPORT_API SolARBundlerCeres : public org::bcom::xpcf::ConfigurableBase,
				public api::solver::map::IBundler {
			public:
				SolARBundlerCeres();
				~SolARBundlerCeres() = default;

				//     org::bcom::xpcf::XPCFErrorCode onConfigured() override final;
				void unloadComponent() override final;

				/// @brief solve a non-linear problem related to bundle adjustement statement expressed as:
				/// minArg(pts3ds,intrinsics,extrinsics) = MIN_cam_i(MIN_3d_j(pts2d_j - reproje(pt3ds_j,intrinsics_i,extrinsics_i)),
				/// @param[in, out] K: camera calibration parameters responsible of 3D points generation.
				/// @param[in, out] D: camera distorsion parameters responsible of 3D points generation
				/// @param[in] selectKeyframes : selected views to bundle following a given strategies. If it is empty then take all keyframes into account to perform global bundle adjustment.
				/// @param[in] useSpanningTree: in the case of the global bundle adjustment, if it true, the optimization is based on a maximal spanning tree.
				/// @return the mean re-projection error after optimization.
				double bundleAdjustment(CamCalibration & K, CamDistortion & D, const std::vector<uint32_t> & selectKeyframes, const bool & useSpanningTree = false) override;

			private:
                /// @brief number of mx iterations number.
                unsigned int m_iterationsNo = 10;
                /// @brief fixing map control.
                unsigned int m_fixedMap = 0;
                /// @brief fixing extrinsic control.
                unsigned int m_fixedKeyframes = 0;
                /// @brief fixing intrinsic control.
                unsigned int m_fixedIntrinsics = 1;
                /// @brief fixing first pose control.
                unsigned int m_fixedFirstPose = 1;
                /// @brief fixing neighbour keyframes control.
                unsigned int m_fixedNeighbourKeyframes = 1;
                /// @brief Maximum number of fixed neighbour keyframes.
                unsigned int m_nbMaxFixedKeyframes = 100;

                /// @brief reference to the storage component use to manage the point cloud.
                SRef<IPointCloudManager>    m_pointCloudManager;
                /// @brief reference to the storage component use to manage the keyframes.
                SRef<IKeyframesManager>     m_keyframesManager;
				/// @brief reference to the storage component use to manage the covisibility graph.
				SRef<ICovisibilityGraph>     m_covisibilityGraph;

				/// @brief transform a rotation matrix to axis-anle representation using Rodrigue's formula.
				/// @param[in]  R:              a pose transform matrix
				/// @param[out] r:             rodrigues angles
				/// keeps translation matrix unchanged
				inline void toRodrigues(Transform3Df &R, Vector3f&r) {
					const double small = 1e-6;
					double th = acos
					(0.5*(fmax(R(0, 0) + R(1, 1) + R(2, 2), -1.0) - 1.0));
					double sth = sin(th);
					double cth = cos(th);
					if (fabs(sth) < small && cth < 0) {
						double W_pt[9], x, y, z;
						W_pt[0] = 0.5*(R(0, 0) + R(0, 0)) - 1.0;
						W_pt[1] = 0.5*(R(1, 0) + R(0, 1));
						W_pt[2] = 0.5*(R(2, 0) + R(0, 2));

						W_pt[3] = 0.5*(R(0, 1) + R(1, 0));
						W_pt[4] = 0.5*(R(1, 1) + R(1, 1)) - 1.0;
						W_pt[5] = 0.5*(R(2, 1) + R(1, 2));

						W_pt[6] = 0.5*(R(0, 2) + R(2, 0));
						W_pt[7] = 0.5*(R(1, 2) + R(2, 1));
						W_pt[8] = 0.5*(R(2, 2) + R(2, 2)) - 1.0;


						x = sqrt(0.5 * (W_pt[0] - W_pt[4] - W_pt[8]));
						y = sqrt(0.5 * (W_pt[4] - W_pt[8] - W_pt[0]));
						z = sqrt(0.5 * (W_pt[8] - W_pt[0] - W_pt[4]));


						if (x >= y && x >= z) {
							y = (W_pt[1] >= 0) ? y : -y;
							z = (W_pt[2] >= 0) ? z : -z;
						}
						else if (y >= x && y >= z) {
							z = (W_pt[5] >= 0) ? z : -z;
							x = (W_pt[1] >= 0) ? x : -x;
						}
						else {
							x = (W_pt[2] >= 0) ? x : -x;
							y = (W_pt[5] >= 0) ? y : -y;
						}

						{
							double scale = th / sqrt(1 - cth);
							r(0) = scale * x;
							r(1) = scale * y;
							r(2) = scale * z;

							return;
						}

					}
					else {
						double a = (fabs(sth) < small) ? 1 : th / sin(th);
						//    double b ;
						r(0) = 0.5*a*(R(2, 1) - R(1, 2));
						r(1) = 0.5*a*(R(0, 2) - R(2, 0));
						r(2) = 0.5*a*(R(1, 0) - R(0, 1));
					}

				}
				/// @brief transform  axis-angles representation to rotation matrix using inverse Rodrigue's formula.
				/// @param[in]  r:             Rodrigue's angles
				/// @param[out] R:             rotation matrix
				/// initialize translation matrix to zero
				inline void iRodrigues(Vector3d &r, Transform3Df&R) {

					const double small = 1e-6;

					double th = sqrt(r(0)*r(0) +
						r(1)*r(1) +
						r(2)*r(2));

					if (th < small) {
						R(0, 0) = 1.0; R(0, 1) = 0.0; R(0, 2) = 0.0;
						R(1, 0) = 0.0; R(1, 1) = 1.0; R(1, 2) = 0.0;
						R(2, 0) = 0.0; R(2, 1) = 0.0; R(2, 2) = 1.0;
						return;
					}

					{
						double x = r(0) / th;
						double y = r(1) / th;
						double z = r(2) / th;

						double xx = x * x;
						double xy = x * y;
						double xz = x * z;
						double yy = y * y;
						double yz = y * z;
						double zz = z * z;

						const double yx = xy;
						const double zx = xz;
						const double zy = yz;

						double sth = sin(th);
						double cth = cos(th);
						double mcth = 1.0 - cth;

						R(0, 0) = 1 - mcth * (yy + zz);
						R(1, 0) = sth * z + mcth * xy;
						R(2, 0) = -sth * y + mcth * xz;

						R(0, 1) = -sth * z + mcth * yx;
						R(1, 1) = 1 - mcth * (zz + xx);
						R(2, 1) = sth * x + mcth * yz;

						R(0, 2) = sth * y + mcth * xz;
						R(1, 2) = -sth * x + mcth * yz;
						R(2, 2) = 1 - mcth * (xx + yy);

					}
				}
			};
		}
	}
}

#endif // SOLARBUNDLERCERES_H
