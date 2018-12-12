#ifndef SOLARBUNDLERCERES_H
#define SOLARBUNDLERCERES_H


#include "api/solver/map/IBundler.h"
#include "xpcf/component/ConfigurableBase.h"
#include "SolARCeresAPI.h"

#include "ceres/ceres.h"
#include "ceres/rotation.h"

#define POINT_DIM 3
#define OBSERV_DIM 2
#define EXT_DIM 6
#define INT_DIM 9



namespace xpcf  = org::bcom::xpcf;

namespace SolAR {
    using namespace datastructure;
    namespace MODULES {
        namespace CERES {
            class SOLARCERES_EXPORT_API SolARBundlerCeres : public org::bcom::xpcf::ConfigurableBase,
                public api::solver::map::IBundler {
            public:
                SolARBundlerCeres();
                ~SolARBundlerCeres() = default;

           //     org::bcom::xpcf::XPCFErrorCode onConfigured() override final;
                void unloadComponent () override final;

                /// @brief solve a non-linear problem related to bundle adjustement statement expressed as:
                /// minArg(pts3ds,intrinsics,extrinsics) = MIN_cam_i(MIN_3d_j(pts2d_j - reproje(pt3ds_j,intrinsics_i,extrinsics_i)),
                /// @param[in] framesToAdjust: contains a set of {2D points, camera extrinsics}.
                /// @param[in] mapToAjust: contains a set of of 3D points .
                /// @param[in] K: camera calibration parameters responsible of 3D points generation.
                /// @param[in] D: camera distorsion parameters responsible of 3D points generation
                /// K, D represent the camera intrinsic parameters
                /// @return[in] selectKeyframes : selected views to bundle following a given strategies (ex: poseGraph).
                /// @return the mean re-projection error after {pts3d, intrinsic, extrinsic} correction.                
               double solve(std::vector<SRef<Keyframe>>&framesToAdjust,
                            std::vector<SRef<CloudPoint>>&mapToAdjust,
                            CamCalibration &K,
                            CamDistortion &D,
                            const std::vector<int>&selectKeyframes) override;


            private :

               /// @brief initialize ceres solver parameters (iterations number, solver type, solver strategy..).
                void initCeresProblem();
                /// @brief solve ceres problem using a pre-definded solver parameters: 03 residual blocks are considered
                ///     a) mutable_extrinsic_for_observation(j): correspondant extrinsic parameters for a given observation.
                ///     b) mutable_intrinsic_for_observation(j): correspondant intrinsic parameters for a given observation.
                ///     c) mutable_point_for_observation(j):     correspondant pose for a given observation.
                /// Each residual encodes a reprojection error ERROR(REPROJ(X_j, EXTR_j, INTR_j),x_) following SolARReprojectionError.
                double solveCeresProblem();
                /// @brief fill ceres internal data
                /// @param[in] framesToAdjust: fills mutable_exteinsic_for_observation buffer.
                /// @param[in] mapToAjust:     fills mutable_point_for_observation buffer .
                /// @param[in] K:              fills mutable_intrinsic_for_observation buffer.
                /// @param[in] D:              fills mutable_intrinsic_for_observation buffer
                void fillCeresProblem(std::vector<SRef<Keyframe>>&framesToAdjust,
                                      std::vector<SRef<CloudPoint>>&mapToAdjust,
                                      CamCalibration &K,
                                      CamDistortion &D,
                                      const std::vector<int>&selectedKeyframes);


                /// @brief update all ceres problem variables.
                /// @param[in] framesToAdjust: takes extrinsic parameters correction.
                /// @param[in] mapToAjust:     takes 3D point correction .
                /// @param[in] K:              takes intrinsic parameters correction.
                /// @param[in] D:              takes intrinsic parameters correction.
                void updateCeresProblem(std::vector<SRef<Keyframe>>&framesToAdjust,
                                        std::vector<SRef<CloudPoint>>&mapToAdjust,
                                        CamCalibration &K,
                                        CamDistortion &D);





                /// @brief update 3D point variable.
                /// @param[in] mapToAjust:     takes 3D point correction .
                void updateMap(std::vector<SRef<CloudPoint>>&mapToAdjust);
                /// @brief update extrinsic parameters variable.
                /// @param[in] framesToAdjust: takes extrinsic parameters correction.
                void updateExtrinsic(std::vector<SRef<Keyframe>>&framesToAdjust);

                /// @brief update intrinsic parameters variable.
                /// @param[in] K:              takes intrinsic parameters correction.
                /// @param[in] D:              takes intrinsic parameters correction.
                void updateIntrinsic(CamCalibration &K,CamDistortion &D);
                /// @brief transform a rotation matrix to axis-anle representation using Rodrigue's formula.
                /// @param[in]  R:              a pose transform matrix
                /// @param[out] r:             rodrigues angles
                /// keeps translation matrix unchanged
                inline void toRodrigues(Transform3Df &R, Vector3f&r){
                    const double small = 1e-6 ;
                    double th = acos
                    (0.5*(fmax(R(0,0)+R(1,1)+R(2,2),-1.0) - 1.0)) ;
                    double sth = sin(th) ;
                    double cth = cos(th) ;
                    if(fabs(sth) < small && cth < 0) {
                        double W_pt [9], x, y, z ;
                        W_pt[0] = 0.5*( R(0,0) + R(0,0) ) - 1.0 ;
                        W_pt[1] = 0.5*( R(1,0) + R(0,1) ) ;
                        W_pt[2] = 0.5*( R(2,0) + R(0,2) );

                        W_pt[3] = 0.5*( R(0,1) + R(1,0) );
                        W_pt[4] = 0.5*( R(1,1) + R(1,1) ) - 1.0;
                        W_pt[5] = 0.5*( R(2,1) + R(1,2) );

                        W_pt[6] =  0.5*( R(0,2) + R(2,0) ) ;
                        W_pt[7] =  0.5*( R(1,2) + R(2,1) ) ;
                        W_pt[8] =  0.5*( R(2,2) + R(2,2) ) - 1.0 ;


                        x = sqrt( 0.5 * (W_pt[0]-W_pt[4]-W_pt[8]) ) ;
                        y = sqrt( 0.5 * (W_pt[4]-W_pt[8]-W_pt[0]) ) ;
                        z = sqrt( 0.5 * (W_pt[8]-W_pt[0]-W_pt[4]) ) ;


                        if( x >= y && x >= z ) {
                            y = (W_pt[1] >=0) ? y : -y ;
                            z = (W_pt[2] >=0) ? z : -z ;
                        } else if( y >= x && y >= z ) {
                            z = (W_pt[5] >=0) ? z : -z ;
                            x = (W_pt[1] >=0) ? x : -x ;
                        } else {
                            x = (W_pt[2] >=0) ? x : -x ;
                            y = (W_pt[5] >=0) ? y : -y ;
                        }

                        {
                            double scale = th / sqrt( 1 - cth ) ;
                            r(0) = scale * x ;
                            r(1) = scale * y ;
                            r(2) = scale * z ;

                            return ;
                        }

                    } else {
                        double a = (fabs(sth) < small) ? 1 : th/sin(th) ;
                    //    double b ;
                        r(0) = 0.5*a*(R(2,1) - R(1,2)) ;
                        r(1) = 0.5*a*(R(0,2) - R(2,0)) ;
                        r(2) = 0.5*a*(R(1,0) - R(0,1)) ;
                    }

                }
                /// @brief transform  axis-angles representation to rotation matrix using inverse Rodrigue's formula.
                /// @param[in]  r:             Rodrigue's angles
                /// @param[out] R:             rotation matrix
                /// initialize translation matrix to zero
                inline void iRodrigues(Vector3d &r, Transform3Df&R){

                    const double small = 1e-6 ;

                    double th = sqrt(r(0)*r(0) +
                                     r(1)*r(1) +
                                     r(2)*r(2) ) ;

                    if( th < small ) {
                        R(0,0) = 1.0 ; R(0,1) = 0.0 ; R(0,2) = 0.0 ;
                        R(1,0) = 0.0 ; R(1,1) = 1.0 ; R(1,2) = 0.0 ;
                        R(2,0) = 0.0 ; R(2,1) = 0.0 ; R(2,2) = 1.0 ;
                        return ;
                    }

                    {
                        double x = r(0) / th ;
                        double y = r(1) / th ;
                        double z = r(2) / th ;

                        double xx = x*x ;
                        double xy = x*y ;
                        double xz = x*z ;
                        double yy = y*y ;
                        double yz = y*z ;
                        double zz = z*z ;

                        const double yx = xy ;
                        const double zx = xz ;
                        const double zy = yz ;

                        double sth  = sin(th) ;
                        double cth  = cos(th) ;
                        double mcth = 1.0 - cth ;

                        R(0,0) = 1          - mcth * (yy+zz) ;
                        R(1,0) =     sth*z  + mcth * xy ;
                        R(2,0) =   - sth*y  + mcth * xz ;

                        R(0,1) =   - sth*z  + mcth * yx ;
                        R(1,1) = 1          - mcth * (zz+xx) ;
                        R(2,1) =     sth*x  + mcth * yz ;

                        R(0,2) =     sth*y  + mcth * xz ;
                        R(1,2) =   - sth*x  + mcth * yz ;
                        R(2,2) = 1          - mcth * (xx+yy) ;

                    }
                }

                /// @return get the number of observations.
                inline int num_observations()       const {
                    return m_observationsNo;
                }
                /// @return get the number of parameters.
                inline int get_parameters()       const {
                    return m_parametersNo;
                }
                /// @return get the number of views.
                inline int get_cameras()       const {
                    return m_camerasNo;
                }
                /// @return get the number of points.
                inline int get_points()       const {
                    return m_pointsNo;
                }
                /// @return get tall the observations of the problem as a buffer.
                inline const double* observations() const {
                    return m_observations;
                }
                /// @return get all the parameters of the problem as a buffer.
                double* mutable_cameras() {
                    return m_parameters;
                }
                /// @return get all the 3D points of the problem as a buffer.
                double* mutable_points() {
                    return m_parameters + (INT_DIM + EXT_DIM) * m_camerasNo;
                }
                /// @return get all intrinsic parameters of the problem as a buffer.
                double* mutable_intrinsic() {
                    return m_parameters + (EXT_DIM) * m_camerasNo;
                }
                /// @return get the extrinsic parameters for a given observation_i as a buffer.
                double * mutable_extrinsic_for_observation(int i) {
                    return mutable_cameras() + m_extrinsicIndex[i] * EXT_DIM;
                }
                /// @return get the intrinsic parameters for a given observation_i as a buffer.
                double * mutable_intrinsic_for_observation(int i) {
                    return mutable_intrinsic() + m_intrinsicIndex[i] * INT_DIM;
                }
                /// @return get the 3D point for a given observation_i as a buffer.
                double* mutable_point_for_observation(int i) {
                    return mutable_points() + m_pointIndex[i] * POINT_DIM;
                }

                /// @brief ceres problem containing residual blocks.
                ceres::Problem m_problem;
                /// @brief ceres problem options containing solver parameters.
                ceres::Solver::Options m_options;
                /// @brief ceres problem summary containing problem minimization evolution.
                ceres::Solver::Summary m_summary;
                /// @brief number of views.
                int m_camerasNo;
                /// @brief number of 3D points.
                int m_pointsNo;
                 /// @brief number of 2D observations.
                int m_observationsNo;
                 /// @brief number of residual parameters.
                int m_parametersNo;
                 /// @brief 3D points indices buffer.
                int* m_pointIndex;
                /// @brief views indices buffer.
                int* m_cameraIndex;
                /// @brief extrinsic parameters buffer.
                int * m_extrinsicIndex;
                /// @brief intrinsic parameters buffer.
                int * m_intrinsicIndex;
                /// @brief observations buffer.
                double* m_observations;
                /// @brief residual parameters buffer.
                double* m_parameters;
                /// @brief number of mx iterations number.
                unsigned int m_iterationsNo = 10;
                /// @brief fixing map control.
                unsigned int m_fixedMap = 0;
                /// @brief fixing extrinsic control.
                unsigned int m_fixedExtrinsics = 0;
                /// @brief fixing intrinsic control.
                unsigned int m_fixedIntrinsics = 1;
                /// @brief fixing first pose control.
                unsigned int m_holdFirstPose = 1;
            };
        }
    }
}

#endif // SOLARBUNDLERCERES_H
