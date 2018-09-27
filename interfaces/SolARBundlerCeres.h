#ifndef SOLARBUNDLERCERES_H
#define SOLARBUNDLERCERES_H


#include "api/solver/map/IBundler.h"
#include "xpcf/component/ComponentBase.h"
#include <vector>
#include "SolARCeresAPI.h"

#include "ceres/ceres.h"
#include "ceres/rotation.h"




#include <string>

#define POINT_DIM 3
#define CAM_DIM 9
#define OBSERV_DIM 2


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

               bool adjustBundle(std::vector<SRef<Keyframe>>&framesToAdjust,
                                 std::vector<SRef<CloudPoint>>&mapToAdjust,
                                 const CamCalibration &K,
                                 const CamDistortion &D,
                                 const std::vector<int>&selectKeyframes) override final;


                void unloadComponent () override final;
                bool saveBundleProblem(std::string&path_ba) override final;

            private :

                void initCeresProbelm();
                bool solveCeresProblem();
                void fillCeresProblem(std::vector<SRef<Keyframe>>&framesToAdjust,
                                      std::vector<SRef<CloudPoint>>&mapToAdjust,
                                      const CamCalibration &K,
                                      const CamDistortion &D,
                                      const std::vector<int>&selectedKeyframes);

                bool fillCeresProblemFromFile(const std::string&path_ba);

                bool updateCeresProblem(std::vector<SRef<Keyframe>>&framesToAdjust,
                                        std::vector<SRef<CloudPoint>>&mapToAdjust,
                                        const std::vector<int>&selectedKeyframes);

                bool updateMap(std::vector<SRef<CloudPoint>>&mapToAdjust);
                bool updateExtrinsic(std::vector<SRef<Keyframe>>&framesToAdjust,
                                     const std::vector<int>&selectedKeyframes);
                bool updateIntrinsic(std::vector<SRef<Keyframe>>&framesToAdjust);



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
                        double b ;
                        r(0) = 0.5*a*(R(2,1) - R(1,2)) ;
                        r(1) = 0.5*a*(R(0,2) - R(2,0)) ;
                        r(2) = 0.5*a*(R(1,0) - R(0,1)) ;
                    }

                }

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


                inline int num_observations()       const {
                    return m_observationsNo;
                }
                inline int get_parameters()       const {
                    return m_parametersNo;
                }
                inline int get_cameras()       const {
                    return m_camerasNo;
                }
                inline int get_points()       const {
                    return m_pointsNo;
                }
                inline const double* observations() const {
                    return m_observations;
                }
                double* mutable_cameras() {
                    return m_parameters;
                }
                double* mutable_points() {
                    return m_parameters + 9 * m_camerasNo;
                }

                double* mutable_camera_for_observation(int i) {
                    return mutable_cameras() + m_cameraIndex[i] * 9;
                }
                double* mutable_point_for_observation(int i) {
                    return mutable_points() + m_pointIndex[i] * 3;
                }

                template<typename T>
                void FscanfOrDie(FILE *fptr, const char *format, T *value) {
                    int num_scanned = fscanf(fptr, format, value);
                    if (num_scanned != 1) {
                       std::cerr << "Invalid UW data file."<<std::endl;
                    }
                }

                ceres::Problem m_problem;
                ceres::Solver::Options m_options;

                double m_pHeight;
                double m_pWidth;

                int m_camerasNo;
                int m_pointsNo;
                int m_observationsNo;
                int m_parametersNo;
                int* m_pointIndex;
                int* m_cameraIndex;
                double* m_observations;
                double* m_parameters;
            };
        }
    }
}



#endif // SOLARBUNDLERCERES_H
