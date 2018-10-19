
#include "SolARBundlerCeres.h"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <utility>
#include <string>


#define POINT_DIM 3
#define CAM_DIM 9
#define OBSERV_DIM 2




using namespace std;
namespace xpcf  = org::bcom::xpcf;


XPCF_DEFINE_FACTORY_CREATE_INSTANCE(SolAR::MODULES::CERES::SolARBundlerCeres)

namespace SolAR {
    using namespace datastructure;
        namespace MODULES {
            namespace CERES {


                struct SnavelyReprojectionError {
                    SnavelyReprojectionError(double observed_x, double observed_y)
                        : observed_x(observed_x), observed_y(observed_y) {}
                    template <typename T>
                    bool operator()(const T* const camera,
                        const T* const point,
                        T* residuals) const {

                        // camera[0,1,2] are the angle-axis rotation.
                        T p[3];

                        ceres::AngleAxisRotatePoint(camera, point, p);
                        // camera[3,4,5] are the translation.
                        p[0] += camera[3];
                        p[1] += camera[4];
                        p[2] += camera[5];
                        // Compute the center of distortion. The sign change comes from
                        // the camera model that Noah Snavely's Bundler assumes, whereby
                        // the camera coordinate system has a negative z axis.
                        T xp = +p[0] / p[2];
                        T yp = +p[1] / p[2];
                        // Apply second and fourth order radial distortion.
                        const T& l1 = camera[7];
                        const T& l2 = camera[8];
                        T r2 = xp*xp + yp*yp;
                        T distortion = 1.0 + r2  * (l1 + l2  * r2);
                        // Compute final projected point position.
                        const T& focal = camera[6];


                  //      std::cout<<"    camera parameters  f: "<<focal<<"  l1: "<<l1<<"  l2: "<<l2<<std::endl;

                        T predicted_x = focal * distortion * xp;
                        T predicted_y = focal * distortion * yp;

                        // The error is the difference between the predicted and observed position.
                        residuals[0] = predicted_x - observed_x;
                        residuals[1] = predicted_y - observed_y;


                        return true;
                    }
                    // Factory to hide the construction of the CostFunction object from
                    // the client code.
                    static ceres::CostFunction* Create(const double observed_x,
                        const double observed_y) {
                        return (new ceres::AutoDiffCostFunction<SnavelyReprojectionError, 2, 9, 3>(
                            new SnavelyReprojectionError(observed_x, observed_y)));
                    }
                    double observed_x;
                    double observed_y;
                };
                struct ceresObserv{
                    int cIdx;
                    int pIdx;
                    Point2Df oPt;
                    void show(){
                        std::cout<<" obervation: "<<std::endl;
                        std::cout<<"    # cam idx: "<<cIdx<<" #3d: "<<pIdx<<" #2d: "<<oPt.getX()<<" "<<oPt.getY()<<std::endl;
                    }
                    ceresObserv(){

                    };
                };



                SolARBundlerCeres::SolARBundlerCeres():ComponentBase(xpcf::toUUID<SolARBundlerCeres>())
                {
                     addInterface<IBundler>(this);
                     LOG_DEBUG(" SolARBundlerCeres constructor");
                }


                bool SolARBundlerCeres::adjustBundle(std::vector<SRef<Keyframe>>&framesToAdjust,
                                                     std::vector<SRef<CloudPoint>>&mapToAdjust,
                                                     const CamCalibration &K,
                                                     const CamDistortion &D,
                                                     const std::vector<int>&selectKeyframes){
                    std::cout<<"0->init ceres problem"<<std::endl;
                    initCeresProblem();
                    std::cout<<"1->fill ceres problem"<<std::endl;
                    fillCeresProblem(framesToAdjust,
                                     mapToAdjust,
                                     K,
                                     D,
                                     selectKeyframes);
                    std::cout<<"2->apply ceres problem"<<std::endl;
                    solveCeresProblem();
                    std::cout<<"3->update ceres problem"<<std::endl;
                    updateCeresProblem(framesToAdjust,
                                       mapToAdjust,
                                       selectKeyframes);
                    return true;
                }
                void SolARBundlerCeres::initCeresProblem(){                    
                    m_options.linear_solver_type = ceres::ITERATIVE_SCHUR;
                    m_options.preconditioner_type = ceres::SCHUR_JACOBI;
                    m_options.dense_linear_algebra_library_type = ceres::LAPACK;
                    m_options.max_num_iterations = 20;
                    m_options.num_threads = 1;
                    m_options.num_linear_solver_threads = 1;
                    m_options.logging_type = ceres::SILENT;
                }
                void SolARBundlerCeres::fillCeresProblem(std::vector<SRef<Keyframe>>&framesToAdjust,
                                                         std::vector<SRef<CloudPoint>>&mapToAdjust,
                                                         const CamCalibration &K,
                                                         const CamDistortion &D,
                                                         const std::vector<int>&selectedKeyframes){



                   bool global_bundling = true;
                    if( selectedKeyframes.size() > 0){
                        global_bundling = false;
                    }
                    std::vector<ceresObserv>observations_temp;
                    if(global_bundling){
                        std::cout<<" global bundler "<<std::endl;
                        for(int i = 0; i < mapToAdjust.size(); ++i){
                            std::map<unsigned int, unsigned int> visibility = mapToAdjust[i]->getVisibility();
                            int idxFrame = 0;
                            for (std::map<unsigned int, unsigned int>::iterator it = visibility.begin(); it != visibility.end(); ++it, ++idxFrame ){
                                if(it->second  != -1){
                                    ceresObserv v;
                                    ++m_observations;
                                   int idx = it->second;
                                    v.oPt  = Point2Df(framesToAdjust[idxFrame]->getKeypoints()[idx]->getX() - K(0,2),
                                                      framesToAdjust[idxFrame]->getKeypoints()[idx]->getY() - K(1,2));
                                    v.cIdx = idxFrame;
                                    v.pIdx = i;

                                    observations_temp.push_back(v);
                                }
                             }
                        }
                    }
                    else{
                        std::cout<<" local bundler: { ";
                        for(int ii = 0; ii < selectedKeyframes.size(); ++ii)
                            std::cout<<selectedKeyframes[ii]<<" ";
                        std::cout<<"}"<<std::endl;

                        for(int i = 0; i < mapToAdjust.size(); ++i){
                            std::map<unsigned int, unsigned int> visibility = mapToAdjust[i]->getVisibility();
                            int idxFrame = 0;
                            for (std::map<unsigned int, unsigned int>::iterator it = visibility.begin(); it != visibility.end(); ++it, ++idxFrame){
                                for(int c = 0; c < selectedKeyframes.size(); ++c){
                                    if(idxFrame == selectedKeyframes[c] && it->second != -1){
                                        ceresObserv v;
                                        ++m_observations;
                                         int idx = it->second;
                                        v.oPt  = Point2Df(framesToAdjust[idxFrame]->getKeypoints()[idx]->getX() - K(0,2),
                                                          framesToAdjust[idxFrame]->getKeypoints()[idx]->getY() - K(1,2));
                                        v.cIdx = idxFrame;
                                        v.pIdx = i;
                                        observations_temp.push_back(v);
                                    }
                                }
                             }
                        }
                    }

                    m_observationsNo = observations_temp.size();
                    m_camerasNo = framesToAdjust.size();
                    m_pointsNo = mapToAdjust.size();

                    std::cout<<"    1.1->obs: "<<m_observationsNo<<" cams: "<<m_camerasNo<<std::endl;
                    m_observations = new double[OBSERV_DIM * m_observationsNo];
                    m_pointIndex   = new int[m_observationsNo];
                    m_cameraIndex  = new int[m_observationsNo];

                    m_parametersNo = CAM_DIM * m_camerasNo + POINT_DIM * m_pointsNo;
                    m_parameters = new double[m_parametersNo];

                    std::cout<<"    1.2->filling observation: "<<std::endl;
                    for (int i = 0; i < m_observationsNo; ++i) {
                        m_cameraIndex[i] = observations_temp[i].cIdx;
                        m_pointIndex[i]  = observations_temp[i].pIdx;

                        m_observations[OBSERV_DIM*i + 0] = observations_temp[i].oPt.getX();
                        m_observations[OBSERV_DIM*i + 1] = observations_temp[i].oPt.getY();

                    }
                    std::cout<<"    1.3->filling parameters (cam): "<<std::endl;

                    for(int  i = 0; i < m_camerasNo; ++i){
                            Vector3f r;
                            Transform3Df pose_temp = framesToAdjust[m_cameraIndex[i]]->getPose();
                            toRodrigues(pose_temp, r);
                            m_parameters[CAM_DIM*i+0] = r[0] ;
                            m_parameters[CAM_DIM*i+1] = r[1] ;
                            m_parameters[CAM_DIM*i+2] = r[2] ;

                            m_parameters[CAM_DIM*i+3] = pose_temp(0,3) ;
                            m_parameters[CAM_DIM*i+4] = pose_temp(1,3) ;
                            m_parameters[CAM_DIM*i+5] = pose_temp(2,3) ;

                            m_parameters[CAM_DIM*i+6] = K(0,0) ;
                            m_parameters[CAM_DIM*i+7] = D(0) ;
                            m_parameters[CAM_DIM*i+8] = D(1) ;

                    }
                    std::cout<<"    1.4->filling parameters (points): "<<std::endl;

                    for(int i = 0; i < mapToAdjust.size(); ++i){
                        m_parameters[POINT_DIM*i + m_camerasNo*CAM_DIM + 0] = mapToAdjust[i]->getX();
                        m_parameters[POINT_DIM*i + m_camerasNo*CAM_DIM + 1] = mapToAdjust[i]->getY();
                        m_parameters[POINT_DIM*i + m_camerasNo*CAM_DIM + 2] = mapToAdjust[i]->getZ();

                    }

                }
                bool SolARBundlerCeres::solveCeresProblem(){
                    for (int i = 0; i < num_observations(); ++i) {
                        ceres::CostFunction* cost_function =
                            SnavelyReprojectionError::Create(m_observations[OBSERV_DIM * i + 0],
                                                             m_observations[OBSERV_DIM * i + 1]);


                        m_problem.AddResidualBlock(cost_function,
                                                 NULL ,
                                                 mutable_camera_for_observation(i),
                                                 mutable_point_for_observation(i));
                    }





                    ceres::Solve(m_options, &m_problem, &m_summary);
                    std::cout << m_summary.FullReport() << "\n";

                    return true;
                }


                bool SolARBundlerCeres::updateMap(std::vector<SRef<CloudPoint>>&mapToAdjust){
                    for (int j = 0; j < get_points(); ++j) {
                        double x = m_parameters[(j * 3 + 0) + (CAM_DIM * get_cameras())];
                        double y = m_parameters[(j * 3 + 1) + (CAM_DIM * get_cameras())];
                        double z = m_parameters[(j * 3 + 2) + (CAM_DIM * get_cameras())];

                        double reprj_err = mapToAdjust[j]->getReprojError();
                        std::map<unsigned int, unsigned int>visibility = mapToAdjust[j]->getVisibility();
                        mapToAdjust[j] = xpcf::utils::make_shared<CloudPoint>(x, y, z,0.0,0.0,0.0,reprj_err,visibility);
                    }

                    return true;
                }
                bool SolARBundlerCeres::updateExtrinsic(std::vector<SRef<Keyframe>>&framesToAdjust,
                                                        const std::vector<int>&selectedKeyframes){
                    for (int j = 0; j < selectedKeyframes.size(); ++j) {
                        int idx =  selectedKeyframes[j];
                        Vector3d r,t, f;
                        Transform3Df pose_temp;
                        r[0] = m_parameters[CAM_DIM * j + 0];
                        r[1] = m_parameters[CAM_DIM * j + 1];
                        r[2] = m_parameters[CAM_DIM * j + 2];

                        iRodrigues(r,pose_temp);

                        pose_temp(0,3) = m_parameters[CAM_DIM * j + 3];
                        pose_temp(1,3) = m_parameters[CAM_DIM * j + 4];
                        pose_temp(2,3) = m_parameters[CAM_DIM * j + 5];

                        f[0] = m_parameters[CAM_DIM * j + 6];
                        f[1] = m_parameters[CAM_DIM * j + 7];
                        f[2] = m_parameters[CAM_DIM * j + 8];

                        framesToAdjust[idx]->setPose(pose_temp);
                    }
                  return true;
                }
                bool SolARBundlerCeres::updateIntrinsic(std::vector<SRef<Keyframe>>&framesToAdjust){
                    // do some stuff !
                    return true;
                }
                bool SolARBundlerCeres::updateCeresProblem(std::vector<SRef<Keyframe>>&framesToAdjust,
                                                           std::vector<SRef<CloudPoint>>&mapToAdjust,
                                                           const std::vector<int>&selectedKeyframes){
                    updateMap(mapToAdjust);
                    updateExtrinsic(framesToAdjust,selectedKeyframes);
                    updateIntrinsic(framesToAdjust);
                    return true;
                }
            }
       }
}
