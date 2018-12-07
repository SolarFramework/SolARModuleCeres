/**
 * @copyright Copyright (c) 2017 B-com http://www.b-com.com/
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#include <iostream>
#include <string>
#include <vector>


//#include <boost/log/core.hpp>

#include "SolARModuleOpengl_traits.h"
#include "SolARModuleCeres_traits.h"
#include "xpcf/xpcf.h"
#include "api/display/I3DPointsViewer.h"
#include "api/solver/map/IBundler.h"

using namespace SolAR;
using namespace SolAR::datastructure;
using namespace SolAR::api;
using namespace SolAR::MODULES::OPENGL;
using namespace SolAR::MODULES::CERES;
namespace xpcf = org::bcom::xpcf;


struct SolARBALoader{
    std::vector<std::vector<SRef<Keypoint>>>m_measurements;
	std::vector<Transform3Df>m_poses;
    std::vector<SRef<CloudPoint>>m_observations;
    std::vector<SRef<Image>> m_views;
    std::vector<SRef<DescriptorBuffer>> m_descriptors;
	CamCalibration  m_intrinsic;
	CamDistortion   m_distorsion;

    bool loadMeasurements(const std::string & path_measures) {
        int N;
        std::ifstream ox(path_measures);
        if (!ox.is_open()) {
            std::cerr << " can't read measurmeents file from: " << path_measures << std::endl;
            return false;
        }
        else {
            std::cout<<" loading measurement: ";
            ox >> N;
            m_measurements.resize(N);
            for (int i = 0; i < N; ++i) {
                std::cout<<i<<" ";
                std::string path_measure;
                ox >> path_measure;
                std::ifstream ox(path_measure);
                if (!ox.is_open()) {
                    std::cerr << " can't find observation file from: " << path_measure << std::endl;
                    return false;
                }
                else {
                    int kp_no;
                    ox >> kp_no;
                    m_measurements[i].resize(kp_no);
                    for (int j = 0; j < kp_no; ++j) {
                        float x,y;
                        ox >>x;
                        ox >>y;
                        m_measurements[i][j] = xpcf::utils::make_shared<Keypoint>(x,y,0.0,0.0,0.0,0.0,0);
                    }
                }   
            }
            std::cout<<" done"<<std::endl;
            return true;
        }
    }
    bool loadObservations(const std::string & path_obs) {
        std::ifstream ox(path_obs);
        if (!ox.is_open()) {
            std::cerr << "can't find cloud from: " << path_obs << std::endl;
            return false;
        }
        else{
            std::cout<<" loading observations: ";
            int obs_no;
            ox >> obs_no;
            m_observations.resize(obs_no);    
            for (int i = 0; i < obs_no; ++i) {
                double x,y,z;
                ox >> x;
                ox >> y;
                ox >> z;

                std::map<unsigned int, unsigned int> visibility_temp;

               m_observations[i] = xpcf::utils::make_shared<CloudPoint>(x, y, z,0.0,0.0,0.0,0.0,visibility_temp);
               int viz_no; ox >> viz_no;
               for(int j = 0; j < viz_no; ++j) {
                   int idxView,idxLoc;
                    ox >>idxView;
                    ox >>idxLoc;
                    m_observations[i]->getVisibility()[idxView] = idxLoc;
                }
            }
            std::cout<<" done"<<std::endl;
        }
        return true;
        
    }
    
    bool loadIntrinsic(const std::string&path_calib){
        std::cout<<"loading intrinsics: ";
        std::ifstream ox(path_calib);
        if(ox.is_open()){
           for(int i = 0; i < 3; ++i){
                for(int j = 0; j < 3; ++j){
                    ox>>m_intrinsic(i,j);
                }
            }
           std::cout<<"done"<<std::endl;
        }
        else{
            LOG_INFO("can't read calirabtion file from", path_calib);
           return false;
        }
        return true;
    }
    bool loadDistorsions(const std::string&path_dist) {
        std::ifstream ox(path_dist);
        if (!ox.is_open()) {
            LOG_INFO("can't read distorsion file from", path_dist);
            return false;
        }
        else {
            std::cout<<"loading intrinsic: ";
                for (int i = 0; i < 5; ++i) {
                    ox >> m_distorsion[i];
            }
                std::cout<<" done"<<std::endl;
        }
        return true;
    }

    bool loadPoses(const std::string & path_poses) {
        std::ifstream ox(path_poses);
        if (!ox.is_open()) {
            std::cerr << "can't find poses file from: " << path_poses << std::endl;
            return false;
        }
        else{
            std::cout<<" loading poses: ";
            int N;
            ox >> N;
            m_poses.resize(N);
            for (unsigned int i = 0; i < N; ++i) {
                for (int ii = 0; ii < 3; ++ii) {
                    for (int jj = 0; jj < 4; ++jj) {
                        ox >> m_poses[i](ii, jj);
                    }
                }
                m_poses[i](3,0) = 0.0;m_poses[i](3,1) = 0.0;m_poses[i](3,2) = 0.0;m_poses[i](3,3) = 1.0;
            }
        }
        std::cout<<" done"<<std::endl;
        return true;
    }
    
    void showPoses()const {
        int idx = 0;
        for (const auto &p : m_poses) {
            std::cout << " pose: " << idx << std::endl;
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 4; ++j) {
                    std::cout << p(i, j) << " ";
                }
                std::cout << std::endl;
            }
            ++idx;
            std::cout << std::endl;
        }
    }
    
    void showObservations()const {
        int idx = 0;
        std::cout << "<cloud>: " << std::endl;
        std::cout<<"    ->size: "<<m_observations.size()<<std::endl;
        for (unsigned int i = 0; i < m_observations.size(); ++i){
            std::cout << "p: " << m_observations[i]->getX() << " " <<  m_observations[i]->getY() << " " <<  m_observations[i]->getZ() << "  ";

            std::map<unsigned int, unsigned int> visibility = m_observations[i]->getVisibility();
            int idxFrame = 0;
            for (std::map<unsigned int, unsigned int>::iterator it = visibility.begin(); it != visibility.end(); ++it){
                std::cout<<it->first<<" "<<it->second<<" ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    
    void showMeasurements()const {
        std::cout << "<Measurements>: " << std::endl;
        for (int i = 0; i < m_measurements.size(); ++i) {
            std::cout << "	<Measurement: " << i << ">:" << std::endl;
            for (int j = 0; j < m_measurements[i].size(); ++j) {
                std::cout << m_measurements[i][j]->getX() << " " <<  m_measurements[i][j]->getY() << std::endl;
            }
        }
    }
    
    void  showIntrinsics()const {
        std::cout << "<Intrinsics>: " << std::endl;
            for (int ii = 0; ii < 3; ++ii) {
                for (int jj = 0; jj < 3; ++jj) {
                    std::cout << m_intrinsic(ii, jj) << " ";
                }
                std::cout << std::endl;
            }
    }
    
    void  showDistorsions()const {
        std::cout << "<distorsions>: " << std::endl;
            for (int ii = 0; ii < 5; ++ii) {
                    std::cout << m_distorsion[ii] << " ";
            }
            std::cout << std::endl;
    }
    
};

int run_bundle(std::string & scene){
    LOG_ADD_LOG_TO_CONSOLE();
    SolARBALoader *ba = new SolARBALoader();

//    std::string scene = "room6";
    
    const std::string path_poses        = "../" + scene + "Bundle/" + scene + "Poses.txt";
    const std::string path_obsrvations  = "../" + scene + "Bundle/" + scene + "Observations.txt";;
    const std::string path_measurements = "../" + scene + "Bundle/" + scene + "Measurements.txt";
    const std::string path_calibration  = "../" + scene + "Bundle/" + scene + "Calibration.txt";
    const std::string path_distorison   = "../" + scene + "Bundle/" + scene + "Distorsion.txt";

    /// Loading ba problem
    ba->loadObservations(path_obsrvations);
    ba->loadMeasurements(path_measurements);
    ba->loadPoses(path_poses);
    ba->loadIntrinsic(path_calibration);
    ba->loadDistorsions(path_distorison);
    /// showing ba problem
//    ba->showObservations();
//    ba->showMeasurements();
    ba->showPoses();
//    ba->showIntrinsics();



    LOG_INFO("-<SolARBALOADER LOADING>-");
    SRef<xpcf::IComponentManager> xpcfComponentManager = xpcf::getComponentManagerInstance();
    if(xpcfComponentManager->load("bundle_config.xml")!=org::bcom::xpcf::_SUCCESS)
    {
        LOG_ERROR("Failed to load the configuration file bundle_config.xml")
        return -1;
    }

    LOG_INFO("-<SolARBundlerCeres LOADING>-");
    auto bundler =xpcfComponentManager->create<SolARBundlerCeres>()->bindTo<api::solver::map::IBundler>();
    LOG_INFO("-<SolAR3DPointsViewerOpengl LOADING>-");
    auto viewer3DPoints =xpcfComponentManager->create<SolAR3DPointsViewerOpengl>()->bindTo<display::I3DPointsViewer>();

    ba->showDistorsions();

    std::vector<SRef<Keyframe>>keyframes;
    keyframes.resize(ba->m_poses.size());
    ba->m_descriptors.resize(ba->m_poses.size());
    ba->m_views.resize(ba->m_poses.size());

    for(unsigned int i = 0; i < keyframes.size(); ++i){
       keyframes[i] = xpcf::utils::make_shared<Keyframe>(ba->m_measurements[i],
                                                         ba->m_descriptors[i],
                                                         ba->m_views[i],
                                                         ba->m_poses[i]);


    }


    std::vector<int>selectedKeyframes;// = {0,1};
    std::vector<SRef<CloudPoint>>cloud, cloud_ba;
    std::vector<Transform3Df>poses , poses_ba;
    for(unsigned i = 0; i < keyframes.size(); ++i){
        Transform3Df tt = keyframes[i]->getPose();
        poses.push_back(tt);
    }

    cloud = ba->m_observations;

    bundler->adjustBundle(keyframes,
                          ba->m_observations,
                          ba->m_intrinsic,
                          ba->m_distorsion,
                          selectedKeyframes);

    for(auto & p: keyframes){
        Transform3Df tt = p->getPose();
        poses_ba.push_back(tt);
    }

    cloud_ba = ba->m_observations;

    std::vector<Transform3Df>framePoses;
    Transform3Df pp = Transform3Df::Identity();
    while(true){
        if (viewer3DPoints->display(cloud,
                                    pp,
                                    poses,
                                    framePoses,
                                    cloud_ba,
                                    poses_ba) == FrameworkReturnCode::_STOP){
                return 0;
            }
        }


    delete ba;
    return 0;
}
int main(int argc, char ** argv){
    if(argc != 2){
        std::cerr<<" ERROR! number of arguments is incorrect, exit.."<<std::endl;
        return -1;
    }
    std::string scene_name = argv[1];
    run_bundle(scene_name);

  return 0;
}



