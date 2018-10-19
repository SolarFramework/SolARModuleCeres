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
    std::vector<SRef<Keyframe>>m_keyframes;
    std::vector<std::vector<SRef<Keypoint>>>m_kps;
    std::vector<Transform3Df> m_poses;
    std::vector<SRef<Image>> m_views; // empty no need for bundle problem
    std::vector<SRef<DescriptorBuffer>> m_descriptors;  // empty no need for bundle problem
    CamCalibration m_intrinsics;
    CamDistortion m_distorsion;
    std::vector<SRef<CloudPoint>>m_cloud;
    const int m_keyframesNo  = 2;

    bool parse_keypoints(std::string&path, std::vector<SRef<Keypoint>>&kp, bool verbose){
        std::cout<<" loading keypoints: "<<std::endl;
        int kp_no;
        std::ifstream ox(path);
        if(ox.is_open()){
            ox>>kp_no;
            kp.resize(kp_no);
            for(int i = 0; i < kp_no; ++i){
              float x; ox>>x;
              float y; ox>>y;
              kp[i] = xpcf::utils::make_shared<Keypoint>(x,y,0.0,0.0,0.0,0.0,0);
              if(verbose){
                  std::cout<<"x: "<<x<<" y: "<<y<<std::endl;
              }
            }
        }
        else{
            LOG_INFO("can't read kp file from",path);
           return false;
        }
        std::cout<<"done"<<std::endl;
        return true;
    }
    bool parse_extrinsic(std::string&path, Transform3Df& pose, bool verbose){
        std::cout<<" loading pose: "<<std::endl;
        std::ifstream ox(path);
        if(ox.is_open()){
           for(int i = 0; i < 4; ++i){
                for(int j = 0; j < 4; ++j){
                    ox>>pose(i,j);
                    if(verbose)
                        std::cout<<pose(i,j)<<" ";
                }
                if(verbose)
                    std::cout<<std::endl;
            }
        }
        else{
            LOG_INFO("can't read pose file from",path);
           return false;
        }
        std::cout<<"done"<<std::endl;
        return true;
    }

    bool parse_intrinsic(std::string&path, CamCalibration& K, bool verbose){
        std::cout<<"loading intrinsics: "<<std::endl;
        std::ifstream ox(path);
        if(ox.is_open()){
           for(int i = 0; i < 3; ++i){
                for(int j = 0; j < 3; ++j){
                    ox>>K(i,j);
                    if(verbose)
                        std::cout<<K(i,j)<<" ";
                }
                if(verbose)
                    std::cout<<std::endl;
            }
        }
        else{
            LOG_INFO("can't read K file from",path);
           return false;
        }
        std::cout<<"done"<<std::endl;
        return true;
    }
    bool parse_distorsion(std::string&path, CamDistortion& D, bool verbose){
        std::cout<<"loading distorsion: "<<std::endl;
        std::ifstream ox(path);
        if(ox.is_open()){
           for(int i = 0; i < 5; ++i){
                ox>>D(i);
                if(verbose)
                    std::cout<<D(i)<<" ";
            }
        }
        else{
            LOG_INFO("can't read D file from",path);
           return false;
        }
        std::cout<<"done"<<std::endl;
        return true;
    }
    bool parse_map(std::string &path, std::vector<SRef<CloudPoint>>&cloud, bool verbose){
        std::cout<<"loading cloud: "<<std::endl;
        int points_no;
        std::ifstream ox(path);
        if(ox.is_open()){
            ox>>points_no;
            cloud.resize(points_no);
            for(int i = 0; i < points_no; ++i){
              double x; ox>>x;
              double y; ox>>y;
              double z; ox>>z;
              std::map<unsigned int, unsigned int> visibility;
              for(int v = 0; v < 2; ++v){
                  int id; ox>>id;
                  visibility[v] = id;
              }
              cloud[i] = xpcf::utils::make_shared<CloudPoint>(x,y,z,0.0,0.0,0.0,0.0, visibility);
              if(verbose){
                  std::cout<<"x: "<<cloud[i]->getX()<<" y: "<<cloud[i]->getY()<<" z: "<<cloud[i]->getZ()<<std::endl;
                  for (std::map<unsigned int, unsigned int>::iterator it = cloud[i]->getVisibility().begin(); it !=cloud[i]->getVisibility().end(); ++it){
                      std::cout<<it->second<<" ";
                  }
                  std::cout<<std::endl;
              }
            }
        }
        else{
            LOG_INFO("can't read kp file from",path);
           return false;
        }
        std::cout<<"done"<<std::endl;
        return true;
    }
    bool load(std::vector<std::string>&paths){
        if(paths.size()!=7){
            std::cerr<<" can't init ba problem, paths must be equal to 7"<<std::endl;
            return false;
        }else{
            m_poses.resize(m_keyframesNo); m_kps.resize(m_keyframesNo); m_keyframes.resize(m_keyframesNo);
            m_descriptors.resize(m_keyframesNo); m_views.resize(m_keyframesNo);
            parse_keypoints(paths[0], m_kps[0], false);
            parse_keypoints(paths[1], m_kps[1], false);

            parse_extrinsic(paths[2], m_poses[0], false);
            parse_extrinsic(paths[3], m_poses[1], false);

            parse_intrinsic(paths[4], m_intrinsics, false);
            parse_distorsion(paths[5], m_distorsion, false);
            parse_map(paths[6], m_cloud, false);

            // fill keyframes (warning   if m_keyframeNo > 2, update correctly datas).
            for(unsigned int c = 0; c < m_keyframesNo; ++c){
                m_keyframes[c] = xpcf::utils::make_shared<Keyframe>(m_kps[c],
                                                                    m_descriptors[c],
                                                                    m_views[c],
                                                                    m_poses[c]);
            }
        }
        return true;
    }
};

int run_bundle(){
    LOG_ADD_LOG_TO_CONSOLE();
    SolARBALoader *ba = new SolARBALoader();
    std::vector<std::string> paths; paths.resize(7);

    paths[0] = "kp0.txt";
    paths[1] = "kp1.txt";
    paths[2] = "pose0.txt";
    paths[3] = "pose1.txt";
    paths[4] = "K.txt";
    paths[5] =  "D.txt";
    paths[6] = "cloud.txt";


    LOG_INFO("-<SolARBALOADER LOADING>-");
    ba->load(paths);
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

    std::vector<int>selectedKeyframes = {0,1};
    std::vector<SRef<CloudPoint>>cloud_before = ba->m_cloud;
    std::vector<Transform3Df>poses_before;
    poses_before.resize(2);
    poses_before[0] = ba->m_keyframes[0]->getPose();
    poses_before[1] = ba->m_keyframes[1]->getPose();

    bundler->adjustBundle(ba->m_keyframes,
                          ba->m_cloud,
                          ba->m_intrinsics,
                          ba->m_distorsion,
                          selectedKeyframes);

    std::vector<SRef<CloudPoint>>cloud_after = ba->m_cloud;
    std::vector<Transform3Df>poses_after;
    poses_after.resize(2);
    poses_after[0] = ba->m_keyframes[0]->getPose();
    poses_after[1] = ba->m_keyframes[1]->getPose();


    std::vector<Transform3Df>framePoses;
    Transform3Df pp = Transform3Df::Identity();
    while(true){
        if (viewer3DPoints->display(cloud_before,
                                    pp,
                                    poses_before,
                                    framePoses,
                                    cloud_after,
                                    poses_after) == FrameworkReturnCode::_STOP){
                return 0;
            }
        }
    delete ba;
    return 0;
}
int main(){
    run_bundle();
  return 0;
}



