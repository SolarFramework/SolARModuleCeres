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

//#define USE_FREE


#include <iostream>
#include <string>
#include <vector>

#include <boost/log/core.hpp>
// ADD COMPONENTS HEADERS HERE

#include "SolARModuleOpencv_traits.h"
#include "SolARModuleOpengl_traits.h"
#include "SolARModuleNonFreeOpencv_traits.h"
#include "SolARModuleCeres_traits.h"
#include "SolARModuleTools_traits.h"

#include "SolARModuleSSBA_traits.h"


#include "xpcf/xpcf.h"

#include "api/input/devices/ICamera.h"
#include "api/features/IKeypointDetector.h"
#include "api/features/IDescriptorsExtractor.h"
#include "api/features/IDescriptorMatcher.h"
#include "api/solver/pose/I3DTransformFinderFrom2D2D.h"
#include "api/solver/map/ITriangulator.h"
#include "api/solver/map/IMapper.h"
#include "api/solver/map/IMapFilter.h"
#include "api/solver/pose/I2D3DCorrespondencesFinder.h"
#include "api/solver/pose/I3DTransformFinderFrom2D3D.h"
#include "api/features/IMatchesFilter.h"
#include "api/display/ISideBySideOverlay.h"
#include "api/display/I2DOverlay.h"
#include "api/display/I3DOverlay.h"
#include "api/display/IImageViewer.h"
#include "api/display/I3DPointsViewer.h"

#include "api/image/IImageLoader.h"
#include "api/solver/map/IBundler.h"

#include "SolAROpenCVHelper.h"
#include "opencv2/highgui.hpp"
#include "opencv2/core.hpp"


using namespace SolAR;
using namespace SolAR::datastructure;
using namespace SolAR::api;
using namespace SolAR::MODULES::OPENCV;
using namespace SolAR::MODULES::NONFREEOPENCV;
using namespace SolAR::MODULES::OPENGL;
using namespace SolAR::MODULES::CERES;
using namespace SolAR::MODULES::TOOLS;
using namespace  SolAR::MODULES::SSBA;

namespace xpcf = org::bcom::xpcf;

struct streamConfig{

  std::string m_dir;
  int m_viewsNo;

  void show(){
      std::cout<<"  ->stream configuration: "<<std::endl;
      std::cout<<"      # views directory: "<<m_dir<<std::endl;
            std::cout<<"      # views number: "<<m_viewsNo<<std::endl;
  }

  bool load(std::string&path){
      std::ifstream ox(path);
      if(!ox.is_open()){
          std::cerr<<" can't read stream config file from: "<<path<<std::endl;
          return false;
      }
      else{
          std::string v[2];
          for(int i = 0; i < 2; ++i)
              ox>>v[i];
          m_dir = v[0];
          m_viewsNo = std::stoi(v[1]);
          show();

          return true;
      }
  }
  streamConfig(){
      m_dir = ""; m_viewsNo = 0;
  }
  streamConfig(std::string a, int b): m_dir(a), m_viewsNo(b){}
};

void basic_mapFiltering( std::vector<SRef<CloudPoint>>&in, double thresh, double&reproj_err){
     std::vector<SRef<CloudPoint>> out;
     out.reserve(in.size());
    reproj_err = 0.0;
    for(const auto &p: in){
        if(p->getReprojError()<thresh){
            out.push_back(p);
            reproj_err+= p->getReprojError();
        }
    }
    in = out;
    reproj_err/=(double)out.size();
}


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
              cv::waitKey(0);
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
bool parse_pose(std::string&path, Transform3Df& pose, bool verbose){
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
       if(verbose)
            cv::waitKey(0);
    }
    else{
        LOG_INFO("can't read pose file from",path);
       return false;
    }
    std::cout<<"done"<<std::endl;
    return true;
}
bool parse_K(std::string&path, CamCalibration& K, bool verbose){
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
       if(verbose)
            cv::waitKey(0);
    }
    else{
        LOG_INFO("can't read K file from",path);
       return false;
    }
    std::cout<<"done"<<std::endl;
    return true;
}
bool parse_D(std::string&path, CamDistortion& D, bool verbose){
    std::cout<<"loading distorsion: "<<std::endl;
    std::ifstream ox(path);
    if(ox.is_open()){
       for(int i = 0; i < 5; ++i){
            ox>>D(i);
            if(verbose)
                std::cout<<D(i)<<" ";
        }
       if(verbose)
        cv::waitKey(0);
    }
    else{
        LOG_INFO("can't read D file from",path);
       return false;
    }
    std::cout<<"done"<<std::endl;
    return true;
}

bool parse_cloud(std::string &path, std::vector<SRef<CloudPoint>>&cloud, bool verbose){
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
              cv::waitKey(0);
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

int run_bundleFromTxt(){
    cv::namedWindow("debug txt");
    std::string path_kp0 = "D:/kp0.txt";
    std::string path_kp1 = "D:/kp1.txt";
    std::string path_pose0 = "D:/pose0.txt";
    std::string path_pose1 = "D:/pose1.txt";
    std::string path_K = "D:/K.txt";
    std::string path_D =  "D:/D.txt";
    std::string path_cloud = "D:/cloud.txt";

    LOG_ADD_LOG_TO_CONSOLE();
    SRef<xpcf::IComponentManager> xpcfComponentManager = xpcf::getComponentManagerInstance();
    if(xpcfComponentManager->load("bundle_config.xml")!=org::bcom::xpcf::_SUCCESS)
    {
        LOG_ERROR("Failed to load the configuration file bundle_config.xml")
        return -1;
    }

    LOG_INFO("SolARBundlerCeres LOADING-------------------------------------------------: ");
    auto bundler =xpcfComponentManager->create<SolARBundlerCeres>()->bindTo<api::solver::map::IBundler>();
    LOG_INFO("SolAR3DPointsViewerOpengl LOADING-------------------------------------------------: ");
    auto viewer3DPoints =xpcfComponentManager->create<SolAR3DPointsViewerOpengl>()->bindTo<display::I3DPointsViewer>();


    std::vector<SRef<Keyframe>>keyframes;        keyframes.resize(2);
    std::vector<std::vector<SRef<Keypoint>>>kps; kps.resize(2);
    std::vector<Transform3Df> poses;    poses.resize(2);
    std::vector<SRef<Image>> views; views.resize(2);
    std::vector<SRef<DescriptorBuffer>> descriptors; descriptors.resize(2);
    CamCalibration K;
    CamDistortion D;
    std::vector<SRef<CloudPoint>>cloud;

    parse_cloud(path_cloud, cloud, false);
    parse_K(path_K, K, false);
    parse_D(path_D, D, false);

    parse_pose(path_pose0, poses[0], false);
    parse_pose(path_pose1, poses[1], false);

    parse_keypoints(path_kp0, kps[0], false);
    parse_keypoints(path_kp1, kps[1], false);

    keyframes[0] = xpcf::utils::make_shared<Keyframe>(kps[0],
                                                     descriptors[0],
                                                     views[0],
                                                     poses[0]);


    keyframes[1] = xpcf::utils::make_shared<Keyframe>(kps[1],
                                                     descriptors[1],
                                                     views[1],
                                                     poses[1]);


    std::vector<int>selectedKeyframes = {0,1};



    std::vector<SRef<CloudPoint>>cloud_before = cloud;
    bundler->adjustBundle(keyframes,
                          cloud,
                          K,
                          D,
                          selectedKeyframes);



    std::vector<SRef<CloudPoint>>cloud_after;

    cloud_after = cloud;
    Transform3Df p0,p1;

    p0 = keyframes[0]->getPose();
    p1 = keyframes[1]->getPose();

    std::vector<Transform3Df>KeyframePoses_after;
    KeyframePoses_after.push_back(p0);
    KeyframePoses_after.push_back(p1);
    std::cout<<" after: "<<cloud_after.size()<<std::endl;

    std::vector<Transform3Df>framePoses;
    Transform3Df pp = Transform3Df::Identity();
    while(true){
        if (viewer3DPoints->display(cloud_before,
                                    pp,
                                    poses,
                                    framePoses,
                                    cloud_after,
                                    KeyframePoses_after) == FrameworkReturnCode::_STOP){
                return 0;
            }
        }
    return 0;
}

int run_onlybundle(){
    // stream config need data folder (ask for it).
    std::string path_stream = "stream_config.txt";
    streamConfig* mStream = new streamConfig();
    mStream->load(path_stream);

    LOG_ADD_LOG_TO_CONSOLE();
    SRef<xpcf::IComponentManager> xpcfComponentManager = xpcf::getComponentManagerInstance();
    if(xpcfComponentManager->load("bundle_config.xml")!=org::bcom::xpcf::_SUCCESS)
    {
        LOG_ERROR("Failed to load the configuration file bundle_config.xml")
        return -1;
    }

    LOG_INFO("Start creating components");
    LOG_INFO("CAMERA LOADING-------------------------------------------------: ");
    auto camera =xpcfComponentManager->create<SolARCameraOpencv>()->bindTo<input::devices::ICamera>();

#ifdef USE_FREE
    LOG_INFO("SolARKeypointDetectorOpencv LOADING-------------------------------------------------: ");
    auto keypointsDetector =xpcfComponentManager->create<SolARKeypointDetectorOpencv>()->bindTo<features::IKeypointDetector>();
    LOG_INFO("SolARDescriptorsExtractorAKAZE2Opencv LOADING-------------------------------------------------: ");
    auto descriptorExtractor =xpcfComponentManager->create<SolARDescriptorsExtractorAKAZE2Opencv>()->bindTo<features::IDescriptorsExtractor>();
#else


   auto  keypointsDetector = xpcfComponentManager->create<SolARKeypointDetectorNonFreeOpencv>()->bindTo<features::IKeypointDetector>();
   auto descriptorExtractor = xpcfComponentManager->create<SolARDescriptorsExtractorSIFTOpencv>()->bindTo<features::IDescriptorsExtractor>();

#endif

    LOG_INFO("SolARDescriptorMatcherKNNOpencv LOADING-------------------------------------------------: ");
    auto matcher =xpcfComponentManager->create<SolARDescriptorMatcherKNNOpencv>()->bindTo<features::IDescriptorMatcher>();
    LOG_INFO("SolARPoseFinderFrom2D2DOpencvD LOADING-------------------------------------------------: ");
    auto poseFinderFrom2D2D =xpcfComponentManager->create<SolARPoseFinderFrom2D2DOpencv>()->bindTo<solver::pose::I3DTransformFinderFrom2D2D>();
    LOG_INFO("SolARSVDTriangulationOpencv LOADING-------------------------------------------------: ");
    auto mapper =xpcfComponentManager->create<SolARSVDTriangulationOpencv>()->bindTo<solver::map::ITriangulator>();
    LOG_INFO("SolARMapper LOADING-------------------------------------------------: ");
    auto poseGraph =xpcfComponentManager->create<SolARMapper>()->bindTo<solver::map::IMapper>();
    LOG_INFO("SolAR3DPointsViewerOpengl LOADING-------------------------------------------------: ");
    auto viewer3DPoints =xpcfComponentManager->create<SolAR3DPointsViewerOpengl>()->bindTo<display::I3DPointsViewer>();
    LOG_INFO("SolARMapFilter LOADING-------------------------------------------------: ");
    auto  mapFilter =xpcfComponentManager->create<SolARMapFilter>()->bindTo<solver::map::IMapFilter>();
    LOG_INFO("SolARBundlerCeres LOADING-------------------------------------------------: ");
    auto bundler =xpcfComponentManager->create<SolARBundlerCeres>()->bindTo<api::solver::map::IBundler>();





    std::vector<SRef<Image>>                            views;
    SRef<Keyframe>                                      keyframe[2];
    std::vector<SRef<Keypoint>>                         keypoints[2];
    SRef<DescriptorBuffer>                              descriptors[2];
    std::vector<DescriptorMatch>                        matches;
    std::vector<SRef<CloudPoint>>                       cloud, filtredCloud;
    std::vector<Transform3Df>                           keyframePoses;





    SRef<Map> map;
    views.resize(mStream->m_viewsNo);
    cv::Mat cvView;
    LOG_INFO("loading views: ");
    int c = 0;
    for(unsigned int i = 0; i < mStream->m_viewsNo; ++i){
        std::string path_temp = mStream->m_dir + std::to_string(i) + ".png";
        cvView = cv::imread(path_temp);
        if(!cvView.empty()) ++c;
        SolAROpenCVHelper::convertToSolar(cvView, views[i]);
    }
    std::cout<<c <<"/"<<mStream->m_viewsNo<<" loaded correctly"<<std::endl;
    keyframePoses.resize(2);



    poseFinderFrom2D2D->setCameraParameters(camera->getIntrinsicsParameters(), camera->getDistorsionParameters());
    mapper->setCameraParameters(camera->getIntrinsicsParameters(), camera->getDistorsionParameters());

    LOG_DEBUG("Intrincic parameters : \n {}", camera->getIntrinsicsParameters());

    if (camera->start() != FrameworkReturnCode::_SUCCESS) // videoFile
    {
        LOG_ERROR("Camera cannot start");
        return -1;
    }

    for(unsigned int v = 0; v < 2; ++v){
        keypointsDetector->detect(views[v], keypoints[v]);
        descriptorExtractor->extract(views[v], keypoints[v], descriptors[v]);
    }
    matcher->match(descriptors[0], descriptors[1], matches);

    int nbOriginalMatches = matches.size();
    keyframePoses[0] = Transform3Df::Identity();

    // Estimate the pose of of the second frame (the first frame being the reference of our coordinate system)
    poseFinderFrom2D2D->estimate(keypoints[0], keypoints[1], keyframePoses[0], keyframePoses[1], matches);

    LOG_INFO("Nb matches for triangulation: {}\\{}", matches.size(), nbOriginalMatches);
    LOG_INFO("Estimate pose of the camera for the frame 2: \n {}", keyframePoses[1].matrix());

    // Triangulate

    keyframe[0] = xpcf::utils::make_shared<Keyframe>(keypoints[0],
                                                     descriptors[0],
                                                     views[0],
                                                     keyframePoses[0]);



    poseGraph->update(map, keyframe[0]);

    SRef<Frame> frame2 = xpcf::utils::make_shared<Frame>(keypoints[1], descriptors[1], views[1], keyframe[0]);

    frame2->setPose(keyframePoses[1]);

    double reproj_error = mapper->triangulate(keypoints[0],
                                              keypoints[1],
                                              matches,
                                              std::make_pair(0, 1),
                                              keyframePoses[0],
                                              keyframePoses[1],
                                              cloud);





    mapFilter->filter(keyframePoses[0], keyframePoses[1], cloud, filtredCloud);


    keyframe[1] = xpcf::utils::make_shared<Keyframe>(frame2);

    poseGraph->update(map,
                      keyframe[1],
                      filtredCloud,
                      matches);

    for(unsigned int c = 0; c < 2; ++c){
        std::ofstream logKP("D:/kp" + std::to_string(c) +".txt");
        logKP<<keypoints[c].size()<<std::endl;
        for(auto & kp: keypoints[c]){
            logKP<<kp->getX()<<" "<<kp->getY()<<std::endl;
        }
        logKP.close();
    }

    for(unsigned int c = 0; c < 2; ++c){
        std::ofstream logPose("D:/pose" + std::to_string(c) +".txt");
        for(int i = 0; i < 4; ++i){
            for(int j = 0; j < 4; ++j){
                logPose<<keyframePoses[c](i,j)<<" ";
            }
            logPose<<std::endl;
        }
        logPose.close();
    }

    std::ofstream logCloud("D:/cloud.txt");
    logCloud<<filtredCloud.size()<<std::endl;
    for(const auto & pp: filtredCloud){
         std::map<unsigned int, unsigned int> visibility = pp->getVisibility();
        logCloud<<pp->getX()<<" "<<pp->getY()<<" "<<pp->getZ()<< " ";
        for (std::map<unsigned int, unsigned int>::iterator it = visibility.begin(); it != visibility.end(); ++it){
            logCloud<<it->second<<" ";
        }
        logCloud<<std::endl;
    }
    logCloud<<std::endl;

    std::ofstream logK("D:/K.txt");
    for(unsigned int i =0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            logK<<camera->getIntrinsicsParameters()(i,j)<<" ";
        }
        logK<<std::endl;
    }
    logK<<std::endl;


    std::ofstream logD("D:/D.txt");
    for(unsigned int i =0; i < 5; ++i){
        logD<<camera->getDistorsionParameters()(i)<<" ";
    }
    logD<<std::endl;


    std::vector<int>selectedKeyframes = {0,1};

    std::vector<SRef<CloudPoint>>cloud_before_ba;
    std::vector<SRef<CloudPoint>>cloud_after_ba;
    std::vector<SRef<Keyframe>> kf=poseGraph->getKeyframes();
    cloud_before_ba = *poseGraph->getMap()->getPointCloud();

    bundler->adjustBundle(kf,
                          *poseGraph->getMap()->getPointCloud(),
                          camera->getIntrinsicsParameters(),
                          camera->getDistorsionParameters(),
                          selectedKeyframes);

    cloud_after_ba = *poseGraph->getMap()->getPointCloud();

    std::vector<Transform3Df>KeyframePoses_after;
    Transform3Df p0,p1;

    p0 = (poseGraph->getKeyframes()[0]->getPose());
    p1 = (poseGraph->getKeyframes()[1]->getPose());

    KeyframePoses_after.push_back(p0);
    KeyframePoses_after.push_back(p1);
    std::cout<<" after: "<<cloud_after_ba.size()<<std::endl;

    std::vector<Transform3Df>framePoses;
    Transform3Df pp = Transform3Df::Identity();
    while(true){
    if (viewer3DPoints->display(cloud_before_ba,
                                    pp,
                                    keyframePoses,
                                    framePoses,
                                    cloud_after_ba,
                                    KeyframePoses_after) == FrameworkReturnCode::_STOP){
            return 0;
        }
    }


    return 0;

}
int run_bundle(){

    // stream config need data folder (ask for it).
    std::string path_stream = "stream_config.txt";
    streamConfig* mStream = new streamConfig();
    mStream->load(path_stream);

    LOG_ADD_LOG_TO_CONSOLE();
    SRef<xpcf::IComponentManager> xpcfComponentManager = xpcf::getComponentManagerInstance();
    if(xpcfComponentManager->load("bundle_config.xml")!=org::bcom::xpcf::_SUCCESS)
    {
        LOG_ERROR("Failed to load the configuration file bundle_config.xml")
        return -1;
    }

    LOG_INFO("Start creating components");
    LOG_INFO("CAMERA LOADING-------------------------------------------------: ");
    auto camera =xpcfComponentManager->create<SolARCameraOpencv>()->bindTo<input::devices::ICamera>();

#ifdef USE_FREE
    LOG_INFO("SolARKeypointDetectorOpencv LOADING-------------------------------------------------: ");
    auto keypointsDetector =xpcfComponentManager->create<SolARKeypointDetectorOpencv>()->bindTo<features::IKeypointDetector>();
    LOG_INFO("SolARDescriptorsExtractorAKAZE2Opencv LOADING-------------------------------------------------: ");
    auto descriptorExtractor =xpcfComponentManager->create<SolARDescriptorsExtractorAKAZE2Opencv>()->bindTo<features::IDescriptorsExtractor>();
#else


   auto  keypointsDetector = xpcfComponentManager->create<SolARKeypointDetectorNonFreeOpencv>()->bindTo<features::IKeypointDetector>();
   auto descriptorExtractor = xpcfComponentManager->create<SolARDescriptorsExtractorSIFTOpencv>()->bindTo<features::IDescriptorsExtractor>();

#endif

    LOG_INFO("SolARDescriptorMatcherKNNOpencv LOADING-------------------------------------------------: ");
    auto matcher =xpcfComponentManager->create<SolARDescriptorMatcherKNNOpencv>()->bindTo<features::IDescriptorMatcher>();
    LOG_INFO("SolARPoseFinderFrom2D2DOpencvD LOADING-------------------------------------------------: ");
    auto poseFinderFrom2D2D =xpcfComponentManager->create<SolARPoseFinderFrom2D2DOpencv>()->bindTo<solver::pose::I3DTransformFinderFrom2D2D>();
    LOG_INFO("SolARSVDTriangulationOpencv LOADING-------------------------------------------------: ");
    auto mapper =xpcfComponentManager->create<SolARSVDTriangulationOpencv>()->bindTo<solver::map::ITriangulator>();
    LOG_INFO("SolARMapper LOADING-------------------------------------------------: ");
    auto poseGraph =xpcfComponentManager->create<SolARMapper>()->bindTo<solver::map::IMapper>();
    LOG_INFO("SolAR3DPointsViewerOpengl LOADING-------------------------------------------------: ");
    auto viewer3DPoints =xpcfComponentManager->create<SolAR3DPointsViewerOpengl>()->bindTo<display::I3DPointsViewer>();
    LOG_INFO("SolARMapFilter LOADING-------------------------------------------------: ");
    auto  mapFilter =xpcfComponentManager->create<SolARMapFilter>()->bindTo<solver::map::IMapFilter>();  
    LOG_INFO("SolARBundlerCeres LOADING-------------------------------------------------: ");
    auto bundler =xpcfComponentManager->create<SolARBundlerCeres>()->bindTo<api::solver::map::IBundler>();

    /*

    std::vector<SRef<CloudPoint>> mapToAdjust;
    std::vector<SRef<Keyframe>> framesToAdjust;

    CamDistortion D;
    CamCalibration K;
    std::vector<int>selectedKeyframes;
    std::cout<<" performing bundle.."<<std::endl;
    bundler->adjustBundle(framesToAdjust,
                          mapToAdjust,
                          K,
                          D,
                          selectedKeyframes);
                          */






    std::vector<SRef<Image>>                            views;
    SRef<Keyframe>                                      keyframe[2];
    std::vector<SRef<Keypoint>>                         keypoints[2];
    SRef<DescriptorBuffer>                              descriptors[2];
    std::vector<DescriptorMatch>                        matches;
    std::vector<SRef<CloudPoint>>                       cloud, filtredCloud;
    std::vector<Transform3Df>                           keyframePoses;





    SRef<Map> map;
    views.resize(mStream->m_viewsNo);
    cv::Mat cvView;
    LOG_INFO("loading views: ");
    int c = 0;
    for(unsigned int i = 0; i < mStream->m_viewsNo; ++i){
        std::string path_temp = mStream->m_dir + std::to_string(i) + ".png";
        cvView = cv::imread(path_temp);
        if(!cvView.empty()) ++c;
        SolAROpenCVHelper::convertToSolar(cvView, views[i]);
    }
    std::cout<<c <<"/"<<mStream->m_viewsNo<<" loaded correctly"<<std::endl;
    keyframePoses.resize(2);



    poseFinderFrom2D2D->setCameraParameters(camera->getIntrinsicsParameters(), camera->getDistorsionParameters());
    mapper->setCameraParameters(camera->getIntrinsicsParameters(), camera->getDistorsionParameters());

    LOG_DEBUG("Intrincic parameters : \n {}", camera->getIntrinsicsParameters());

    if (camera->start() != FrameworkReturnCode::_SUCCESS) // videoFile
    {
        LOG_ERROR("Camera cannot start");
        return -1;
    }

    for(unsigned int v = 0; v < 2; ++v){
        keypointsDetector->detect(views[v], keypoints[v]);
        descriptorExtractor->extract(views[v], keypoints[v], descriptors[v]);
    }
    matcher->match(descriptors[0], descriptors[1], matches);

    int nbOriginalMatches = matches.size();
    keyframePoses[0] = Transform3Df::Identity();

    // Estimate the pose of of the second frame (the first frame being the reference of our coordinate system)
    poseFinderFrom2D2D->estimate(keypoints[0], keypoints[1], keyframePoses[0], keyframePoses[1], matches);

    LOG_INFO("Nb matches for triangulation: {}\\{}", matches.size(), nbOriginalMatches);
    LOG_INFO("Estimate pose of the camera for the frame 2: \n {}", keyframePoses[1].matrix());

    // Triangulate

    keyframe[0] = xpcf::utils::make_shared<Keyframe>(keypoints[0],
                                                     descriptors[0],
                                                     views[0],
                                                     keyframePoses[0]);

    poseGraph->update(map, keyframe[0]);

    SRef<Frame> frame2 = xpcf::utils::make_shared<Frame>(keypoints[1], descriptors[1], views[1], keyframe[0]);

    frame2->setPose(keyframePoses[1]);

    double reproj_error = mapper->triangulate(keypoints[0],
                                              keypoints[1],
                                              matches,
                                              std::make_pair(0, 1),
                                              keyframePoses[0],
                                              keyframePoses[1],
                                              cloud);





    mapFilter->filter(keyframePoses[0], keyframePoses[1], cloud, filtredCloud);


    keyframe[1] = xpcf::utils::make_shared<Keyframe>(frame2);

    poseGraph->update(map,
                      keyframe[1],
                      filtredCloud,
                      matches);


    std::vector<int>selectedKeyframes = {0,1};

    std::vector<SRef<CloudPoint>>cloud_before_ba;
    std::vector<SRef<CloudPoint>>cloud_after_ba;
    std::vector<SRef<Keyframe>> kf=poseGraph->getKeyframes();
    cloud_before_ba = *poseGraph->getMap()->getPointCloud();

    bundler->adjustBundle(kf,
                          *poseGraph->getMap()->getPointCloud(),
                          camera->getIntrinsicsParameters(),
                          camera->getDistorsionParameters(),
                          selectedKeyframes);

    cloud_after_ba = *poseGraph->getMap()->getPointCloud();

    std::vector<Transform3Df>KeyframePoses_after;
    Transform3Df p0,p1;

    p0 = (poseGraph->getKeyframes()[0]->getPose());
    p1 = (poseGraph->getKeyframes()[1]->getPose());

    KeyframePoses_after.push_back(p0);
    KeyframePoses_after.push_back(p1);
    std::cout<<" after: "<<cloud_after_ba.size()<<std::endl;

    std::vector<Transform3Df>framePoses;
    Transform3Df pp = Transform3Df::Identity();
    while(true){
    if (viewer3DPoints->display(cloud_before_ba,
                                    pp,
                                    keyframePoses,
                                    framePoses,
                                    cloud_after_ba,
                                    KeyframePoses_after) == FrameworkReturnCode::_STOP){
            return 0;
        }
    }


    return 0;
}
int main(){
    cv::namedWindow("debug window", 0);
    run_bundleFromTxt();
 //   run_onlybundle();
//    run_bundle();
    cv::waitKey(0);
  return 0;
}



