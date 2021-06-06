// g2o - General Graph Optimization
// Copyright (C) 2011 H. Strasdat
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <iostream>
#include <stdint.h>

#include <unordered_set>

#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/solver.h"
#include "g2o/core/robust_kernel_impl.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/solvers/dense/linear_solver_dense.h"
#include "g2o/types/sba/types_six_dof_expmap.h"
#include "g2o/solvers/structure_only/structure_only_solver.h"
#include "g2o/stuff/sampler.h"

#if defined G2O_HAVE_CHOLMOD
#include "g2o/solvers/cholmod/linear_solver_cholmod.h"
#else
#include "g2o/solvers/eigen/linear_solver_eigen.h"
#endif

using namespace Eigen;
using namespace std;

class Sample {
 public:
  static int uniform(int from, int to) { return static_cast<int>(g2o::Sampler::uniformRand(from, to)); }
};

//BA的入口函数。这个demo是同时估计相机的SE(3)和路标点的3D坐标
//随机生成500个路标点true_points
int main(int argc, const char* argv[]){
  if (argc<2)
  {
    cout << endl;
    cout << "Please type: " << endl;
    cout << "ba_demo [PIXEL_NOISE] [OUTLIER RATIO] [ROBUST_KERNEL] [STRUCTURE_ONLY] [DENSE]" << endl;
    cout << endl;
    cout << "PIXEL_NOISE: noise in image space (E.g.: 1)" << endl;
    cout << "OUTLIER_RATIO: probability of spuroius observation  (default: 0.0)" << endl;
    cout << "ROBUST_KERNEL: use robust kernel (0 or 1; default: 0==false)" << endl;
    cout << "STRUCTURE_ONLY: performe structure-only BA to get better point initializations (0 or 1; default: 0==false)" << endl;
    cout << "DENSE: Use dense solver (0 or 1; default: 0==false)" << endl;
    cout << endl;
    cout << "Note, if OUTLIER_RATIO is above 0, ROBUST_KERNEL should be set to 1==true." << endl;
    cout << endl;
    exit(0);
  }

  double PIXEL_NOISE = atof(argv[1]);
  double OUTLIER_RATIO = 0.0;

  if (argc>2)  {
    OUTLIER_RATIO = atof(argv[2]);
  }

  bool ROBUST_KERNEL = false;
  if (argc>3){
    ROBUST_KERNEL = atoi(argv[3]) != 0;
  }
  bool STRUCTURE_ONLY = false;
  if (argc>4){
    STRUCTURE_ONLY = atoi(argv[4]) != 0;
  }

  bool DENSE = false;
  if (argc>5){
    DENSE = atoi(argv[5]) != 0;
  }

  cout << "PIXEL_NOISE: " <<  PIXEL_NOISE << endl;
  cout << "OUTLIER_RATIO: " << OUTLIER_RATIO<<  endl;
  cout << "ROBUST_KERNEL: " << ROBUST_KERNEL << endl;
  cout << "STRUCTURE_ONLY: " << STRUCTURE_ONLY<< endl;
  cout << "DENSE: "<<  DENSE << endl;
  
  //Step1:初始化优化器，设置线性求解器。LM的优化方法，线性求解器根据输入DENSE的值以及系统是否安装了CHOLMOD来决定。
  g2o::SparseOptimizer optimizer;
  optimizer.setVerbose(true);
  std::unique_ptr<g2o::BlockSolver_6_3::LinearSolverType> linearSolver; //这里的6,3代表了相机姿态是6自由度，3是路标点的姿态
  if (DENSE) { //如果是稠密，那么使用LinearSolverDense进行线性方程的求解
    linearSolver = g2o::make_unique<g2o::LinearSolverDense<g2o::BlockSolver_6_3::PoseMatrixType>>();
  } else {
#ifdef G2O_HAVE_CHOLMOD
    using BaLinearSolver = g2o::LinearSolverCholmod<g2o::BlockSolver_6_3::PoseMatrixType>; //否则，使用LinearSolverCholmod
#else
    using BaLinearSolver = g2o::LinearSolverEigen<g2o::BlockSolver_6_3::PoseMatrixType>;
#endif
    linearSolver = g2o::make_unique<BaLinearSolver>();
  }

  g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(
    g2o::make_unique<g2o::BlockSolver_6_3>(std::move(linearSolver))
  );//传入BlockSolver创建一个OptimizationAlgorithmLevenberg类型的OptimizationAlgorithm
  optimizer.setAlgorithm(solver);

  //Step2:生成500个路标点true_points和15个相机的位姿true_poses。
  //
  vector<Vector3d> true_points; //路标点的3D坐标。500个路标点
  for (size_t i=0;i<500; ++i)
  {
    true_points.push_back(Vector3d((g2o::Sampler::uniformRand(0., 1.)-0.5)*3,
                                   g2o::Sampler::uniformRand(0., 1.)-0.5,
                                   g2o::Sampler::uniformRand(0., 1.)+3));
  }

  double focal_length= 1000.; //焦距
  Vector2d principal_point(320., 240.); //成像平面光心的像素坐标

  vector<g2o::SE3Quat,
      aligned_allocator<g2o::SE3Quat> > true_poses; //相机的位姿
  //Step3:利用相机的参数生成一个相机参数的顶点CameraParameters
  g2o::CameraParameters * cam_params
      = new g2o::CameraParameters (focal_length, principal_point, 0.); //相机参数
  cam_params->setId(0);

  if (!optimizer.addParameter(cam_params)) {
    assert(false);
  }

  int vertex_id = 0; //给顶点编号用的
  //Step4:相机的真实位姿添加到优化器的顶点
  for (size_t i=0; i<15; ++i) { //共有15个顶点
    Vector3d trans(i*0.04-1.,0,0); //平移递增，旋转不变

    Eigen:: Quaterniond q;
    q.setIdentity();
    g2o::SE3Quat pose(q,trans);
    g2o::VertexSE3Expmap * v_se3
        = new g2o::VertexSE3Expmap();
    v_se3->setId(vertex_id);
    if (i<2){
      v_se3->setFixed(true);
    }
    v_se3->setEstimate(pose);
    optimizer.addVertex(v_se3);
    true_poses.push_back(pose);
    vertex_id++;
  }
  int point_id=vertex_id; //路标点的id是从相机id之后开始的
  int point_num = 0;
  double sum_diff2 = 0;

  cout << endl;
  unordered_map<int,int> pointid_2_trueid;
  unordered_set<int> inliers;
  //Step5:遍历每一个地图点，首先基于真实值添加一定的白噪声，产生带噪声的3D坐标点，然后，检查这个点能够被15个相机位姿
  //的多少个看到，保证至少两个位姿能够看到，符合这样要求的点被添加到VertexSBAPointXYZ类型的顶点。然后，为了模拟真实环境
  //的数据，根据OUTLIER_RATIO来产生一些外点，外点的坐标是随机游走在640*480的任何一个格子上。这一步也生成了边。
  for (size_t i=0; i<true_points.size(); ++i){//遍历每一个路标点
    g2o::VertexSBAPointXYZ * v_p
        = new g2o::VertexSBAPointXYZ();
    v_p->setId(point_id);
    v_p->setMarginalized(true); //边缘化路标点
    v_p->setEstimate(true_points.at(i)
                     + Vector3d(g2o::Sampler::gaussRand(0., 1),
                                g2o::Sampler::gaussRand(0., 1),
                                g2o::Sampler::gaussRand(0., 1))); //设定初始估计值，是真值加高斯噪声
    int num_obs = 0;
    for (size_t j=0; j<true_poses.size(); ++j){//遍历每一个相机顶点，看看当前路标点能不能被这个相机位姿看到，如果看到，就加一下num_obs
      Vector2d z = cam_params->cam_map(true_poses.at(j).map(true_points.at(i))); // CameraParameters提供了计算3D点投影到像素的函数
      if (z[0]>=0 && z[1]>=0 && z[0]<640 && z[1]<480){ //能不能看到就是说投影完成之后在不在640*480范围内
        ++num_obs;
      }
    }
    if (num_obs>=2){ //被超过1个顶点看到就算是该路标是有效路标，加入路标的顶点中
      optimizer.addVertex(v_p);
      bool inlier = true;
      //生成边
      for (size_t j=0; j<true_poses.size(); ++j){
        //首先，测量值是根据真实的位姿和路标点计算出来的，是真实的测量值。
        Vector2d z
            = cam_params->cam_map(true_poses.at(j).map(true_points.at(i))); //和上面一样，做投影

        if (z[0]>=0 && z[1]>=0 && z[0]<640 && z[1]<480){
          double sam = g2o::Sampler::uniformRand(0., 1.);
          if (sam<OUTLIER_RATIO){ //这里是把当前的观测值做成错误观测.OUTLIER_RATIO默认是0，如果启动程序增加这个值，就会增加outliers的概率。
            z = Vector2d(Sample::uniform(0,640),
                         Sample::uniform(0,480));
            inlier= false;
          }
          z += Vector2d(g2o::Sampler::gaussRand(0., PIXEL_NOISE),
                        g2o::Sampler::gaussRand(0., PIXEL_NOISE)); //这里是正常的添加高斯噪声的观测
          g2o::EdgeProjectXYZ2UV * e
              = new g2o::EdgeProjectXYZ2UV();
          e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(v_p)); //路标点
          e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>
                       (optimizer.vertices().find(j)->second)); //当前观测到第i个路标点的相机顶点
          e->setMeasurement(z); // 设置测量值
          e->information() = Matrix2d::Identity(); //信息矩阵设置为单位阵
          if (ROBUST_KERNEL) {
            g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
            e->setRobustKernel(rk); //设置鲁邦核函数
          }
          e->setParameterId(0, 0);
          optimizer.addEdge(e);//添加边到优化器
        }
      }
  //Step6:保存路标内点的id,初始估计的误差平方和等需要优化之后做对比的参数信息。
      if (inlier){//如果当前的这个路标点的观测是正确的，也就是投影出来的像素坐标没问题
        inliers.insert(point_id); //把路标点的序号保存下来
        Vector3d diff = v_p->estimate() - true_points[i];// 记录路标点初始估计和真实值的差值

        sum_diff2 += diff.dot(diff); //求平方
      }
      pointid_2_trueid.insert(make_pair(point_id,i)); //保存图优化器里面的顶点编号和路标点在路标点vector的配对信息
      ++point_id; //路标点id加1
      ++point_num;//记录添加到优化器里面的路标点数量
    }
  }

  //Step7:开始优化，根据STRUCTURE_ONLY,决定做纯地图点的优化还是同时优化。
  cout << endl;
  optimizer.initializeOptimization();
  optimizer.setVerbose(true);
  //只优化地图点
  if (STRUCTURE_ONLY){
    g2o::StructureOnlySolver<3> structure_only_ba;
    cout << "Performing structure-only BA:"   << endl;
    g2o::OptimizableGraph::VertexContainer points;
    for (g2o::OptimizableGraph::VertexIDMap::const_iterator it = optimizer.vertices().begin(); it != optimizer.vertices().end(); ++it) {
      g2o::OptimizableGraph::Vertex* v = static_cast<g2o::OptimizableGraph::Vertex*>(it->second);
      if (v->dimension() == 3)
        points.push_back(v);
    }
    structure_only_ba.calc(points, 10);
  }
  optimizer.save("ba_test.g2o");
  cout << endl;
  cout << "Performing full BA:" << endl;
  //开始full optimize
  optimizer.optimize(10);
  cout << endl;
  cout << "Point error before optimisation (inliers only): " << sqrt(sum_diff2/inliers.size()) << endl;
  point_num = 0;
  sum_diff2 = 0;
  //Step8:遍历优化后的地图顶点，计算优化后内点的均方误差。
  for (unordered_map<int,int>::iterator it=pointid_2_trueid.begin();
       it!=pointid_2_trueid.end(); ++it){
    g2o::HyperGraph::VertexIDMap::iterator v_it
        = optimizer.vertices().find(it->first);
    if (v_it==optimizer.vertices().end()){
      cerr << "Vertex " << it->first << " not in graph!" << endl;
      exit(-1);
    }
    g2o::VertexSBAPointXYZ * v_p
        = dynamic_cast< g2o::VertexSBAPointXYZ * > (v_it->second);
    if (v_p==0){
      cerr << "Vertex " << it->first << "is not a PointXYZ!" << endl;
      exit(-1);
    }
    Vector3d diff = v_p->estimate()-true_points[it->second];
    if (inliers.find(it->first)==inliers.end())
      continue;
    sum_diff2 += diff.dot(diff);
    ++point_num;
  }

  optimizer.save("ba_test_afterOpt.g2o");
  cout << "Point error after optimisation (inliers only): " << sqrt(sum_diff2/inliers.size()) << endl;
  cout << endl;
}
