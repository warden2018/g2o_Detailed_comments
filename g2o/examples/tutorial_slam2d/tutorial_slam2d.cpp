// g2o - General Graph Optimization
// Copyright (C) 2011 R. Kuemmerle, G. Grisetti, W. Burgard
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
#include <cmath>

#include "simulator.h"

#include "vertex_se2.h"
#include "vertex_point_xy.h"
#include "edge_se2.h"
#include "edge_se2_pointxy.h"
#include "types_tutorial_slam2d.h"

#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/factory.h"
#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/solvers/eigen/linear_solver_eigen.h"

using namespace std;
using namespace g2o;
using namespace g2o::tutorial;

int main()
{
  // TODO simulate different sensor offset
  // simulate a robot observing landmarks while travelling on a grid
  //Step1: 通过仿真器生成优化所需要的机器人位姿，地图点坐标和机器人的观测值
  SE2 sensorOffsetTransf(0.2, 0.1, -0.1); //是传感器到机器人的外参数
  int numNodes = 300; //机器人顶点数
  Simulator simulator;
  simulator.simulate(numNodes, sensorOffsetTransf);

  /*********************************************************************************
   * creating the optimization problem
   ********************************************************************************/
  //Step2:创建优化器，设置优化算法为GN,设置线性求解器为动态的
  typedef BlockSolver< BlockSolverTraits<-1,-1> >  SlamBlockSolver; //这里不明白-1 -1 代表什么？是动态的定义顶点的维度吗？ -- 这里是的-1，-1是Eigen的Dynamic定义出来的值
  /*
  template <>
  struct BlockSolverTraits<Eigen::Dynamic, Eigen::Dynamic>
  {
    static const int PoseDim = Eigen::Dynamic;
    static const int LandmarkDim = Eigen::Dynamic;
    typedef MatrixX PoseMatrixType;
    typedef MatrixX LandmarkMatrixType;
    typedef MatrixX PoseLandmarkMatrixType;
    typedef VectorX PoseVectorType;
    typedef VectorX LandmarkVectorType;

    typedef SparseBlockMatrix<PoseMatrixType> PoseHessianType;
    typedef SparseBlockMatrix<LandmarkMatrixType> LandmarkHessianType;
    typedef SparseBlockMatrix<PoseLandmarkMatrixType> PoseLandmarkHessianType;
    typedef LinearSolver<PoseMatrixType> LinearSolverType;
  };
  */
  typedef LinearSolverEigen<SlamBlockSolver::PoseMatrixType> SlamLinearSolver; //这里通过BlockSolver的PoseMatrixType类型定义了SlamLinearSolver的类型

  // allocating the optimizer
  SparseOptimizer optimizer; // 实例化一个SparseOptimizer
  auto linearSolver = g2o::make_unique<SlamLinearSolver>(); // 实例化线性求解器SlamLinearSolver
  linearSolver->setBlockOrdering(false);
  OptimizationAlgorithmGaussNewton* solver = new OptimizationAlgorithmGaussNewton(
    g2o::make_unique<SlamBlockSolver>(std::move(linearSolver))); // SlamBlockSolver类型的solver作为参数初始化OptimizationAlgorithmGaussNewton

  optimizer.setAlgorithm(solver); // SparseOptimizer 设置优化器

  // add the parameter representing the sensor offset
  //Step3:创建传感器的外参数，代表了相机在机器人上的安装位置
  ParameterSE2Offset* sensorOffset = new ParameterSE2Offset;
  sensorOffset->setOffset(sensorOffsetTransf);
  sensorOffset->setId(0);
  optimizer.addParameter(sensorOffset);

  // adding the odometry to the optimizer
  // first adding all the vertices
  cerr << "Optimization: Adding robot poses ... ";

  //Step4:添加机器人运动顶点，机器人里程计的边，地图点顶点，机器人对地图点的观测边到图中

  /*
   * 添加机器人运动顶点到优化器
   */
  for (size_t i = 0; i < simulator.poses().size(); ++i) {
    const Simulator::GridPose& p = simulator.poses()[i];
    const SE2& t = p.simulatorPose;
    VertexSE2* robot =  new VertexSE2;
    robot->setId(p.id); // 设置机器人顶点的id,和仿真里面生成的机器人pose id相同
    robot->setEstimate(t); // 设置估计值
    optimizer.addVertex(robot); //添加顶点到优化器当中
  }
  cerr << "done." << endl;

  // second add the odometry constraints
  cerr << "Optimization: Adding odometry measurements ... ";
  /*
   * 添加里程计边到图优化器当中
   */
  for (size_t i = 0; i < simulator.odometry().size(); ++i) {
    const Simulator::GridEdge& simEdge = simulator.odometry()[i];

    EdgeSE2* odometry = new EdgeSE2;
    odometry->vertices()[0] = optimizer.vertex(simEdge.from); //设置边的第一个顶点id
    odometry->vertices()[1] = optimizer.vertex(simEdge.to); //设置边的第二个顶点id
    odometry->setMeasurement(simEdge.simulatorTransf); //设置前后顶点的转换
    odometry->setInformation(simEdge.information); // 设置信息矩阵
    optimizer.addEdge(odometry); // 添加边
  }
  cerr << "done." << endl;

  // add the landmark observations
  cerr << "Optimization: add landmark vertices ... ";
  /*
   * 添加路标到图优化器
   */
  for (size_t i = 0; i < simulator.landmarks().size(); ++i) {
    const Simulator::Landmark& l = simulator.landmarks()[i];
    VertexPointXY* landmark = new VertexPointXY;
    landmark->setId(l.id);
    landmark->setEstimate(l.simulatedPose);
    optimizer.addVertex(landmark);
  }
  cerr << "done." << endl;

  cerr << "Optimization: add landmark observations ... ";
  /*
   * 添加机器人对路标的观测边到优化器
   */
  for (size_t i = 0; i < simulator.landmarkObservations().size(); ++i) {
    const Simulator::LandmarkEdge& simEdge = simulator.landmarkObservations()[i];
    EdgeSE2PointXY* landmarkObservation =  new EdgeSE2PointXY;
    landmarkObservation->vertices()[0] = optimizer.vertex(simEdge.from);
    landmarkObservation->vertices()[1] = optimizer.vertex(simEdge.to);
    landmarkObservation->setMeasurement(simEdge.simulatorMeas);
    landmarkObservation->setInformation(simEdge.information);
    landmarkObservation->setParameterId(0, sensorOffset->id());
    optimizer.addEdge(landmarkObservation);
  }
  cerr << "done." << endl;


  /*********************************************************************************
   * optimization
   ********************************************************************************/
  //Step5:开始图优化
  // dump initial state to the disk
  optimizer.save("tutorial_before.g2o");

  // prepare and run the optimization
  // fix the first robot pose to account for gauge freedom
  VertexSE2* firstRobotPose = dynamic_cast<VertexSE2*>(optimizer.vertex(0));
  firstRobotPose->setFixed(true);
  optimizer.setVerbose(true);

  cerr << "Optimizing" << endl;
  optimizer.initializeOptimization();
  optimizer.optimize(10);
  cerr << "done." << endl;
  //Step6:将优化结果写到文件中，并且删除图，清理内存。结束程序。
  optimizer.save("tutorial_after.g2o");
  
  // freeing the graph memory
  optimizer.clear();

  return 0;
}
