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

#include "simulator.h"

#include "g2o/stuff/sampler.h"

#include <map>
#include <iostream>
#include <cmath>
using namespace std;

namespace g2o {
  namespace tutorial {

    using namespace Eigen;

#  ifdef _MSC_VER
    inline double round(double number)
    {
      return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
    }
#  endif

    typedef std::map<int, std::map<int, Simulator::LandmarkPtrVector> > LandmarkGrid;

    Simulator::Simulator()
    {
      time_t seed = time(0);
      Sampler::seedRand(static_cast<unsigned int>(seed));
    }

    Simulator::~Simulator()
    {
    }

    void Simulator::simulate(int numNodes, const SE2& sensorOffset)
    {
      // simulate a robot observing landmarks while travelling on a grid
      int steps = 5; //走5步直线，转弯一次
      double stepLen = 1.0; //机器人相邻时间的运动步长
      int boundArea = 50; //机器人运动的范围

      double maxSensorRangeLandmarks = 2.5 * stepLen; // 机器人步长的2.5倍。

      int landMarksPerSquareMeter = 1; //路标点在每一平方米内包含的数量
      double observationProb = 0.8; //机器人能够看到某个路标点的概率

      int landmarksRange=2;

      Vector2d transNoise(0.05, 0.01);//机器人运动时候的平移噪声
      double rotNoise = DEG2RAD(2.);//机器人旋转的噪声
      Vector2d landmarkNoise(0.05, 0.05);

      Vector2d bound(boundArea, boundArea);

      VectorXd probLimits;
      probLimits.resize(MO_NUM_ELEMS);//转弯的方向，左转或者右转
      cout << "MO_NUM_ELEMS: " << MO_NUM_ELEMS <<endl;
      //采样左转或者右转,如果左转，那么limit范围是0~0.5，否则右转就是0.5~1.0
      for (int i = 0; i < probLimits.size(); ++i) {
        probLimits[i] = (i + 1) / (double) MO_NUM_ELEMS;
        cout << "probLimits: " << probLimits[i] << endl;
      }

      Matrix3d covariance;
      covariance.fill(0.);
      covariance(0, 0) = transNoise[0]*transNoise[0]; //传感器的测量误差x方向的方差
      covariance(1, 1) = transNoise[1]*transNoise[1];//传感器的测量误差y方向的方差
      covariance(2, 2) = rotNoise*rotNoise;//传感器的测量误差角度的方差
      Matrix3d information = covariance.inverse(); //协方差矩阵的逆就是信息矩阵，越大，说明越可信

      SE2 maxStepTransf(stepLen * steps, 0, 0); //x方向单次直线的最长的距离
      Simulator::PosesVector& poses = _poses;
      poses.clear();
      LandmarkVector& landmarks = _landmarks;
      landmarks.clear();
      Simulator::GridPose firstPose;
      firstPose.id = 0;
      firstPose.truePose = SE2(0,0,0);
      firstPose.simulatorPose = SE2(0,0,0);
      poses.push_back(firstPose);
      cerr << "Simulator: sampling nodes ...";
        /*********************
        * Step1：创建机器人真实运动的轨迹（离散的2D坐标和朝向角度）
        */
      while ((int)poses.size() < numNodes) {
        // add straight motions
        //在五步之内，添加当前方向的直行位姿
        for (int i = 1; i < steps && (int)poses.size() < numNodes; ++i) {
          Simulator::GridPose nextGridPose = generateNewPose(poses.back(), SE2(stepLen,0,0), transNoise, rotNoise);
          poses.push_back(nextGridPose);
        }
        if ((int)poses.size() == numNodes)
          break;

        // sample a new motion direction
        //通过随机数来决定是左转还是右转
        double sampleMove = Sampler::uniformRand(0., 1.);
        int motionDirection = 0;
        while (probLimits[motionDirection] < sampleMove && motionDirection+1 < MO_NUM_ELEMS) {
          motionDirection++;
        }

        SE2 nextMotionStep = getMotion(motionDirection, stepLen);
        Simulator::GridPose nextGridPose = generateNewPose(poses.back(), nextMotionStep, transNoise, rotNoise);

        // check whether we will walk outside the boundaries in the next iteration
        //根据转弯之后的角度，计算下一个五步，会不会超越规定的机器人运动范围，如果两个方向都没办法保证在范围内，就停止机器人的运动。
        SE2 nextStepFinalPose = nextGridPose.truePose * maxStepTransf;
        if (fabs(nextStepFinalPose.translation().x()) >= bound[0] || fabs(nextStepFinalPose.translation().y()) >= bound[1]) {
          //cerr << "b";
          // will be outside boundaries using this
          for (int i = 0; i < MO_NUM_ELEMS; ++i) {
            nextMotionStep = getMotion(i, stepLen);
            nextGridPose = generateNewPose(poses.back(), nextMotionStep, transNoise, rotNoise);
            nextStepFinalPose = nextGridPose.truePose * maxStepTransf;
            if (fabs(nextStepFinalPose.translation().x()) < bound[0] && fabs(nextStepFinalPose.translation().y()) < bound[1])
              break;
          }
        }

        poses.push_back(nextGridPose);
      }
      cerr << "done." << endl;
      //print poses for understanding the model
      for(int k = 0; k < poses.size(); k++) {
        cout << "The " << k << "th simulated SE(2) pose is: " << 
              poses[k].simulatorPose[0] << ", " << 
              poses[k].simulatorPose[1] << ", " << 
              poses[k].simulatorPose[2]<< endl;

        cout << "The " << k << "th true SE(2) pose is: " << 
              poses[k].truePose[0] << ", " << 
              poses[k].truePose[1] << ", " << 
              poses[k].truePose[2]<< endl;
      }

      // creating landmarks along the trajectory
      cerr << "Simulator: Creating landmarks ... ";
        /*********************
        * Step2:依据poses创建LandmarkGrid
        */
      LandmarkGrid grid; //路标的二维数组，二维格子里面包含的路标信息
      for (PosesVector::const_iterator it = poses.begin(); it != poses.end(); ++it) {
        int ccx = (int)round(it->truePose.translation().x());//机器人x坐标取整数
        int ccy = (int)round(it->truePose.translation().y());//机器人y坐标取整数
        for (int a=-landmarksRange; a<=landmarksRange; a++) //landmarksRange是2，代表了机器人周围边长是2范围内的大正方形
          for (int b=-landmarksRange; b<=landmarksRange; b++){
            int cx=ccx+a;//每一个格子的x方向索引号,也是格子最靠近坐标原点的顶点坐标->第一象限就是左下角，第二象限就是右下角，第三象限就是右上角，第四象限就是左上角。
            int cy=ccy+b;//每一个格子的y方向索引号
            LandmarkPtrVector& landmarksForCell = grid[cx][cy]; //按照行、列从grid里面取出。下面是如何更新这一个格子里面的路标点
            if (landmarksForCell.size() == 0) {
              for (int i = 0; i < landMarksPerSquareMeter; ++i) {//landMarksPerSquareMeter==1
                Landmark* l = new Landmark(); //在堆中创建出来的路标点 后面会整体从内存中删除
                double offx, offy;
                do {
                  offx = Sampler::uniformRand(-0.5*stepLen, 0.5*stepLen); //-0.5 ~ +0.5
                  offy = Sampler::uniformRand(-0.5*stepLen, 0.5*stepLen); //-0.5 ~ +0.5
                } while (hypot_sqr(offx, offy) < 0.25*0.25); //随机数要在0.25的范围内。是一个圆形范围
                l->truePose[0] = cx + offx;
                l->truePose[1] = cy + offy;
                landmarksForCell.push_back(l);
              }
            }
          }
      }
      cerr << "done." << endl;
      


      cerr << "Simulator: Simulating landmark observations for the poses ... ";
      /*********************
       * Step3:
      * 根据创建出来的机器人顶点和路标点，用传感器的感知范围过滤出来每一个机器人顶点观测到的路标点
      *
      */
      double maxSensorSqr = maxSensorRangeLandmarks * maxSensorRangeLandmarks; //机器人传感器最大扫描范围
      int globalId = 0;
      for (PosesVector::iterator it = poses.begin(); it != poses.end(); ++it) { //遍历每一个机器人点
        Simulator::GridPose& pv = *it;
        int cx = (int)round(it->truePose.translation().x());
        int cy = (int)round(it->truePose.translation().y());
        int numGridCells = (int)(maxSensorRangeLandmarks) + 1;

        pv.id = globalId++;
        SE2 trueInv = pv.truePose.inverse();

        for (int xx = cx - numGridCells; xx <= cx + numGridCells; ++xx)
          for (int yy = cy - numGridCells; yy <= cy + numGridCells; ++yy) { //站在机器人当前的位置，通过grid遍历每一个格子上面的路标向量
            LandmarkPtrVector& landmarksForCell = grid[xx][yy];
            if (landmarksForCell.size() == 0)
              continue;
            for (size_t i = 0; i < landmarksForCell.size(); ++i) { //判断机器人到路标的距离是不是在传感器的范围maxSensorSqr
              Landmark* l = landmarksForCell[i];
              double dSqr = hypot_sqr(pv.truePose.translation().x() - l->truePose.x(), pv.truePose.translation().y() - l->truePose.y());
              if (dSqr > maxSensorSqr)
                continue;
              double obs = Sampler::uniformRand(0.0, 1.0);
              if (obs > observationProb) // we do not see this one...
                continue; //在这里可以根据随机生成的obs值，来确定要不要加入这个路标点
              if (l->id < 0)
                l->id = globalId++;
              if (l->seenBy.size() == 0) {//之前没有被看到过
                Vector2d trueObservation = trueInv * l->truePose;//机器人对该路标的观测值，机器人坐标系下
                Vector2d observation = trueObservation;//机器人的真实观测需要添加一定的高斯噪声
                observation[0] += Sampler::gaussRand(0., landmarkNoise[0]);
                observation[1] += Sampler::gaussRand(0., landmarkNoise[1]);
                l->simulatedPose = pv.simulatorPose * observation; //更新路标在世界坐标系下的坐标（噪声添加好了）
              }//如果之前被看到了，那么，就不需要再更新
              l->seenBy.push_back(pv.id); //更新被当前相机看到了
              pv.landmarks.push_back(l); //该相机看到了这个路标，
            }
          }

      }
      cerr << "done." << endl;
      // add the odometry measurements
      /*********************
       * Step4:
      * 这里是添加前后两帧机器人顶点的转换
      */
      _odometry.clear();
      cerr << "Simulator: Adding odometry measurements ... ";
      for (size_t i = 1; i < poses.size(); ++i) {
        const GridPose& prev = poses[i-1]; //当前相机顶点的前一个顶点
        const GridPose& p = poses[i]; //当前顶点

        _odometry.push_back(GridEdge());
        GridEdge& edge = _odometry.back();

        edge.from = prev.id; //赋值前一个顶点id
        edge.to = p.id; //赋值当前顶点id
        edge.trueTransf = prev.truePose.inverse() * p.truePose; //旧姿态转换到新姿态的变换。后面根据g2o怎么用的确定一下
        edge.simulatorTransf = prev.simulatorPose.inverse() * p.simulatorPose;
        edge.information = information;
      }
      cerr << "done." << endl;

      _landmarks.clear();
      _landmarkObservations.clear();
      // add the landmark observations
      //Step5:
      //这里是添加相机到路标之前的观测值
      {
        cerr << "Simulator: add landmark observations ... ";
          /*********************
          * 从机器人顶点里面抽取不重复的路标点
          */
        Matrix2d covariance; covariance.fill(0.);
        covariance(0, 0) = landmarkNoise[0]*landmarkNoise[0];
        covariance(1, 1) = landmarkNoise[1]*landmarkNoise[1];
        Matrix2d information = covariance.inverse();

        for (size_t i = 0; i < poses.size(); ++i) {
          const GridPose& p = poses[i];
          for (size_t j = 0; j < p.landmarks.size(); ++j) {
            Landmark* l = p.landmarks[j];
            if (l->seenBy.size() > 0 && l->seenBy[0] == p.id) { //如果该相机顶点下的第一个路标是本顶点，那么，保存路标到landmarks
              landmarks.push_back(*l); //这里的landmarks不会重复
            }
          }
        }

         /*********************
          * 遍历每一个机器人节点，节点中保存了该节点看到的路标
          * 如果这个路标是第一次被这个节点看到，观测值的计算是
          * 通过带噪声机器人位姿和带噪声的路标点得到。再往后的顶点看到
          * 的坐标就是真实值加高斯噪声
         */
        for (size_t i = 0; i < poses.size(); ++i) {
          const GridPose& p = poses[i];
          SE2 trueInv = (p.truePose * sensorOffset).inverse(); //Tw->r(世界坐标系到机器人坐标系转换矩阵) * Tr->s(机器人坐标系到传感器坐标系转换)
          for (size_t j = 0; j < p.landmarks.size(); ++j) {
            Landmark* l = p.landmarks[j];
            Vector2d observation;
            Vector2d trueObservation = trueInv * l->truePose; // Tw->s * p_landmark = p_land_in_sensor 这里计算的是路标点在传感器中的坐标
            observation = trueObservation; //传感器测量到的路标点的真值
            if (l->seenBy.size() > 0 && l->seenBy[0] == p.id) { // write the initial position of the landmark
              observation = (p.simulatorPose * sensorOffset).inverse() * l->simulatedPose; // 路标点第一次被看到时候的顶点下的坐标当做初始值
            } else {//如果不是第一次被看到，那么，真值加高斯噪声
              // create observation for the LANDMARK using the true positions
              observation[0] += Sampler::gaussRand(0., landmarkNoise[0]);
              observation[1] += Sampler::gaussRand(0., landmarkNoise[1]);
            }

            _landmarkObservations.push_back(LandmarkEdge()); //路标观测边
            LandmarkEdge& le = _landmarkObservations.back();

            le.from = p.id;//顶点id
            le.to = l->id; //路标点id
            le.trueMeas = trueObservation; //机器人上传感器测量路标点的真实值
            le.simulatorMeas = observation; //添加了高斯白噪声的观测值
            le.information = information; //信息矩阵
          }
        }
        cerr << "done." << endl;
      }


      // cleaning up
        /*********************
         * Step6:
        * 从堆内存中删除掉LandmarkGrid grid
           */
      for (LandmarkGrid::iterator it = grid.begin(); it != grid.end(); ++it) {
        for (std::map<int, Simulator::LandmarkPtrVector>::iterator itt = it->second.begin(); itt != it->second.end(); ++itt) {
          Simulator::LandmarkPtrVector& landmarks = itt->second;
          for (size_t i = 0; i < landmarks.size(); ++i)
            delete landmarks[i];
        }
      }

    }

    Simulator::GridPose Simulator::generateNewPose(const Simulator::GridPose& prev, const SE2& trueMotion, const Eigen::Vector2d& transNoise, double rotNoise)
    {
      Simulator::GridPose nextPose;
      nextPose.id = prev.id + 1;
      nextPose.truePose = prev.truePose * trueMotion;
      SE2 noiseMotion = sampleTransformation(trueMotion, transNoise, rotNoise);
      nextPose.simulatorPose = prev.simulatorPose * noiseMotion;
      return nextPose;
    }

    SE2 Simulator::getMotion(int motionDirection, double stepLen)
    {
      switch (motionDirection) {
        case MO_LEFT:
          return SE2(stepLen, 0, 0.5*M_PI);
        case MO_RIGHT:
          return SE2(stepLen, 0, -0.5*M_PI);
        default:
          cerr << "Unknown motion direction" << endl;
          return SE2(stepLen, 0, -0.5*M_PI);
      }
    }

    SE2 Simulator::sampleTransformation(const SE2& trueMotion_, const Eigen::Vector2d& transNoise, double rotNoise)
    {
      Vector3d trueMotion = trueMotion_.toVector();
      SE2 noiseMotion(
          trueMotion[0] + Sampler::gaussRand(0.0, transNoise[0]),
          trueMotion[1] + Sampler::gaussRand(0.0, transNoise[1]),
          trueMotion[2] + Sampler::gaussRand(0.0, rotNoise));
      return noiseMotion;
    }

  } // end namespace
} // end namespace
