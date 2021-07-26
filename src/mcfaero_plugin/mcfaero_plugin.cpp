/*
 * Copyright (C) 2014-2016 Open Source Robotics Foundation
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
 *
*/

#include <algorithm>
#include <string>

#include "common.h"
#include "gazebo/common/Assert.hh"
#include "gazebo/physics/physics.hh"
#include "gazebo/sensors/SensorManager.hh"
#include "gazebo/transport/transport.hh"
#include "gazebo/msgs/msgs.hh"
#include "mcfaero_plugin/mcfaero_plugin.h"
#include "mcfoamy_v2/McFoamy_FM_v2.h"

using namespace gazebo;
using namespace ignition::math;


GZ_REGISTER_MODEL_PLUGIN(McFAeroPlugin)

/////////////////////////////////////////////////
McFAeroPlugin::McFAeroPlugin()
{
}

/////////////////////////////////////////////////
McFAeroPlugin::~McFAeroPlugin()
{
}

/////////////////////////////////////////////////
void McFAeroPlugin::Load(physics::ModelPtr _model,
                     sdf::ElementPtr _sdf)
{
  GZ_ASSERT(_model, "McFAeroPlugin _model pointer is NULL");
  GZ_ASSERT(_sdf, "McFAeroPlugin _sdf pointer is NULL");
  this->model = _model;
  this->sdf = _sdf;

  this->world = this->model->GetWorld();
  GZ_ASSERT(this->world, "McFAeroPlugin world pointer is NULL");

  this->physics = this->world->Physics();

  GZ_ASSERT(this->physics, "McFAeroPlugin physics pointer is NULL");

  GZ_ASSERT(_sdf, "McFAeroPlugin _sdf pointer is NULL");


  if (_sdf->HasElement("link_name"))
  {
    sdf::ElementPtr elem = _sdf->GetElement("link_name");
    // GZ_ASSERT(elem, "Element link_name doesn't exist!");
    std::string linkName = elem->Get<std::string>();
    this->link = this->model->GetLink(linkName);
    // GZ_ASSERT(this->link, "Link was NULL");

    if (!this->link)
    {
      gzerr << "Link with name[" << linkName << "] not found. "
        << "The McFAeroPlugin will not generate forces\n";
    }
    else
    {
      this->updateConnection = event::Events::ConnectWorldUpdateBegin(
          boost::bind(&McFAeroPlugin::OnUpdate, this));
    }
  }

  if (_sdf->HasElement("raileron_joint_name"))
  {
    std::string raileronJointName = _sdf->Get<std::string>("raileron_joint_name");
    this->raileronJoint = this->model->GetJoint(raileronJointName);
    if (!this->raileronJoint)
    {
      gzerr << "Joint with name[" << raileronJointName << "] does not exist.\n";
    }
  }

  if (_sdf->HasElement("laileron_joint_name"))
  {
    std::string laileronJointName = _sdf->Get<std::string>("laileron_joint_name");
    this->laileronJoint = this->model->GetJoint(laileronJointName);
    if (!this->laileronJoint)
    {
      gzerr << "Joint with name[" << laileronJointName << "] does not exist.\n";
    }
  }

  if (_sdf->HasElement("elevator_joint_name"))
  {
    std::string elevatorJointName = _sdf->Get<std::string>("elevator_joint_name");
    this->elevatorJoint = this->model->GetJoint(elevatorJointName);
    if (!this->elevatorJoint)
    {
      gzerr << "Joint with name[" << elevatorJointName << "] does not exist.\n";
    }
  }

  if (_sdf->HasElement("rudder_joint_name"))
  {
    std::string rudderJointName = _sdf->Get<std::string>("rudder_joint_name");
    this->rudderJoint = this->model->GetJoint(rudderJointName);
    if (!this->rudderJoint)
    {
      gzerr << "Joint with name[" << rudderJointName << "] does not exist.\n";
    }
  }

  if (_sdf->HasElement("motor_joint_name"))
  {
    std::string motorJointName = _sdf->Get<std::string>("motor_joint_name");
    this->motorJoint = this->model->GetJoint(motorJointName);

    if (!this->motorJoint)
    {
      gzerr << "Joint with name[" << motorJointName << "] does not exist.\n";
    }
  }






}

/////////////////////////////////////////////////
void McFAeroPlugin::OnUpdate()
{
  GZ_ASSERT(this->link, "Link was NULL");

  // get linear velocity at cp in inertial frame
  //Vector3d vel = this->link->WorldLinearVel();
//  Vector3d velI = vel;
//  velI.Normalize();
  //Vector3d V_g = vel;



  // if (vel.Length() <= 0.01)
  //   return;

  // pose of body
  Pose3d pose = this->link->WorldPose();
  Quaterniond q_gr = pose.Rot(); // body att. in Gazebo frames (from FLU to ENU)
  Quaterniond q_gb = q_gr*q_FLU_to_FRD.Inverse(); // from FRD to ENU
  Quaterniond q_nb = q_ENU_to_NED*q_gb; // form FRD to NED
  Vector3d Eul = q_nb.Euler();
  Vector3d omega_nb_b = q_FLU_to_FRD.RotateVector(this->link->RelativeAngularVel());

  // Body velocity in PX4 frames (NED, FRD)
  //Vector3d V_b = q_gb.Inverse().RotateVector(V_g);
  //Vector3d V_w = Vector3d(0,0,0); //placeholder wind velocity for later

  //Vector3d Vel_b = V_b-V_w;

  Vector3d vel_cg0_b = q_FLU_to_FRD.RotateVector(this->link->RelativeLinearVel());
  Vector3d vel_wind = Vector3d(0,0,0); //placeholder wind velocity for later
  Vector3d vel_aspd = vel_cg0_b - vel_wind;

  double vel_M = vel_aspd.Length();
  // Testing function
  double LAilDef = -1 * this->laileronJoint->Position(0);
  double ElevDef = -1 * this->elevatorJoint->Position(0);
  double RudDef = -1 * this->rudderJoint->Position(0);
  double mtr_spd = this->motorJoint->GetVelocity(0);

  mtr_spd = floorf(10*fabs(mtr_spd));


  double wIn = mtr_spd;
  double u_mcf = vel_aspd.X();
  double v_mcf = vel_aspd.Y();
  double w_mcf = vel_aspd.Z();
  double p_mcf = omega_nb_b.X();
  double q_mcf = omega_nb_b.Y();
  double r_mcf = omega_nb_b.Z();

  //Output variables
  double Fx;
  double Fy;
  double Fz;
  double Mx;
  double My;
  double Mz;
  McFoamy_FM_v2(LAilDef, ElevDef, RudDef, wIn, u_mcf, v_mcf, w_mcf, p_mcf,
                q_mcf, r_mcf, &Fx, &Fy, &Fz, &Mx, &My, &Mz);


  Vector3d Forces_mcf = Vector3d(Fx,Fy,Fz);
  Vector3d Moments_mcf = Vector3d(Mx,My,Mz);

//  Vector3d Forces_mcf = Vector3d(0,0,0);
//  Vector3d Moments_mcf = Vector3d(0,0,0);

  Vector3d force_gzb = q_gb*Forces_mcf;
  Vector3d moment_gzb = q_gb*Moments_mcf;

  // Correct for nan or inf
  force_gzb.Correct();
  moment_gzb.Correct();




  if (0) //mtr_spd > 0 vel_M > 5 vel_M > 205
  {
    gzdbg << "=============================\n";
//     gzdbg << "sensor: [" << this->GetHandle() << "]\n";
//     gzdbg << "Link: [" << this->link->GetName() << "]\n" << "Pose: [" << pose << "]\n";
//
    // gzdbg << "R_aileron control: [" << this->raileronJoint->Position(0) << "]\n";
    // gzdbg << "L_aileron control: [" << this->laileronJoint->Position(0) << "]\n";
    // gzdbg << "Elevator control: [" << this->elevatorJoint->Position(0) << "]\n";
    gzdbg << "Rudder control: [" << this->rudderJoint->Position(0) << "]\n";
    // gzdbg << "Rotor control: [" << this->motorJoint->GetVelocity(0) << "]\n";
// //    gzdbg << "Euler angles: [" << Eul << "]\n";
//     gzdbg << "Vel_b: [" << vel_cg0_b << "]\n";
//     gzdbg << "Angular velocity: [" << omega_nb_b << "]\n";
    // gzdbg << "Forces: [" << Forces_mcf << "]\n";
    // gzdbg << "Moments: [" << Moments_mcf << "]\n";
    // gzdbg << "Test: [" << mtr_spd << "]\n";
  }
  // apply forces at cg (with torques for position shift)
  this->link->AddForce(force_gzb);
  this->link->AddTorque(moment_gzb);

}
