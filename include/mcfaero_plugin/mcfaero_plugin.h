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
#ifndef _GAZEBO_MCFAERO_PLUGIN_HH_
#define _GAZEBO_MCFAERO_PLUGIN_HH_

#include <string>
#include <vector>


#include "gazebo/common/Plugin.hh"
#include "gazebo/physics/physics.hh"
#include "gazebo/transport/TransportTypes.hh"
#include <ignition/math.hh>


namespace gazebo
{

  /// \brief A plugin that simulates lift and drag.
  class GAZEBO_VISIBLE McFAeroPlugin : public ModelPlugin
  {
    /// \brief Constructor.
    public: McFAeroPlugin();

    /// \brief Destructor.
    public: ~McFAeroPlugin();

    // Documentation Inherited.
    public: virtual void Load(physics::ModelPtr _model, sdf::ElementPtr _sdf);

    /// \brief Callback for World Update events.
    protected: virtual void OnUpdate();

    /// \brief Connection to World Update events.
    protected: event::ConnectionPtr updateConnection;

    /// \brief Pointer to world.
    protected: physics::WorldPtr world;

    /// \brief Pointer to physics engine.
    protected: physics::PhysicsEnginePtr physics;

    /// \brief Pointer to model containing plugin.
    protected: physics::ModelPtr model;

    /// \brief Pointer to link currently targeted by mud joint.
    protected: physics::LinkPtr link;

    /// \brief Pointer to right aileron joint
    protected: physics::JointPtr raileronJoint;

    /// \brief Pointer to left aileron joint
    protected: physics::JointPtr laileronJoint;

    /// \brief Pointer to elevator joint
    protected: physics::JointPtr elevatorJoint;

    /// \brief Pointer to rudder joint
    protected: physics::JointPtr rudderJoint;

    /// \brief Pointer to motor joint
    protected: physics::JointPtr motorJoint;

    /// \brief SDF for this plugin;
    protected: sdf::ElementPtr sdf;
    ///




  };
}
#endif
