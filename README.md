# Path Planning

[![Udacity - Self-Driving Car NanoDegree](https://s3.amazonaws.com/udacity-sdc/github/shield-carnd.svg)](http://www.udacity.com/drive)

---

The goals of this project are the following:

* Implement Prediction, Behavior Planner and Trajectory modules for highway driving.
* Use a Finite State Machine for the Behavior Planning.
* Test implementation by driving a car on a highway (in a simulator), maintaining a speed close to 50 mph and safely changing lanes to go around slower traffic.

---

## 1. Files

My project includes the following files:

- [<b>C++</b> - The source code](https://github.com/ArjaanBuijk/CarND-Path-Planning-Project/tree/master/src)
- [<b>README.md</b> - A summary of the project](https://github.com/ArjaanBuijk/CarND-Path-Planning-Project/blob/master/README.md)
- [<b>video_MPC.mp4</b> - Video showing the car driving on the highway for 10 miles without incidents](https://github.com/ArjaanBuijk/CarND-Path-Planning-Project/blob/master/videoProject.mp4)<br>(short version showing first 1.75 miles is embedded here as an animated gif)
    ![](https://github.com/ArjaanBuijk/CarND-Path-Planning-Project/blob/master/images/videoProject.gif?raw=true)



## 2. Note on naming convention
In this document and in the code, the following naming convention is used to distinguis between the car are we are driving and the surrounding traffic:

|Name|Refers to|
|-|-|
|car|The car we are driving|
|object|The other vehicles on the road|

The software is written using a functional design, not with classes. To group the functions that belong to each module, they are named with an common prefix, as follows:

|Module| Function name prefix|
|-|-|
|Prediction| predict\___ |
|Behavior Planner| behavior\___ |
|Trajectory| trajectory\___ |


## 3. Project Description
The software modules in a self driving car and the data flow between them are:
<div style="text-align:center"><img src="https://github.com/ArjaanBuijk/CarND-Path-Planning-Project/blob/master/images/modules.gif?raw=true" style="width:500px;"/></div>

The responsibility of each module, and their inputs and outputs are summarized in this table.

|Module|Responsibility|Input|Output|
|-|-|-|-|
|Motion Control|Drive car along provided trajectory|Trajectory that car should follow|Actuator controls (Steering, Throttle, Brake,...)|
|Sensor Fusion|Detect stationary and moving objects around car|Sensor data|Objects location & movement|
|Localization|Determine location, orientation and speed of car|Maps, GPS & Sensor Fused data|Location, orientation and speed of car|
|Prediction|Predict future location and motion of objects|Maps, Sensor Fusion & Localization|Trajectories of objects with associated probability|
|Behavior Planner|Suggest maneuvers which are Feasible, Safe, Legal & Efficient|Maps, Route, Localization & Prediction|Suggested maneuvers|
|Trajectory|For each suggested maneuver, define a trajectory and select the best trajectory that the car should follow|Maps, Localization, Prediction & Behavior |Trajectory that car must follow|


In this project, a highway map is provided as waypoints. The highway is 6945.554m long and looks like this in global xy space:
<div style="text-align:center"><img src="https://github.com/ArjaanBuijk/CarND-Path-Planning-Project/blob/master/images/highway-map-xy.gif?raw=true" style="width:500px;"/></div>

The objective of our project is to write the <b>Behavior</b>, <b>Prediction</b> and <b>Trajectory</b> modules, and provide the simulator with a trajectory for the car to follow. The trajectory is defined by points along a path with a certain spacing. The spacing indicates the distance the car must drive in 0.02 seconds and thus defines the speed along the trajectory.

The simulator uses this trajectory and takes care of the <b>Motion Control</b>, so that the car indeed follows the trajectory at the defined speed.

The simulator also takes care of <b>Sensor Fusion</b> and <b>Localization</b>, and each cycle it provides the results back, so we can update the trajectory.

<b>Localization data provided by simulator</b>

|name|Description|
|-|-|
|car_x| The car's x position in map coordinates|
|car_y| The car's y position in map coordinates|
|car_s| The car's s position in frenet coordinates|
|car_d| The car's d position in frenet coordinates|
|car_yaw| The car's yaw angle in the map. Provided in degrees, but we store it in radians.|
|car_speed| The car's speed. Provided in MPH, but we store it in m/s.|

<b>Sensor Fusion data provided by simulator</b>

A list of objects on the same side of the highway. For each object, the following data is provided:

|name|Description|
|-|-|
|object_id|The object's unique ID|
|object_x|The object's x position in map coordinates|
|object_y|The object's y position in map coordinates|
|object_vx|The object's x velocity in m/s|
|object_vy|The object's y velocity in m/s|
|object_s|The object's s position in frenet coordinates|
|object_d|The object's d position in frenet coordinates|

<b>Previous path data provided by simulator</b>

The simulator returns the remainder of the trajectory that the car has not yet traversed. It provides the following data:

|name|Description|
|-|-|
|previous_path_x|x of remainder of previous trajectory that car has not yet traversed|
|previous_path_y|y of remainder of previous trajectory that car has not yet traversed|
|end_path_s|Frenet s value at end of previous trajectory|
|end_path_d|Frenet d value at end of previous trajectory|


## 4. Description of Prediction Module

The Prediction Module is implemented in a single function:  <b>predict___objects</b>

Since all objects behave very consistent, our prediction logic could be kept extremely simple. The objects are assumed to continue to drive with constant velocity, staying in their lane.

This is a simplification, but it works well for this project.

If objects are introduced that exhibit more erratic behavior, then our Prediction module must be updated, for example with a Gaussian Naive Bayes predictor.


## 5. Description of Behavior Module


We implemented a conservative driving behavior, making sure to keep our distance to objects in front and that there is sufficient space before making a lane change. This allowed us to use the simple prediction approach described above.

The Behavior Module is implemented as a Finite State Machine, with 3 states:

|State|Description|
|-|-|
|KL| Keep Lane|
|LCL| Lane Change Left|
|LCR| Lane Change Right|


The states and the available transitions are summarized in this diagram:
<div style="text-align:center"><img src="https://github.com/ArjaanBuijk/CarND-Path-Planning-Project/blob/master/images/Finite-State-Machine-diagram.gif?raw=true" style="width:500px;"/></div>

Each state supports self transition, and there are no accepting states. A transition directly from LCL to LCR is not allowed. It always has to go over the KL state.

The following functions take care of the behavior planning.

|Function| Description|
|-|-|
|behavior___possible_successor_states| If told to <b>stay_in_lane (*)</b>, the only successor state is KL.<br> Else the possible successor states are as defined in the Finite State Machine diagram above.<br> The function returns a vector of possible next states|
|behavior___cost| Calculate total cost to transition to next_state, as the sum of several cost functions|
|behavior___cost_collission_danger|We strongly penalize lane changes that have a collission danger.<br> <b>WEIGHT_SAFETY_COLLISSION_COST</b> is set so high, at 1.E10, that lane changes with collission danger are blocked.<br>If the next state is KL, the function does NOT check or penalize for collission danger , because collission avoidance while staying in the lane is fully handled by the Trajectory module.|
|behavior___cost_speed| We penalize lanes that drive slower than the speed limit:<br><br> <b>WEIGHT_EFFICIENCY_SPEED_COST</b>\*<br>(SPEED_LIMIT-lane_speed)/SPEED_LIMIT<br><br> with:<br><b>WEIGHT_EFFICIENCY_SPEED_COST</b>=1.0|
|behavior___cost_target_lane|If all equal, we prefer to drive in the middle lane. This reflects driving behavior on our Michigan roads. It has the advantage that it is less likely to get caught in a pocket of slow traffic with inability to get out. The cost function for this is given by:<br><br> <b>WEIGHT_TARGET_LANE_COST</b>\*<br>(TARGET_LANE-next_lane)^2<br><br>This cost function must give an extremely small contribution to the overall cost, achieved by setting:<br><b>WEIGHT_TARGET_LANE_COST</b>=1.E-6|


(*) A logic was implemented to avoid that the car makes a rapid succession of lane changes. Too quick a succession of lane changes would lead to excessive jerk and must be avoided. It must also be avoided to ensure predictable driving behavior for the surrounding traffic. A simple but effective logic is used. Each time the state changes to LCL or LCR, indicating a lane change is starting, a counter is initialized with the number of trajectory points that must be consumed by the simulator before a new lane change will be permitted. The name of this counter is <b>lane_change_count_down</b>, and it is initialized to <b>150</b> points. Each cycle when the simulator returns the data, we calculate how many points of the trajectory it consumed, and subtract this from the counter. The boolean <b>stay_in_lane</b> will remain true until the counter has reached zero.


## 6. Description of Trajectory Module

The following functions take care of the Trajectory generation.

|Function| Description|
|-|-|
|trajectory___next_lane| Determines the next lane the car will move too (or stay in), using the maneuver recommended by the Bavior module|
|trajectory___ref_vel|Determines the new safe reference velocity that car can drive at, in the next lane.<br> Strict collision avoidance with the closest object in front of the car in the next lane is enforced.<br>If we are far enough (<b>GAP_TO_OBJECT_AHEAD_FOLLOW</b>=30m), then we will follow it with the object's speed, but if we are too close (<b>GAP_TO_OBJECT_AHEAD_BREAK</b>=15m) we will brake to increase the distance.|
|trajectory___update| Defines the new trajectory:<br>We re-use that part of the trajectory from the previous cycle that the simulator has not yet consumed, and extend it with new points to reflect the maneuver recommended by the Behavior module. We always define a trajectory with a fixed number of points (<b>TRAJ_NPOINTS</b>=50). <br>The points we add to the trajectory, typically between 4-10, are defined along a spline that is defined in the car's local coordinate system and 'merges' with the next lane's center line at 30m distance ahead of the previous path's end. The spline also goes through the next lanes's center line at 60m and 90m distance.<br>The points we add to the remainder of the previous path are spaced with a distance that reflects the new reference velocity.<br> This approach works really nice to avoid jerk and completely eliminates the need to use Polynial Trajectory Generation (PTG).


## 7. Result and Summary


The rubric points are fulfilled:
- The car drives more than 4.32 miles without incident. The video recording provided with the project submission was interupted at 10 miles, and in other tests the car drove over 30 miles without incident.
- The car does not go over the speed limit of 50 mph.
- The car goes around slower traffic whenever it is safe to do and maintains a speed close to the speed limit.
- Max Acceleration and Jerk are not exceeded.
- The car does not have collissions
- The car stays in its lane, <b>except that I have given it a preference to drive in the middle lane</b>!
- The car makes lane changes that are smooth and are completed before new lane changes are considered.


A part of the logic that could be further improved is that the car does not look backwards far enough, and relies a little bit on the objects behind it to adjust their speed. If a lane change is done at very low speeds, and objects from behind are coming at high speed, they sometimes come quite close. No collission was observed though.


---

### Dependencies
* You can download the Term3 Simulator which contains the Path Planning Project from the releases tab (https://github.com/udacity/self-driving-car-sim/releases)

* cmake >= 3.5
 * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets)
  * Run either `install-mac.sh` or `install-ubuntu.sh`.
  * If you install from source, checkout to commit `e94b6e1`, i.e.
    ```
    git clone https://github.com/uWebSockets/uWebSockets 
    cd uWebSockets
    git checkout e94b6e1
    ```

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./path_planning`.
