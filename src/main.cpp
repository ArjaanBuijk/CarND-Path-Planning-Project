#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// For converting back and forth between mph and m/s
double mph2mps(double x) { return x * 0.44704; }
double mps2mph(double x) { return x / 0.44704; }

// For converting back and forth between inch and m
double inch2m(double x) { return x * 0.0254; }
double m2inch(double x) { return x / 0.0254; }

// Provide print-out on costs
const bool VERBOSE = true;

// Highway has 6 lanes- 3 heading in each direction.
// Each lane is 4 m wide.
const double LANES_AVAILABLE = 3;
const double LANE_WIDTH = 4.0;
const vector<double> D_LANES = {0.5*LANE_WIDTH, 1.5*LANE_WIDTH, 2.5*LANE_WIDTH}; // d of left, middle, right lane center line

// Use typical width & length in USA for both car and objects
// See: https://www.reference.com/vehicles/width-length-average-car-9eb7b00283fb1bd8
const double CAR_WIDTH     = inch2m(75.0);
const double CAR_LENGTH    = inch2m(200.0);
const double OBJECT_WIDTH  = inch2m(75.0);
const double OBJECT_LENGTH = inch2m(200.0);

// Speedlimit
const double SPEED_LIMIT = mph2mps(50.0);

// Target velocity a little bit below the speedlimit
const double SPEED_TARGET = mph2mps(49.5);

// Flag danger when the gap with object in front is too close
const double GAP_TO_OBJECT_AHEAD_BREAK  = 15.0; // if gap is less then this, break!
const double GAP_TO_OBJECT_AHEAD_FOLLOW = 30.0; // if gap is less than this, but higher than break gap, just follow at same speed

// Change in ref_vel per cycle to use when it can accellerate or when it needs to slow down
// This value is chosen to stay within allowed longitudinal jerk
//const double MAX_VEL_CHANGE = 0.224;
const double MAX_VEL_CHANGE = 0.2;

// Characteristics of the trajectory that we give back to the simulator
const double TRAJ_DT       = 0.02; // time between two points of trajectory
const int    TRAJ_NPOINTS  = 50;   // number of points in trajectory
const double TRAJ_DURATION = TRAJ_DT*TRAJ_NPOINTS; // duration of trajectory
const int    TRAJ_NPOINTS_TO_EAT_BETWEEN_LANE_CHANGES  = 150;   // number of trajectory points to eat between lane changes

// In our Finite State Machine, we consider 3 states:
//  KL  = Keep Lane
//  LCL = Lane Change Left
//  LCR = Lane Change Right
enum State {KL, LCL, LCR};

// Weights for the cost functions
const double WEIGHT_SAFETY_COLLISSION_COST     = 1.E10;
const double WEIGHT_EFFICIENCY_SPEED_COST      = 1.0;
const double WEIGHT_TARGET_LANE_COST           = 1.E-6;

// Regions behind and ahead for collision cost calculation
const double COLLISSION_COST_REGION_BEHIND = 10; //m
const double COLLISSION_COST_REGION_AHEAD  = 5; //m beyond end of previous path

// Regions behind and ahead for speed cost calculation
const double SPEED_COST_REGION_BEHIND =  0; //m
const double SPEED_COST_REGION_AHEAD  = 50; //m

// We like to drive in the middle lane, if all else is equal
const int TARGET_LANE = 1;



// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(size_t i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

// Use this one if ClosestWaypoint is behind the car
int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
// maps_x, maps_y are the map waypoints that are determined at the start of the simulation
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
// maps_s is calculated at beginning.
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}


// Helper function to check if an object is in a certain lane
bool is_object_in_lane(double d, int lane){
  if ( d+0.5*OBJECT_WIDTH > D_LANES[lane]-0.5*LANE_WIDTH &&
       d-0.5*OBJECT_WIDTH < D_LANES[lane]+0.5*LANE_WIDTH ){
    return true;
  }
  else{
    return false;
  }
}

// Helper function to get the lane of an object
int lane_of_object(double d){
  int lane=0;
  for (int i=0; i<LANES_AVAILABLE; i++){
    if ( d <= (i+1)*LANE_WIDTH)
      return lane;
    lane++;
  }
  return LANES_AVAILABLE; // when we get here, it appears object is off the road...
}

// Possible successor states for car
vector<int> possible_successor_states(int lane, int state, int lane_change_count_down) {
  vector<int> next_states;
  if (lane_change_count_down>0){
    next_states.push_back(KL);           // first finish a lane change maneuver
  }
  else {
    next_states.push_back(state);         // each state supports self transition
    if (state==KL){
      if (lane==0){
        next_states.push_back(LCR);
      }
      else if (lane==LANES_AVAILABLE - 1){
        next_states.push_back(LCL);
      }
      else{
        next_states.push_back(LCL);
        next_states.push_back(LCR);
      }
    }
    else {
      next_states.push_back(KL);
    }
  }
  return next_states;
}

// Safety collission danger cost:
// We penalize lane changes have collission danger
double safety_collission_danger_cost(int lane, int state, int next_state,
                   double car_x, double car_y, double car_s, double end_path_s,
                   vector<double> previous_path_x, vector<double> previous_path_y,
                   vector<vector<double>> sensor_fusion, vector<vector<double>> predictions,
                   vector<double> lane_speeds){

  double cost = 0.0;

  int next_lane = lane;
  if (next_state != KL){
    if (next_state == LCL) {
      next_lane -= 1;
    }
    if (next_state == LCR) {
      next_lane += 1;
    }

    // when switching lane, is there an object on it's path ?
    bool too_close = false;
    for (size_t i=0; i<sensor_fusion.size(); ++i){
      double object_s   = sensor_fusion[i][5];
      double object_d   = sensor_fusion[i][6];

      if (is_object_in_lane(object_d,next_lane)){
        if ( object_s > car_s-COLLISSION_COST_REGION_BEHIND &&
             object_s < end_path_s+COLLISSION_COST_REGION_AHEAD) {
          too_close = true;
          break;
        }
      }
    }
    if (too_close) {
      cost += WEIGHT_SAFETY_COLLISSION_COST;
    }
  }

  if (VERBOSE)
    cout<<"lane, next lane, collision cost = "<<lane<<", "<<next_lane<<", "<<cost<<'\n';
  return cost;
}

// Efficiency Speed cost:
// We penalize lanes that drive slower than the speedlimit
double efficiency_speed_cost(int lane, int state, int next_state,
               vector<double> lane_speeds){

  int next_lane = lane;
  if (next_state == LCL) {
    next_lane -= 1;
  }
  if (next_state == LCR) {
    next_lane += 1;
  }

  double cost = 0.0;
  if (lane_speeds[next_lane] <= SPEED_LIMIT) {
    cost += WEIGHT_EFFICIENCY_SPEED_COST*(SPEED_LIMIT - lane_speeds[next_lane])/SPEED_LIMIT;
  }

  if (VERBOSE){
    cout<<"lane speed in lane              = "<<lane_speeds[lane]<<'\n';
    cout<<"lane speed in next lane         = "<<lane_speeds[next_lane]<<'\n';
    cout<<"lane, next lane, speed cost     = "<<lane<<", "<<next_lane<<", "<<cost<<'\n';
  }
  return cost;
}

// Target lane cost:
// We penalize lanes that are not the preferred target lane
double cost_target_lane(int lane, int state, int next_state ){

  int next_lane = lane;
  if (next_state == LCL) {
    next_lane -= 1;
  }
  if (next_state == LCR) {
    next_lane += 1;
  }

  double cost = 0.0;
  cost += WEIGHT_TARGET_LANE_COST*(TARGET_LANE-next_lane)*(TARGET_LANE-next_lane);

  if (VERBOSE){
    cout<<"lane, next lane, target cost    = "<<lane<<", "<<next_lane<<", "<<cost<<'\n';
  }
  return cost;
}

// Calculate total cost to transition to next_state
double calculate_cost(int lane, int state, int next_state,
                   double car_x, double car_y, double car_s, double end_path_s,
                   vector<double> previous_path_x, vector<double> previous_path_y,
                   vector<vector<double>> sensor_fusion, vector<vector<double>> predictions,
                   vector<double> lane_speeds){
  double cost = 0.0;

  cost+= safety_collission_danger_cost( lane, state, next_state,
                                      car_x, car_y, car_s, end_path_s,
                                      previous_path_x, previous_path_y,
                                      sensor_fusion, predictions,
                                      lane_speeds);

  cost+= efficiency_speed_cost( lane,  state,  next_state,
                     lane_speeds);


  cost+= cost_target_lane( lane,  state,  next_state );

  return cost;
}




int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  const double max_s = 6945.554; // This is length of the track in meters.


  // The d vector is pointing in the direction of the right-hand side of the road.

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  // start in lane 1
  int prev_lane = 1;
  int lane = 1;

  // finish lane change once started
  int lane_change_count_down = 0;

  // start at 0.0 velocity, and let logic apply acceleration to target velocity
  double ref_vel = 0.0;

  // start in state KL
  int state = KL;


  h.onMessage([&prev_lane,&lane,&lane_change_count_down,&ref_vel,&state,&max_s,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object
          if (VERBOSE)
            cout<<"------------------------------------------------------------------------\n";

        	// car's localization Data, converted to SI units
          double car_x = j[1]["x"]; // m
          double car_y = j[1]["y"]; // m
          double car_s = j[1]["s"]; // m
          double car_d = j[1]["d"]; // m
          double car_yaw = deg2rad(j[1]["yaw"]); // rad
          //double car_speed = mph2mps(j[1]["speed"]); // m/s

          // calculate the lane the car is in now
          //lane      = lane_of_object(car_d);

          // Remainder of previous trajectory that car has not yet traversed
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          int prev_size = previous_path_x.size();

          // Frenet s and d values at end of previous trajectory
          double end_path_s = j[1]["end_path_s"];
          //double end_path_d = j[1]["end_path_d"];
          if (prev_size==0){
            end_path_s = car_s;
            //end_path_d = car_d;
          }
          else if (lane_change_count_down>0){
            lane_change_count_down -= (TRAJ_NPOINTS - prev_size);
          }
          if (VERBOSE)
            cout<<"lane_change_count_down = "<<lane_change_count_down<<'\n';

          // Sensor Fusion Data, a list of all objects on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

          json msgJson;

          vector<double> next_x_vals;
          vector<double> next_y_vals;

          // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds

          //******************************************************************************************************
          // PREDICTION - Start

          // simple prediction of object location at end of next trajectory
          // -> assume they move with constant velocity

          vector<vector<double>> predictions;
          for (size_t i=0; i<sensor_fusion.size(); ++i){
            //int    object_id  = sensor_fusion[i][0];
            double object_x   = sensor_fusion[i][1];
            double object_y   = sensor_fusion[i][2];
            double object_vx  = sensor_fusion[i][3];
            double object_vy  = sensor_fusion[i][4];
            double object_s   = sensor_fusion[i][5];
            double object_d   = sensor_fusion[i][6];

            double object_speed = sqrt(object_vx*object_vx + object_vy*object_vy);

            object_x += object_vx*TRAJ_DURATION;
            object_y += object_vy*TRAJ_DURATION;
            object_s += object_speed*TRAJ_DURATION;

            vector<double> prediction;
            prediction.push_back(object_x );
            prediction.push_back(object_y );
            prediction.push_back(object_s );
            prediction.push_back(object_d );
            prediction.push_back(object_speed );

            predictions.push_back(prediction);
          }

          // PREDICTION - End
          //******************************************************************************************************


          //******************************************************************************************************
          // BEHAVIOR - Start

          // determine possible successor states from current state
          vector<int> next_states;
          next_states = possible_successor_states(lane, state, lane_change_count_down);

          // for all the lanes:
          // - calculate the miniumum speed of objects in a region around the car
          // - calculate the open space in front of the car
          // right now, and at end of previous trajectory.
          vector<double> lane_speeds    = {1000.0, 1000.0, 1000.0};
          /*
          for (size_t i=0; i<sensor_fusion.size(); ++i){
            double object_vx  = sensor_fusion[i][3];
            double object_vy  = sensor_fusion[i][4];
            double object_s   = sensor_fusion[i][5];
            double object_d   = sensor_fusion[i][6];
            double object_speed = sqrt(object_vx*object_vx + object_vy*object_vy);
            int object_lane = lane_of_object(object_d);

            //if (object_s > car_s      - SPEED_COST_REGION_BEHIND &&
            //    object_s < end_path_s + SPEED_COST_REGION_AHEAD) {
            if (object_s > car_s      - SPEED_COST_REGION_BEHIND &&
                object_s < end_path_s + SPEED_COST_REGION_AHEAD) {

              lane_speeds[object_lane] = min(lane_speeds[object_lane], object_speed);
            }
          }
          */

          // speed of closest car in front of us
          vector<double> gap_closests = {1000.0,1000.0,1000.0}; // we look for the closest car in each lane that is in front of us
          for (size_t i=0; i<sensor_fusion.size(); ++i){
            //int    object_id  = sensor_fusion[i][0];
            //double object_x   = sensor_fusion[i][1];
            //double object_y   = sensor_fusion[i][2];
            double object_vx  = sensor_fusion[i][3];
            double object_vy  = sensor_fusion[i][4];
            double object_s   = sensor_fusion[i][5];
            double object_d   = sensor_fusion[i][6];

            double object_speed = sqrt(object_vx*object_vx + object_vy*object_vy);

            int object_lane = lane_of_object(object_d);

            if (object_s > car_s &&
                object_s < car_s + SPEED_COST_REGION_AHEAD) { // if it is in front of us now and not too far away

              // Calculate where the object will be at end time of our previous path
              object_s += ((double)prev_size*TRAJ_DT*object_speed);

              if ( object_s > end_path_s &&
                   object_s < end_path_s + SPEED_COST_REGION_AHEAD) { // if it is in front of us at end and not too far away
                double gap = object_s - end_path_s;
                if (gap < gap_closests[object_lane]){
                  gap_closests[object_lane] = gap;
                  lane_speeds[object_lane] = object_speed;
                }
              }
            }
          }

          // calculate the costs for each possible successor state
          vector<double> costs;
          for (int next_state : next_states){
            double cost = calculate_cost(lane, state, next_state,
                                      car_x, car_y, car_s, end_path_s,
                                      previous_path_x, previous_path_y,
                                      sensor_fusion, predictions,
                                      lane_speeds);
            costs.push_back(cost);
          }


          // select successor state with lowest cost
          int index_best = min_element(costs.begin(), costs.end()) - costs.begin();
          int next_state = next_states[index_best];



          // BEHAVIOR - End
          //******************************************************************************************************


          //******************************************************************************************************
          // TRAJECTORY - Start


          // Check on objects in front of previous path's end

          bool close_break  = false;
          bool close_follow = false;
          double gap_closest = 1000.0; // we look for the closest car that is in our lane, in front of us
          double speed_closest = 0.0;

          for (size_t i=0; i<sensor_fusion.size(); ++i){
            //int    object_id  = sensor_fusion[i][0];
            //double object_x   = sensor_fusion[i][1];
            //double object_y   = sensor_fusion[i][2];
            double object_vx  = sensor_fusion[i][3];
            double object_vy  = sensor_fusion[i][4];
            double object_s   = sensor_fusion[i][5];
            double object_d   = sensor_fusion[i][6];

            if (is_object_in_lane(object_d,lane)){

              double object_speed = sqrt(object_vx*object_vx + object_vy*object_vy);

              // Calculate where the car will be at end time of our previous path
              object_s += ((double)prev_size*TRAJ_DT*object_speed);

              if ( object_s > end_path_s) {       // is it in front of us at end of previous trajectory ?
                double gap = object_s - end_path_s;
                if (gap < gap_closest){
                  gap_closest = gap;
                  if ( gap < GAP_TO_OBJECT_AHEAD_FOLLOW ){
                    close_follow  = true;
                    speed_closest = object_speed;
                    if ( gap < GAP_TO_OBJECT_AHEAD_BREAK ){
                      close_break = true;
                    }
                  }
                }
              }
            }
          }

          // Check if we're switching lanes, and start count_down
          if (next_state == LCL || next_state == LCR ){
            prev_lane = lane;
            lane_change_count_down = TRAJ_NPOINTS_TO_EAT_BETWEEN_LANE_CHANGES;
            if (VERBOSE)
              cout<<"SWITCHING LANE - lane_change_count_down = "<<lane_change_count_down<<'\n';
            if (next_state == LCL){
              lane -= 1;
            }
            else if (next_state == LCR) {
              lane += 1;
            }
          }

          // set new reference velocity
          if (close_break) {                      // really close to next car --> Break!
            ref_vel -= MAX_VEL_CHANGE;
          }
          else if (close_follow){                 // close, but not too close to next car --> Follow!
            if (ref_vel < speed_closest){
              ref_vel += MAX_VEL_CHANGE;
            }
            else if (ref_vel > speed_closest){
              ref_vel -= MAX_VEL_CHANGE;
            }
          }
          else {                                  // all clear --> Go speed limit!
            if (ref_vel < SPEED_TARGET){
              ref_vel += MAX_VEL_CHANGE;
            }
            else if (ref_vel > SPEED_TARGET){
              ref_vel -= MAX_VEL_CHANGE;
            }
          }


          // Create a list of widely spaced (x,y) waypoints, evenly spaced at 30m
          // later we will interpolate these waypoints with a spline and fill it in with more points
          vector<double> ptsx;
          vector<double> ptsy;

          // reference x,y, yaw states
          // either we will reference the starting point as where the car is or at the previous paths end point
          double ref_x = car_x;
          double ref_y = car_y;
          double ref_yaw = car_yaw;

          // if previous size is almost empty, use the car as starting reference
          if (prev_size < 2) {
            // Use two points that make the path tangent to the car
            double prev_car_x = car_x - cos(car_yaw);
            double prev_car_y = car_y - sin(car_yaw);

            ptsx.push_back(prev_car_x);
            ptsx.push_back(car_x);

            ptsy.push_back(prev_car_y);
            ptsy.push_back(car_y);
          }
          // use the previous path's end oint as starting reference
          else {
            // Redefine reference state as previous path end points
            ref_x = previous_path_x[prev_size-1];
            ref_y = previous_path_y[prev_size-1];

            double ref_x_prev = previous_path_x[prev_size-2];
            double ref_y_prev = previous_path_y[prev_size-2];
            ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

            // Use two points that make the path tangent to the previous path's end point
            ptsx.push_back(ref_x_prev);
            ptsx.push_back(ref_x);

            ptsy.push_back(ref_y_prev);
            ptsy.push_back(ref_y);
          }

          // In Frenet coordinates, add evenly 30m spaced points ahead of the starting reference
          vector<double> next_wp0 = getXY(end_path_s+30,D_LANES[lane],map_waypoints_s,map_waypoints_x, map_waypoints_y);
          vector<double> next_wp1 = getXY(end_path_s+60,D_LANES[lane],map_waypoints_s,map_waypoints_x, map_waypoints_y);
          vector<double> next_wp2 = getXY(end_path_s+90,D_LANES[lane],map_waypoints_s,map_waypoints_x, map_waypoints_y);

          ptsx.push_back(next_wp0[0]);
          ptsx.push_back(next_wp1[0]);
          ptsx.push_back(next_wp2[0]);

          ptsy.push_back(next_wp0[1]);
          ptsy.push_back(next_wp1[1]);
          ptsy.push_back(next_wp2[1]);

          // Transform x,y of points from global into local car coordinate system
          for (size_t i=0; i < ptsx.size(); ++i){
            double dx = ptsx[i] - ref_x;
            double dy = ptsy[i] - ref_y;

            ptsx[i] =  dx*cos(ref_yaw) + dy*sin(ref_yaw);
            ptsy[i] = -dx*sin(ref_yaw) + dy*cos(ref_yaw);
          }

          // Create a spline
          tk::spline s;

          // Set (x,y) points to the spline
          s.set_points(ptsx, ptsy);

          double x_end_l = 0.0;
          double y_end_l = 0.0;


          // Start with all of the remaining previous path points from last time, if we have any
          // Do not re-create them each time. This helps with the transitions.
          if (prev_size > 0){
            for (int i=0; i<prev_size; ++i){
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }
            double x_end = previous_path_x[prev_size-1];
            double y_end = previous_path_y[prev_size-1];
            // Transform last point into local car coordinate system
            double dx = x_end - ref_x;
            double dy = y_end - ref_y;
            x_end_l =  dx*cos(ref_yaw) + dy*sin(ref_yaw);
            y_end_l = -dx*sin(ref_yaw) + dy*cos(ref_yaw);
          }

          // Add new points to previous trajectory until we have TRAJ_NPOINTS points total
          // Here we will always output TRAJ_NPOINTS points.
          for (int i=0; i< TRAJ_NPOINTS-prev_size; ++i) {
            x_end_l += TRAJ_DT*ref_vel; // To drive at new ref_vel we need to space the points
                                     // along the trajectory with TRAJ_DT*ref_vel.
                                     // Because we work in the car local coordinate system,
                                     // it is ok to space the x value with this amount.
            y_end_l = s(x_end_l);

            // Transform it back to global coordinates from the car local coordinates
            double x_end = ref_x + x_end_l*cos(ref_yaw) - y_end_l*sin(ref_yaw);
            double y_end = ref_y + x_end_l*sin(ref_yaw) + y_end_l*cos(ref_yaw);

            next_x_vals.push_back(x_end);
            next_y_vals.push_back(y_end);
          }

          // TRAJECTORY - End
          //******************************************************************************************************

          // END OF TODO
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"control\","+ msgJson.dump()+"]";

          //this_thread::sleep_for(chrono::milliseconds(1000));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
