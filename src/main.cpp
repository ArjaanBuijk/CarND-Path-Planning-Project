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

// Highway has 6 lanes- 3 heading in each direction.
// Each lane is 4 m wide.
const double lane_width = 4.0;
const vector<double> d_lane = {0.5*lane_width, 1.5*lane_width, 2.5*lane_width}; // d of left, middle, right lane center line

// Typical car width & length in USA
// See: https://www.reference.com/vehicles/width-length-average-car-9eb7b00283fb1bd8
const double car_width  = inch2m(75.0);
const double car_length = inch2m(200.0);

// Speedlimit
const double speed_limit = mph2mps(50);

// Target velocity a little bit below the speedlimit
const double speed_target = mph2mps(49.5);

// Flag danger when the gap with car in front is too close
const double gap_too_close = 30.0;

// Change in ref_vel per cycle to use when it can accellerate or when it needs to slow down
// This value is chosen to stay within allowed longitudinal jerk
const double ideal_dvel = 0.224;

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


// Helper function to check if a car is in a certain lane
bool car_in_lane(double d,int lane){

  if ( d+0.5*car_width > d_lane[lane]-0.5*lane_width &&
       d-0.5*car_width < d_lane[lane]+0.5*lane_width ){
    return true;
  }
  else{
    return false;
  }

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
  int lane = 1;

  // start at 0.0 velocity, and let logic apply acceleration to target velocity
  double ref_vel = 0.0;

  h.onMessage([&lane,&ref_vel,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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

        	// Main car's localization Data
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
          double car_d = j[1]["d"];
          double car_yaw = j[1]["yaw"];
          double car_speed = mph2mps(j[1]["speed"]);

          cout<<"-----------------------------------------------------\n"
              <<"Car's localization data provided by simulator:\n"
              <<"car_x    = "<<car_x<<'\n'
              <<"car_y    = "<<car_y<<'\n'
              <<"car_s    = "<<car_s<<'\n'
              <<"car_d    = "<<car_d<<'\n'
              <<"car_yaw  = "<<car_yaw<<'\n'
              <<"car_speed= "<<car_speed<<'\n';

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

          int prev_size = previous_path_x.size();

          json msgJson;

          vector<double> next_x_vals;
          vector<double> next_y_vals;

          // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds

          // test1: Move car forward in straight line at 50 mph
          /*
          double dist_inc = 0.5;
          for(int i = 0; i < 50; i++)
          {
            next_x_vals.push_back(car_x+(dist_inc*i)*cos(deg2rad(car_yaw)));
            next_y_vals.push_back(car_y+(dist_inc*i)*sin(deg2rad(car_yaw)));
          }
          */

          // test2: Move car at 50 mph, while staying in lane, and drive smoothly
          //        and avoids collisions with other cars in front of it

          //******************************************************************************************************
          // Check on other cars in front
          if (prev_size > 0){
            car_s = end_path_s;
          }

          bool too_close = false;
          int index_closest=-1;
          double gap_closest = 1000.0; // we look for the closest car that is in our lane, in front of us
          double speed_closest = 0.0;

          for (size_t i=0; i<sensor_fusion.size(); ++i){
            float d = sensor_fusion[i][6];

            if (car_in_lane(d,lane)){
              double vx = sensor_fusion[i][3];
              double vy = sensor_fusion[i][4];
              double check_speed = sqrt(vx*vx+vy*vy);
              double check_car_s = sensor_fusion[i][5];
              // Calculate where the car will be at end time of our previous path
              check_car_s += ((double)prev_size*0.02*check_speed);

              if ( check_car_s > car_s) {       // is it in front of us?
                double gap = check_car_s - car_s;
                if (gap < gap_closest){
                  index_closest = i;
                  gap_closest = gap;
                  if ( gap < gap_too_close ){
                    cout<<"Found a car in front, in our lane that is too close!\n";
                    too_close = true;
                    speed_closest = check_speed;
                    // TO DO: add further logic, like:
                    // (-) flag it for slow-down
                    // (-) flag it for potential passing
                    if (lane > 0){
                        lane -= 1; // blindly pass on the left...
                    }
                  }
                }
              }
            }
          }

          // set new reference velocity
          if (too_close){
            if (ref_vel < speed_closest){
              ref_vel += ideal_dvel;
            }
            else if (ref_vel > speed_closest){
              ref_vel -= ideal_dvel;
            }
          }
          else {
            if (ref_vel < speed_target){
              ref_vel += ideal_dvel;
            }
            else if (ref_vel > speed_target){
              ref_vel -= ideal_dvel;
            }
          }

          //******************************************************************************************************
          // Create a list of widely spaced (x,y) waypoints, evenly spaced at 30m
          // later we will interpolate these waypoints with a spline and fill it in with more points
          vector<double> ptsx;
          vector<double> ptsy;

          // reference x,y, yaw states
          // either we will reference the starting point as where the car is or at the previous paths end point
          double ref_x = car_x;
          double ref_y = car_y;
          double ref_yaw = deg2rad(car_yaw);

          // if previous size is almost empty, use the car as starting reference
          if (prev_size < 2) {
            // Use two points that make the path tangent to the car
            double prev_car_x = car_x - cos(deg2rad(car_yaw));
            double prev_car_y = car_y - sin(deg2rad(car_yaw));

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
          vector<double> next_wp0 = getXY(car_s+30,d_lane[lane],map_waypoints_s,map_waypoints_x, map_waypoints_y);
          vector<double> next_wp1 = getXY(car_s+60,d_lane[lane],map_waypoints_s,map_waypoints_x, map_waypoints_y);
          vector<double> next_wp2 = getXY(car_s+90,d_lane[lane],map_waypoints_s,map_waypoints_x, map_waypoints_y);

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

          // Add more points, driving at new ref_vel, until we have a path of 30m long
          cout<<"ref_vel = "<<ref_vel<<'\n';
          cout<<"x_end_l = "<<x_end_l<<'\n';
          // Here we will always output 50 points.
          for (int i=0; i< 50-prev_size; ++i) {
            x_end_l += 0.02*ref_vel;
            y_end_l = s(x_end_l);

            cout<<"x_end_l = "<<x_end_l<<'\n';

            // Transform it back to global coordinates from the car local coordinates
            double x_end = ref_x + x_end_l*cos(ref_yaw) - y_end_l*sin(ref_yaw);
            double y_end = ref_y + x_end_l*sin(ref_yaw) + y_end_l*cos(ref_yaw);

            next_x_vals.push_back(x_end);
            next_y_vals.push_back(y_end);
          }


          /*
          // Calculate how to break up spline points so that we travel at our desired reference velocity
          double target_x = 30.0; // going out 30 m to the horizon
          double target_y = s(target_x);
          double target_dist = sqrt((target_x)*(target_x) + (target_y)*(target_y));
          double N = (target_dist/(0.02*ref_vel)); // Number of points we need to put on the spline to travel with target velocity
          double dx = (target_x)/N; // x-distance between points


          // Fill up the rest of our pth planner after filling it with previous points.
          // Here we will always output 50 points.
          for (int i=0; i< 50-prev_size; ++i){
            double x_l = (i+1)*dx;
            double y_l = s(x_l);

            // Transform it back to global coordinates from the car local coordinates
            double x_point = ref_x + x_l*cos(ref_yaw) - y_l*sin(ref_yaw);
            double y_point = ref_y + x_l*sin(ref_yaw) + y_l*cos(ref_yaw);

            next_x_vals.push_back(x_point);
            next_y_vals.push_back(y_point);

          }
          */

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
















































































