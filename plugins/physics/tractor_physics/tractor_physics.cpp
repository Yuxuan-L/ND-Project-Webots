/*
 * File: tractor_physics.cpp
 * Date: 2019.8.1
 * Description: physics plugin for Webots simulation
 * Author: Yonglin Jing
 * E-mail: 11712605@mail.sustech.edu.cn
 */

#include <ode/ode.h>
#include <plugins/physics.h>
#include <math.h>

//Green - x
//Red   - y
//Blue  - z

//Using Wismer Luth model
typedef struct {
  double slip_ratio;
  double wheel_numeric;
  double towed_force;
  double torque;
  double drawbar_pull;
  double trct_effc;   //tractive efficiency
}mechanics_t;

typedef struct{
  double theta;       //dependent on load 100N
  double angular_vel; //angular velocity
  double horiz_vel;   //horizontal velcity
  double wheel_w;     //wheel width
  double wheel_r;     //wheel radius
  double load;        //normal force form the ground

  dBodyID wheel_body;
  dGeomID wheel_geom;

  mechanics_t mechanics;
}wheel_t;

/*-------------Algrithm constant-------------*/
// Soil constant
#define DEFORMATION_COHESION_MODULUS   0.0f
#define DEFORMATION_FRICTIONAL_MODULUS 410400.0f  //4.104*10^5
#define SINKAGE_EXP 0.8f

// Vehivle specification
#define VEHICLE_LOAD 103f  //kg
#define FRONT_WHEEL_RADIUS 0.38
#define FRONT_WHEEL_WIDTH  0.19
#define REAR_WHEEL_RADIUS  0.6
#define REAR_WHEEL_WIDTH   0.37
#define THETA_CONSTANT 0.523598  // pi/6

static wheel_t front_right_w, front_left_w, rear_right_w, rear_left_w;
static double cone_index = 0.0;

//static pthread_mutex_t mutex; // needed to run with multi-threaded version of ODE

//Need to put Body and Geom into wheel_t, set the wheel instance into a vector
//static dBodyID frontRightWheelBody = NULL, frontLeftWheelBody = NULL, 
//                rearRightWheelBody = NULL, rearLeftWheelBody = NULL;

//static dGeomID frontRightWheelGeom = NULL, frontLeftWheelGeom = NULL,
//                rearRightWheelGeom = NULL, rearLeftWheelGeom  = NULL;

static dBodyID floorBody = NULL;
static dGeomID floorGeom = NULL;

static pthread_mutex_t mutex;

dVector3 f;

/*
 * Note: This plugin will become operational only after it was compiled and associated with the current world (.wbt).
 * To associate this plugin with the world follow these steps:
 *  1. In the Scene Tree, expand the "WorldInfo" node and select its "physics" field
 *  2. Then hit the [Select] button at the bottom of the Scene Tree
 *  3. In the list choose the name of this plugin (same as this file without the extention)
 *  4. Then save the .wbt by hitting the "Save" button in the toolbar of the 3D view
 *  5. Then reload the world: the plugin should now load and execute with the current simulation
 */

inline double cone_index_calc(double r,   //wheel radius  which??????????
                              double t,   //theta1, dependent on load, or use pre-calculated THETA_CONSTANT
                              double kc,  //cohesion modulus of deformation
                              double se,  //sinkage exponent
                              double kp){ //frictional modulus of deformation
  double zetam = 0.0;
  zetam = r*(1-cos(t));

  return 1.625*(kc/(se+1)*(pow((zetam+1.5),(se+1))-pow(zetam, (se+1)))+0.517*kp*(pow((zetam+1.5), (se+2))/(se+1)/(se+2)+pow(zetam, (se+2))/(se+2)-(zetam+1.5)*pow(zetam, (se+1))/(se+1)));
}

inline double slip_ratio_calc(double r,  //wheel radius 
                              double w,  //wheel angular velocity
                              double vx){  //wheel horizontal velocity
  return (r*w-vx)/(r*w);
}

inline double wheel_numeric_calc(double b,  //wheel width 
                                 double r,  //wheel radius
                                 double L){ //wheel load
  return cone_index*b*2*r/L;
}

inline double towed_force_calc(double L,   //wheel load
                               double wn){ //wheel numeric
  return L*(1.2/wn+0.04);
}

inline double torque_calc(double r,   //wheel radius
                          double L,   //wheel load
                          double wn,  //wheel numeric
                          double sr){ //slip ratio
  
  return r*L*0.75*(1-exp(-0.3*wn*sr));
}

inline double drawbar_pull_calc(double t,   //wheel torque
                                double r,   //wheel radius
                                double tf){ //towed force
  return t/r-tf;
}

inline double tractive_efficiency_calc(double dp,  //drawbar pull
                                       double vx,  //horizontal velocity
                                       double t,   //torque
                                       double w){  //angular velocity
  return dp*vx/t/w;
}

void parameter_init() {
  front_right_w.wheel_r = FRONT_WHEEL_RADIUS;
  front_right_w.wheel_w = FRONT_WHEEL_WIDTH;

  front_left_w.wheel_r  = FRONT_WHEEL_RADIUS;
  front_left_w.wheel_w  = FRONT_WHEEL_WIDTH;
  
  rear_right_w.wheel_r  = REAR_WHEEL_RADIUS;
  rear_right_w.wheel_w  = REAR_WHEEL_WIDTH;
  
  rear_left_w.wheel_r   = REAR_WHEEL_RADIUS;
  rear_left_w.wheel_w   = REAR_WHEEL_WIDTH;

  cone_index = cone_index_calc();
}

void wheel_mechanics_calc(wheel_t* wheel){
  //Note that claculation order should not be changed

  wheel->mechanics.slip_ratio = slip_ratio_calc(wheel->wheel_r, 
                                                wheel->wheel_w, 
                                                wheel->horiz_vel);

  wheel->mechanics.wheel_numeric = wheel_numeric_calc(wheel->wheel_w, 
                                                      wheel->wheel_r, 
                                                      wheel->load);

  wheel->mechanics.towed_force = towed_force_calc(wheel->load, 
                                                  wheel->mechanics.wheel_numeric);

  wheel->mechanics.torque = torque_calc(wheel->wheel_r, 
                                        wheel->load, 
                                        wheel->mechanics.wheel_numeric, 
                                        wheel->mechanics.slip_ratio);

  wheel->mechanics.drawbar_pull = drawbar_pull_calc(wheel->mechanics.torque, 
                                                    wheel->wheel_r, 
                                                    wheel->mechanics.towed_force);

  wheel->mechanics.trct_effc = tractive_efficiency_calc(wheel->mechanics.drawbar_pull, 
                                                        wheel->horiz_vel, 
                                                        wheel->mechanics.torque, 
                                                        wheel->angular_vel);
}

void webots_physics_init() {
  pthread_mutex_init(&mutex, NULL);
  /*
   * Get ODE object from the .wbt model, e.g.
   *   dBodyID body1 = dWebotsGetBodyFromDEF("MY_ROBOT");
   *   dBodyID body2 = dWebotsGetBodyFromDEF("MY_SOLID");
   *   dGeomID geom2 = dWebotsGetGeomFromDEF("MY_SOLID");
   * If an object is not found in the .wbt world, the function returns NULL.
   * Your code should correcly handle the NULL cases because otherwise a segmentation fault will crash Webots.
   *
   * This function is also often used to add joints to the simulation, e.g.
   *   dWorldID world = dBodyGetWorld(body1);
   *   pthread_mutex_lock(&mutex);
   *   dJointID joint = dJointCreateBall(world, 0);
   *   dJointAttach(joint, body1, body2);
   *   pthread_mutex_unlock(&mutex);
   *   ...
   */

  floorBody = dWebotsGetBodyFromDEF("FLOOR");
  floorGeom = dWebotsGetGeomFromDEF("FLOOR");
  
  if(floorBody == NULL){
    dWebotsConsolePrintf("!!! error : could not get floorBody.\r\n");
  }
  else if(floorGeom == NULL){
    dWebotsConsolePrintf("!!! error : could not get floorGeom.\r\n");
  }else {
    dWebotsConsolePrintf("floor Body&Geom found !\r\n");
  }

  /*
   * TFRW - Tractor Front Right Wheel
   * TFLW - Tractor Front Left  Wheel
   * TRRW - Tractor Rear  Right Wheel
   * TRLW - Tractor Rear  Left  Wheel
   */

  front_right_w.wheel_body = dWebotsGetBodyFromDEF("TFRW");
  front_right_w.wheel_geom = dWebotsGetGeomFromDEF("TFRW");

  if (front_right_w.wheel_body == NULL){
    dWebotsConsolePrintf("!!! error : could not get frontRightWheelBody.\r\n");
  }else if (front_right_w.wheel_geom == NULL){
    dWebotsConsolePrintf("!!! error : could not get frontRightWheelGeom.\r\n");
  }else {
    dWebotsConsolePrintf("frontRightWheel Body&Geom found !\r\n");
  }

  front_left_w.wheel_body = dWebotsGetBodyFromDEF("TFLW");
  front_left_w.wheel_geom = dWebotsGetGeomFromDEF("TFLW");

  if (front_left_w.wheel_body == NULL){
    dWebotsConsolePrintf("!!! error : could not get frontLeftWheelBody.\r\n");
  }else if(front_left_w.wheel_geom == NULL){
    dWebotsConsolePrintf("!!! error : could not get frontLeftWheelGeom.\r\n");
  }else {
    dWebotsConsolePrintf("frontLeftWheel Body&Geom found !\r\n");
  }
  
  rear_right_w.wheel_body = dWebotsGetBodyFromDEF("TRRW");
  rear_right_w.wheel_geom = dWebotsGetGeomFromDEF("TRRW");

  if (rear_right_w.wheel_body == NULL){
    dWebotsConsolePrintf("!!! error : could not get rearRightWheelBody.\r\n");
  }else if(rear_right_w.wheel_geom == NULL){
    dWebotsConsolePrintf("!!! error : could not get rearRightWheelGeom.\r\n");
  }else {
    dWebotsConsolePrintf("rearRightWheel Body&Geom found !\r\n");
  }
  
  rear_left_w.wheel_body = dWebotsGetBodyFromDEF("TRLW");
  rear_left_w.wheel_geom = dWebotsGetGeomFromDEF("TRLW");

  if (rear_left_w.wheel_body == NULL){
    dWebotsConsolePrintf("!!! error : could not get rearLeftWheelBody.\r\n");
  }else if(rear_left_w.wheel_geom == NULL){
    dWebotsConsolePrintf("!!! error : could not get rearLeftWheelGeom.\r\n");
  }else {
    dWebotsConsolePrintf("rearLeftWheel Body&Geom found !\r\n");
  }

}

void webots_physics_step() {
  /*
   * Do here what needs to be done at every time step, e.g. add forces to bodies
   *   dBodyAddForce(body1, f[0], f[1], f[2]);
   *   ...
   */
  pthread_mutex_lock(&mutex);
  //const dReal* tmp = dBodyGetForce(frontRightWheelBody);
  //dWebotsConsolePrintf("get load: %f %f %f\r\n", tmp[0], tmp[1], tmp[2]);
  pthread_mutex_unlock(&mutex);


  //dBodyAddRelForce(frontRightWheelBody, f[0], f[1], f[2]);
  //dBodyAddRelForce(frontLeftWheelBody, f[0], f[1], f[2]);
}

void webots_physics_step_end(){
  pthread_mutex_lock(&mutex);
  //const dReal* tmp = dBodyGetForce(frontRightWheelBody);
  //dWebotsConsolePrintf("get load: %f %f %f\r\n", tmp[0], tmp[1], tmp[2]);
  pthread_mutex_unlock(&mutex);
}

int webots_physics_collide(dGeomID g1, dGeomID g2) {
  /*
   * This function needs to be implemented if you want to overide Webots collision detection.
   * It must return 1 if the collision was handled and 0 otherwise.
   * Note that contact joints should be added to the contact_joint_group which can change over the time, e.g.
   *   n = dCollide(g1, g2, MAX_CONTACTS, &contact[0].geom, sizeof(dContact));
   *   dJointGroupID contact_joint_group = dWebotsGetContactJointGroup();
   *   dWorldID world = dBodyGetWorld(body1);
   *   ...
   *   pthread_mutex_lock(&mutex);
   *   dJointCreateContact(world, contact_joint_group, &contact[i])
   *   dJointAttach(contact_joint, body1, body2);
   *   pthread_mutex_unlock(&mutex);
   *   ...
   */
  //dContact contact[4];
  
  //if (dAreGeomsSame(g1, frontRightWheelGeom) || dAreGeomsSame(g2, frontRightWheelGeom)){
  //  if((dAreGeomsSame(g1,frontRightWheelGeom) && dAreGeomsSame(g2,floorGeom)) || (dAreGeomsSame(g1,floorGeom) && dAreGeomsSame(g2,frontRightWheelGeom))){
      
  //  }
  //}

  return 0;
}

void webots_physics_cleanup() {
  /*
   * Here you need to free any memory you allocated in above, close files, etc.
   * You do not need to free any ODE object, they will be freed by Webots.
   */
  pthread_mutex_destroy(&mutex);
}
