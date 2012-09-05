#include "WProgram.h"
void homingInterruptHandler();
void sensingInterruptHandler();
int getOpenLoopStepSetpoint(float omega);
void takeStep(int dir, int nStep);
void updateSetpoint();
void updateAngularVelocity();
void updateError();
void setup();
void loop();
unsigned long cnt = 0;
/////////////// pin declarations ////////////////////////
const int limitPin = 2;   //this pin goes low when the pot is limited to the left
const int setPin = 0;     // select the input pin the setpoint
const int sensorInterruptPin = 3;
const int ledPin = 13;    // select the pin for the LED
const int stepperPins[2] = {
  4, 5}; //stepper for pot controlled on these pins
///////////////////////////////////////////////////////////

/////////////// declare data array, its length, and indexing variable//////////////
#define bufferLen  10 		//number of points in sensing input buffer
volatile unsigned long T[bufferLen] ; 	// sensing input buffer of periods in microseconds
volatile unsigned long t[bufferLen] ; 	// sensing input timestamp buffer in microseconds
volatile unsigned long idx = 0;		//counter for T
volatile unsigned long now, then; //hold clock values between sensing interrupts
////////////////////////////////////////////////////////////////////////////////////

////////////////////// boolean state variables ////////////////////////
boolean bufferFilled = false;		//sensor buffer full of valid readings?
volatile boolean homed = false;			//has the pot been homed?
boolean setpointChanged = true;	        //has the setpoint changed since last update?
const boolean limitState = LOW;         //state of limitPin when the limit switch is tripped
//////////////////////////////////////////////////////////////////////////

//\\\\\\\\\\\\\\\\\\\\\\\\\\\ control state variables \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

/////////stepping state, setpoint, and parameters////////////
const boolean CW = LOW; //write CW to stepperPins[0] to drive clockwise
volatile int n = 0;				//number of steps CW from the home position
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! calibrate!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
const int maxStep = 1000;  		//max number of steps we can take from the home position
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
int openLoopStepSetpoint = 0;		//open loop setpoint for the step
boolean DIR = CW;				//temp variable to hold direction for stepping
const unsigned long stepDelayLoopCount = 10000;         //nr of loops to burn indirectly sets the length of time to wait for the motor to step
///////////////////////////////////////////////////////

///////////angular velocity state, setpoint, and parameters//////////////
float angularVelocity = 0;				//current angular velocity computed from the measured periods
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! calibrate!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
const float angVelLim[2] = {
  70, 150}; 	                // angular velocity at {n=0, n = maxStep}
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
float angularVelocitySetpoint = 100;	        //initialize setpoint to middle speed
float angVelSetpointRez;				//resolution of setpoint in rad/s
float setpointChangeThres = 10;		//consider angular velocity setpoint changed if the change exceeds this threshold
const bool constantSetpoint = true; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!never update the setpoint, must hardcode desired setpoint into setup()!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
///////////////////////////////////////////////////////////////////


///////error buffer, derivative, integral, and proportional components////////////
#define errLen  3			//number of points in the error buffer 
float E[errLen];				//buffer of errors in angular velocity
float Ei=0, Ep=0, Ed= 0;  			//integral, error, and derivative
const float iThres = 5, pThres = 10, dThres=1000;	//thresholds for setpoint error triggering integral, proportional, and derivative action, respectively
float timestamps[errLen];                        //timestamp the samples in E
//////////////////////////////////////////////////////////////////////////

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


//:::::::::::::::::::::: INTERRUPTS:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::	
void homingInterruptHandler(){
  unsigned long dummyIdx,dumb=0;
  noInterrupts();
  // when the limitPin goes low, actuate this homing routine for the pot
  Serial.println("homing interrupt");
  for(dummyIdx=0;dummyIdx<stepDelayLoopCount;dummyIdx++ )
    dumb+=1;

  while(digitalRead(limitPin)==limitState){
    //Serial.println("taking homing steps");
    takeStep(CW,10);
    for(dummyIdx=0;dummyIdx<stepDelayLoopCount;dummyIdx++ )
      dumb+=1;

  }
  takeStep(CW,10);//go an additional 10 steps to make sure
  homed = true;
  Serial.println(homed,BIN);
  n = 0;
  Serial.println("end homing interrupt");
  interrupts();
}
void sensingInterruptHandler(){
  // when the reflectance goes high, record the time since the last spike
  noInterrupts();
  now = micros();
  Serial.print("sensing interrupt, now - then = ");
  Serial.println(now-then);

  if(now>then){
    T[idx%bufferLen] = now-then;
    t[idx%bufferLen] = now;
    idx++;
  }
  then = now;	
  interrupts();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::						

//-------------------------------helper functions-----------------------------

int getOpenLoopStepSetpoint(float omega){
  //map the angular velocity omega to an open loop step number nOpen using 
  //linear interpolation between angVelLim and (0, stepMax)
  return maxStep*(omega - angVelLim[0])/(angVelLim[1]-angVelLim[0]);
}

void takeStep(int dir, int nStep){
  int j, dumb;
  //drive the pot nSteps in the dir direction
  // dir = -1 : step CCW
  // dir =  1 : step CW
  //Serial.println("taking step"); 
  //set direction ///
  if(dir)
    digitalWrite(stepperPins[0],HIGH);
  else
    digitalWrite(stepperPins[0],LOW);
  //////////////////

  ///take the specified number of steps, incrementing n as necessary///
  digitalWrite(stepperPins[1],LOW);
  for(int i = 0; i<nStep;i++){
     //pulse the drive pin
    digitalWrite(ledPin,HIGH);//indicate when we're pulsing
    digitalWrite(stepperPins[1],HIGH);		


    for(j=0;j<(pow(2,8)-1);j++ ){
        dumb+=1;
    }
    //Serial.println(micros()-t);
    
    digitalWrite(stepperPins[1],LOW);
    //keep track of which step we are on
    if(dir==1){
      n++;
  //    Serial.print("+");
    }
    else if(dir ==-1){
      n--;
 //     Serial.print("-");
    }
        digitalWrite(ledPin,LOW);

  }
  //////////////////////////////////////////////////////////////////////////////////
 
}
//-------------------------------------------------------------------------------------------


//>>>>>>>>>>>>>>> getters and updaters <<<<<<<<<<<<<<<<
void updateSetpoint(){
  //get the new setpoint and set the setpointChanged flag accordingly
  int val = 0;
  float newSetpoint;

  //read the setpoint from the designated pin
  val = analogRead(setPin);

  //compute new setpoint by linear interpolation between the limits
  newSetpoint = angVelLim[0]+(angVelLim[1]-angVelLim[0])*float(val)/1023.0;

  //angularVelocitySetpoint = 500;//comment for real deal

  /////////////ncomment for real deal/////////
  if (abs(angularVelocitySetpoint-newSetpoint)>setpointChangeThres){
    setpointChanged = true;
    Serial.print("new setpoint set at ");
    Serial.print(newSetpoint);
    Serial.println(" Hz");
    delay(10);
  }
  else
    setpointChanged = false;

  angularVelocitySetpoint = newSetpoint;
  //////////////////////////////////////////
}

void updateAngularVelocity(){
  //compute the angular velocity in rad/s from the array of periods T
  float val = 0;
  int nValid = 0;
  for(int i = 0;i<bufferLen;i++){
    if(T[i]>4200){
      val+= T[i]; 
      nValid++;
    } 
  }

  angularVelocity = 1e6*float(nValid)/float(val);
}

void updateError(){
  //push the error buffer
  for(int i = 0; i<errLen-1;i++){
    E[i] = E[i+1];
  }
  //compute current error sample
  E[errLen-1] = angularVelocity - angularVelocitySetpoint;

  //integral error
  Ei+=E[errLen-1]; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  needs delta t to be used with Ep !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  //proportional error
  Ep= E[errLen-1];
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! derivative computation goes here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  //				pass for now							//
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

}
//>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<

//::::::::::::::::::::            SETUP             :::::::::::::::::::::::::::::::
void setup() {

  //////////////// initialization ///////////////////////////
  pinMode(ledPin, OUTPUT);
  pinMode(stepperPins[0], OUTPUT);
  pinMode(stepperPins[1], OUTPUT);
  angVelSetpointRez = (angVelLim[1] - angVelLim[0])/256.0;
  n = maxStep;//init to maxStep so we can guarantee to step ccw until the home switch hits
  //n = maxStep/2;//initialize to maxStep for real deal


  Serial.begin(9600);
  Serial.print("homed = ");
  Serial.println(homed,BIN);
  Serial.println("initializing");
  delay(5);
  /////////////////////////////////////////////////////////

  //////////////// home the pot /////////////////////////

  if(limitState==LOW)
    digitalWrite(limitPin,HIGH);

  if (digitalRead(limitPin)==limitState){
    Serial.print("limit pin is ");
    Serial.print(limitState,BIN);
    Serial.println(", turning CW");
    while(digitalRead(limitPin)==limitState){
      takeStep(CW,10);
      delay(100);
    }
    takeStep(CW,50);
    n=0;
    homed = true;
    Serial.println("homed");  
  }
  else{
    Serial.print("limit pin is ");
    Serial.print(!limitState);
    Serial.println(", turning CCW");

    Serial.print("homed = ");
    Serial.println(homed,BIN);	
    while(true){
      takeStep(!CW, 5);
      delay(200);

      if(digitalRead(limitPin)==limitState){
          delay(1000);
        while(true){
          takeStep(CW, 5);
          delay(1000);
          if(digitalRead(limitPin)==!limitState)
            break;
          }
          takeStep(CW, 5);
          homed = true;
          n=0;
          break;
      }


    }
    if(homed)
      detachInterrupt(0);
    else
      Serial.println("THIS SHOULD NEVER EXECUTE: something went wrong and the device is not homed");

  }
  /////////////////////////////////////////////////////////

  //////////////// initialize the sensing buffer ////////////////
  then = micros();
  attachInterrupt(1, sensingInterruptHandler, RISING);

  while(idx<bufferLen){
    Serial.print("filling buffer idx = ");
    Serial.println(idx);
    digitalWrite(ledPin, HIGH);
    delay(100);
    digitalWrite(ledPin, LOW);
    delay(100);
  }
  ///////////////////////////////////////////////////////////////	
}
//::::::::::::::::::::            END SETUP             :::::::::::::::::::::::::::::::


void loop() {

  //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  // 								OPEN LOOP CONTROL
  //             if there is significant change in the velocity setpoint, 
  //				drive to the corresponding openLoop step setpoint
  if(!constantSetpoint)
    updateSetpoint();
/*
  if (setpointChanged){
    Serial.print("setpoint change detected: ");
    Serial.println(angularVelocitySetpoint);
    openLoopStepSetpoint = getOpenLoopStepSetpoint(angularVelocitySetpoint);
    if(n<openLoopStepSetpoint){
      DIR = CW;
    }
    else if(n>openLoopStepSetpoint){
      DIR = !CW;
    }
    takeStep(DIR, abs(openLoopStepSetpoint-n));
    setpointChanged = false; 
  }
  */
  //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\     

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // 								CLOSED LOOP CONTROL						//
  //             using the thresholds on proportional, integral, and 									//
  //		derivative action, decide whether to step the motor									//
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  updateAngularVelocity();
  updateError(); 

  //using only I control to start
  if(Ei>iThres){
    takeStep(!CW, 1);
  }
  else if(Ei<-iThres){
    takeStep(CW, 1);
  }

  delay(100);

  if(cnt%10){
   // Serial.print("wheel frequency: ");
    Serial.print(angularVelocity);
    Serial.print("\n");
 //   Serial.println(" Hz");
    delay(10);	
  }
  cnt++;
  //debug statements
  //Serial.println(T[idx%bufferLen]);
  //Serial.println(Ei);
  //Serial.println(angularVelocity/(2*3.14));
  //delay(1);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

}

int main(void)
{
	init();

	setup();
    
	for (;;)
		loop();
        
	return 0;
}

