
/////////////// pin declarations ////////////////////////
const int sensorInterruptPin = 3;
const int ledPin = 13;    // select the pin for the LED
const int stepperPins[2] = {
  4, 5}; 
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
//////////////////////////////////////////////////////////////////////////

const boolean CW = LOW; //write CW to stepperPins[0] to drive clockwise
volatile int n = 256;				//number of steps CW from the home position
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! calibrate!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
const int maxStep = 500;  		//max number of steps we can take from the home position
///////////angular velocity state, setpoint, and parameters//////////////
float angularVelocity = 0;				//current angular velocity computed from the measured periods

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

unsigned long cnt = 0;
//:::::::::::::::::::::: INTERRUPTS:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::	

void sensingInterruptHandler(){
  // when the reflectance goes high, record the time since the last spike
  noInterrupts();
  now = micros();
  //Serial.print("sensing interrupt, now - then = ");
  //Serial.println(now-then);

  if(now>then){
    T[idx%bufferLen] = now-then;
    t[idx%bufferLen] = now;
    idx++;
  }
  then = now;	
  interrupts();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::						



//>>>>>>>>>>>>>>> getters and updaters <<<<<<<<<<<<<<<<


void updateAngularVelocity(){
  //compute the angular velocity in rad/s from the array of periods T
  float val = 0;
  for(int i = 0;i<bufferLen;i++){
    val+= 1e6*2*3.14/(float)T[i];  
  }
  val/=bufferLen;
  angularVelocity = val;
}

//>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<
void takeStep(boolean dir, int nStep){
  //drive the pot nSteps in the dir direction

  //set direction ///
  unsigned long dummyIdx, dumb=0;

  //LOW for CW, HIGH for CCW
  digitalWrite(stepperPins[0],dir);

  //////////////////

  //\\\\\\\\\\\\\\\\\take the specified number of steps, incrementing n as necessary\\\\\\\\\\\\\\
  digitalWrite(stepperPins[1],LOW);
  for(int i = 0; i<nStep;i++){
    //check for end conditions
    if(n==0 && dir==!CW){
      Serial.println("ccw limit!");
      break;
    }
    else if(n==maxStep-1 && dir==CW){
      Serial.println("cw limit!");
      break;
    }

    //////////pulse the drive pin//////////////
    digitalWrite(stepperPins[1],HIGH);		

    //hardcode a delay by adding since delays won't work in interrupts!
    for(dummyIdx=0;dummyIdx<stepDelayLoopCount;dummyIdx++ )
      dumb+=1;

    digitalWrite(stepperPins[1],LOW);
    /////////////////////////////////////////


    //keep track of which step we are on
    if(dir==CW){
      n++;
      //    Serial.print("+");//debug statement
    }
    else if(dir ==!CW){
      n--;
      //     Serial.print("-");//debug statement
    }

  }
//::::::::::::::::::::            SETUP             :::::::::::::::::::::::::::::::
void setup() {

  //////////////// initialization ///////////////////////////
  pinMode(ledPin, OUTPUT);
  
  Serial.begin(9600);
  Serial.println("initializing");
  delay(5);
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
    
   updateAngularVelocity();
   updateError(); 
   
   //using only I control to start
   if(Ei>iThres){
   takeStep(!CW, 1);
   }
   else if(Ei<-iThres){
   takeStep(CW, 1);
   }
   
   
   
   if(cnt%100){
	Serial.print("wheel frequency: ");
	Serial.print(angularVelocity/(2*3.14));
	Serial.println(" Hz");
	}
	delay(10);
	else
		delay(10)
  
  cnt++;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

}

