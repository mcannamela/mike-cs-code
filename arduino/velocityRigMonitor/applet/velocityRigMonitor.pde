char d;
char c;
int N;
int i,j,dumb;
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
const int maxStep = 1600;  		//max number of steps we can take from the home position
unsigned long stepDelayLoopCount = 255;
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
  int nValid = 0;
  for(int i = 0;i<bufferLen;i++){
    if(T[i]>4200){
      val+= T[i]; 
      nValid++;
    } 
  }

  angularVelocity = 1e6*2*3.14*float(nValid)/float(val);
}

//>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<
void takeStep(int dir, int nStep){
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
//::::::::::::::::::::            SETUP             :::::::::::::::::::::::::::::::
void setup() {

  //////////////// initialization ///////////////////////////

  pinMode(ledPin, OUTPUT);
  pinMode(stepperPins[0], OUTPUT);
  pinMode(stepperPins[1], OUTPUT);
  Serial.begin(9600);
  Serial.println("initializing");
  delay(5);
  /////////////////////////////////////////////////////////

  //////////////// initialize the sensing buffer ////////////////
  then = micros();
  attachInterrupt(1, sensingInterruptHandler, RISING);

  /*
  while(idx<bufferLen){
   Serial.print("filling buffer idx = ");
   Serial.println(idx);
   digitalWrite(ledPin, HIGH);
   delay(100);
   digitalWrite(ledPin, LOW);
   delay(100);
   }*/
  ///////////////////////////////////////////////////////////////	
}
//::::::::::::::::::::            END SETUP             :::::::::::::::::::::::::::::::


void loop() {

  updateAngularVelocity();



  if(cnt%100){
    // Serial.print("wheel frequency: ");
    Serial.print(angularVelocity/(2*3.14));
    Serial.print("\n");
    //   Serial.println(" Hz");

    delay(10);	
  }

  else
    delay(10);

  if(Serial.available()>0){

    d = Serial.read();

    N = 0;
    while(Serial.available()>0){
      c = Serial.read();
      N+=  int(c);
    }
    Serial.print(d);
    Serial.println(N);
    if( d=='+')
      takeStep(CW,N);
    else if (d=='-')
      takeStep(!CW,N);

  }
  delay(100);
  cnt++;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

}

