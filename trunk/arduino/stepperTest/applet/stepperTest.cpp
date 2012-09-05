#include "WProgram.h"
void takeStep(int dir, int nStep);
void setup();
void loop();
int ledPin = 13;    // select the pin for the LED
int stepperPins[2] = {4, 5 }; //stepper for pot controlled on these 
int n = 128;
int maxStep = 256;
long stepDelay = 10000;
char d;
char c;
int N;
int dumb=0;
unsigned long t;
unsigned long j,i;


void takeStep(int dir, int nStep){
  //drive the pot nSteps in the dir direction
  // dir = -1 : step CCW
  // dir =  1 : step CW
  //Serial.println("taking step"); 
  //set direction ///
  if(dir==1)
    digitalWrite(stepperPins[0],HIGH);
  else if(dir ==-1)
    digitalWrite(stepperPins[0],LOW);
  //////////////////

  ///take the specified number of steps, incrementing n as necessary///
  digitalWrite(stepperPins[1],LOW);
  for(int i = 0; i<nStep;i++){
    //check for end conditions
    /*if(n==0 && dir==-1){
     // Serial.println("ccw limit!");
      break;
    }
    else if(n==maxStep-1 && dir==1){
     // Serial.println("cw limit!");
      break;
    }*/

    //pulse the drive pin
    digitalWrite(ledPin,HIGH);//indicate when we're pulsing
    digitalWrite(stepperPins[1],HIGH);		
    /*
    //wait for the step to occur before sending another pulse
    if (stepDelay>10000)
      delay(stepDelay/1000); //change to delay if stepDelay is large
    else
      delayMicroseconds(stepDelay); //change to delay if stepDelay is large
    */
   // t = micros();
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



void setup() {

  //////////////// initialization ///////////////////////////
  pinMode(ledPin, OUTPUT);
  pinMode(stepperPins[0], OUTPUT);
  pinMode(stepperPins[1], OUTPUT);
  Serial.begin(9600);
  ///////////////////////////////////////////////////////////////	
}
//::::::::::::::::::::            END SETUP             :::::::::::::::::::::::::::::::


void loop() {


  if(Serial.available()>0){
    
   d = Serial.read();
   i= 0;
   N = 0;
   while(Serial.available()>0){
     c = Serial.read();
     N+=  int(c);
   }
     
   Serial.println(N);
   if( d=='+')
     takeStep(1,N);
   else if (d=='-')
     takeStep(-1,N);
   
  }
  delay(50);

  ///////////////////////////////////////////////////////////

}

int main(void)
{
	init();

	setup();
    
	for (;;)
		loop();
        
	return 0;
}

