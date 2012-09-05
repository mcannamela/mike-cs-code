/////////////// pin declarations ////////////////////////
#include "WProgram.h"
void sensingInterruptHandler();
void setup();
void loop();
int inPin = 0;
const int sensorInterruptPin = 2;
const int ledPin = 13;
///////////////////////////////////////////////////////////

/////////////// declare data array, its length, and indexing variable//////////////
volatile unsigned long i;//counting variable
double fr; //Water flowrate in gpm
volatile unsigned int led; //led toggle variable
volatile unsigned long time ; 	// sensing input buffer of periods in microseconds
volatile unsigned long now, then; //hold clock values between sensing interrupts
////////////////////////////////////////////////////////////////////////////////////

//:::::::::::::::::::::: INTERRUPTS:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::	

void sensingInterruptHandler(){
  // when the pulse from water, record the time since the last spike
  noInterrupts();
  now = micros();
  if(now>then+3000000){
    time = now-then;
    then = now;	
  if(led==1)
    {
      digitalWrite(ledPin, HIGH);   // sets the LED on
      led = 0;
    }
    else
    {
       digitalWrite(ledPin, LOW);    // sets the LED off
       led = 1;
    }
  }
  interrupts();
}

void setup() {
  // initialize serial communication:
  Serial.begin(9600); 
  // initialize the LED pins:
  then = micros();
  pinMode(sensorInterruptPin,INPUT);
  attachInterrupt(0, sensingInterruptHandler, FALLING);
  pinMode(ledPin, OUTPUT);
  analogReference(DEFAULT);
}


void loop() {
  if (Serial.available() > 0) {
    int inByte = Serial.read();

    switch (inByte) {
    case 'w':    
    case 'W':    
      for(int i = 3;i<6;i++){
        Serial.print(map(analogRead(i),0,1024,0,1024));
        Serial.print("  ");
      }
      Serial.print('\n');
      break;
    
    case 'f':
    case 'F':
      fr = (1.0/time)*60000000;
      Serial.print(fr);
      Serial.print('\n');
      break;
      
    case 't':    
    case 'T':    
        Serial.print(map(analogRead(0),0,1024,0,5000));
        Serial.print("  ");
      Serial.print('\n');
      break;
      
    default:
      Serial.println('bad command: valid strings are "t" and "T"');
      break;

    }
  }
  delay(10);
}


int main(void)
{
	init();

	setup();
    
	for (;;)
		loop();
        
	return 0;
}

