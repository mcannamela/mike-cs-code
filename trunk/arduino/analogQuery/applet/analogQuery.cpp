#include "WProgram.h"
void setup();
void loop();
int inPin = 0;
void setup() {
  // initialize serial communication:
  Serial.begin(9600); 
  // initialize the LED pins:
}


void loop() {
  if (Serial.available() > 0) {
    int inByte = Serial.read();

    switch (inByte) {
    case 't':    

    case 'T':    
      for(int i = 0;i<2;i++){
        Serial.print(map(analogRead(i),0,1023,0,5040));
        Serial.print("  ");
      }
      Serial.print('\n');
      break;
    default:
      Serial.println('bad command: valid strings are "t" and "T"');
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

