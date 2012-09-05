//:::::::::::::::::::::::::::::initialize index:::::::::::::::::::::::::::::::::::::::::::::::
#include "WProgram.h"
void setup();
void loop();
int sensorPin = 2;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void setup() 
{
  Serial.begin(9600); // start serial communication
}

void loop() 
{
  if (Serial.available() > 0) 
  {
    int inByte = Serial.read();
    switch (inByte) {
    case 't':    
      Serial.print(map(analogRead(sensorPin),0,1023,0,5000));
      Serial.print("  ");
      Serial.print("\n");
    }
  }
}

int main(void)
{
	init();

	setup();
    
	for (;;)
		loop();
        
	return 0;
}

