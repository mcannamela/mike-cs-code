//:::::::::::::::::::::::::::::initialize index:::::::::::::::::::::::::::::::::::::::::::::::
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
      Serial.print(analogRead(sensorPin));
      Serial.print("  ");
      Serial.print("\n");
    }
  }
}
