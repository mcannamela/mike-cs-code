int pulsePin = 4;
int pwmPin = 5;
int triggerPin = 5;
int onTime = 10;
int offTime = 46;



void setup() {
  // initialize pulse pin
  pinMode(pulsePin, OUTPUT);
  pinMode(pwmPin, OUTPUT);
  analogWrite(pwmPin, 128);
}


void loop() {
  digitalWrite(pulsePin,LOW);
  delay(offTime);
//delay(200);
  digitalWrite(pulsePin, HIGH);
  delay(onTime);
  //delay(200);
}
