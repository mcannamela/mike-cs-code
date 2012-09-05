int pulsePin = 4;
void setup() {
  // initialize pulse pin
  pinMode(pulsePin, OUTPUT);
  digitalWrite(pulsePin, LOW);
}


void loop() {
  pinMode(pulsePin, OUTPUT);
  delayMicroseconds(100);
  pinMode(pulsePin, INPUT);
  delayMicroseconds(10);
}
