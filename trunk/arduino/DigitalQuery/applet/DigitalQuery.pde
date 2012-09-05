int inPin = 3;
void setup() {
  // initialize serial communication:
  Serial.begin(9600); 
  // initialize the LED pins:
}


void loop() {
      digitalWrite(13,digitalRead(inPin));
      delay(100);
}
