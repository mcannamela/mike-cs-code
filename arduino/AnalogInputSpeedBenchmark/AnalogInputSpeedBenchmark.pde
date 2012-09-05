int ledPin = 13;
int potPin = 2;    // select the input pin for the potentiometer
int controlPin = 11;

int val = 0;       // variable to store the value coming from the sensor
long start = 0;
long stop = 0;
long duration = 0;
void setup() {
  pinMode(potPin, INPUT);
  pinMode(ledPin, OUTPUT);  // declare the ledPin as an OUTPUT
  Serial.begin(9600); 
}

void loop() {
  start = millis();
  digitalWrite(ledPin, HIGH);
  for(int i = 0; i<10000; i++){
    val = analogRead(potPin);    // read the value from the sensor
    analogWrite(controlPin, 125);           // sets the value (range from 0 to 255) 
  }
  stop = millis();
  duration = stop-start;
  Serial.print("10000 samples read in ");
  Serial.print(duration);
  Serial.println("milliseconds");
  digitalWrite(ledPin, LOW);
  delay(300);
}
