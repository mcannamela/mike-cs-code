void sensingInterruptHandler(){
  // when the reflectance goes high, record the time since the last spike
  noInterrupts();
  Serial.println(micros());
  interrupts();
}
void setup(){
Serial.begin(9600);
attachInterrupt(1, sensingInterruptHandler, RISING);
}
void loop(){
  delay(10);
}
