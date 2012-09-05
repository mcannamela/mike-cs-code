
int nTCs = 2;//number of thermocouples
int SOpin[] = {2,4};//slave out
int CSpin[] = {3,5};//control 
int SCKpin = 6;//clock
int binaryArray[16]; 

void setup() {
  // initialize serial communication:
  Serial.begin(9600);
  for(int i = 0;i < nTCs;i++){
    pinMode(SOpin[i], INPUT);
    pinMode(CSpin[i], OUTPUT);
  }
  pinMode(SCKpin, OUTPUT);
}

float cookBinaryArray(void){
  int val = 0;
  int thePower = 0;
  float T = 0;
  for(int i = 12; i > 0; i--){
    val += binaryArray[i]*pow(2,thePower);
    thePower++;
  }
  T = float(val)/4;
  return T;
} 

float queryTC(int SO, int CS){
  digitalWrite(CS, LOW);
  //delay(1);
  for(int i = 0; i <16; i++){ 
    digitalWrite(SCKpin, HIGH);
    //delay(1);
    binaryArray[i] = digitalRead(SO);
    digitalWrite(SCKpin, LOW);
    //delay(1);


    //Serial.print(binaryArray[i]);
    //Serial.print(" ");

  }
  digitalWrite(CS, HIGH);


  //Serial.print('\n');

  return cookBinaryArray();
}


void loop() {
  if (Serial.available() > 0) {
    int inByte = Serial.read();

    switch (inByte) {
    case 't':    

    case 'T':    
      for(int i = 0;i < nTCs;i++){
        Serial.print(queryTC(SOpin[i], CSpin[i]));
        Serial.print("  ");
        delay(10);
      }
      Serial.print('\n');
      break;
    default:
      Serial.println('bad command: valid strings are "t" and "T"');
    }
  }
  delay(10);
}


