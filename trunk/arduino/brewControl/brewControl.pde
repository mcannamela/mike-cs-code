
//************************************************************
//***********************matt's code**************************
//************************************************************
/////////////////////////////////////////////////////////////////
//                        Pin Declarations                     //
/////////////////////////////////////////////////////////////////
int id0 = 7;
int id1 = 8;
int id2 = 9;
int id3 = 10;
int id7 = 11;
int wrt = 12;
int mod = 13;

int mashDisplay = 4;
int hltDisplay = 0;
int setpointSwitch = A3;
int showSetpoint = LOW;
/////////////////////////////////////////////////////////////////
//                        Global Variables                     //
/////////////////////////////////////////////////////////////////
//Digit holder array for each digit of the number
int _n[4] = {8,8,8,8};
//Decimal holder array for the state of each decimal point
int DPs[4] = {LOW, LOW, LOW, LOW};
int mattsCnt = 0;

/*
float debug_T[] = {00.00,
                0.001,
                1.111,
                2.222,
                33.33,
                44.44,
                55.55,
                666.6,
                777.7,
                888.8,
                9999.,
                99999};
*/

/*
/////////////////////////////////////////////////////////////////
//                         Main Loop                           //
/////////////////////////////////////////////////////////////////
void loop()
{
  dispNumber(debug_T[mattsCnt],0);
  dispNumber(debug_T[mattsCnt],4);
  delay(1000);
  mattsCnt = (mattsCnt+1)%12;
}
*/

/////////////////////////////////////////////////////////////////
//                      dispNumber Function                    //
/////////////////////////////////////////////////////////////////
//this function will display the number passed on the 4 digit  //
//7 segment display. the number must be greater then 0.001 and //
//less than 10000 or else 8.8.8.8. will be displayed. The int  //
//passed for disp should be either a 0 or a 4 to address the   //
//first or second 4 digit seven segment.                       //
/////////////////////////////////////////////////////////////////
void dispNumber(double num, int disp)
{
  //sets the decimal point and separates digits
  if(num<10000)             //
  {                         //
  if(num<1000)              //
    {                       //these checks determine where the
      if(num<100)           //decimal point should be placed by
      {                     //determining how big the number is
        if(num<10)          //
        {                   //
          if(num<0.001)     //
          {
            //error number too small display will be 0.0.0.7.
            _n[0] = 0;
            _n[1] = 7;
            _n[2] = 0;
            _n[3] = 0;
            DPs[0] = LOW;
            DPs[1] = LOW;
            DPs[2] = LOW;
            DPs[3] = LOW;
          }
          else
          {
            //0.001<=num<10 decimal set as X.XXX
            sepDigs(num);
            //num is sent as is because sepDigs accepts numbers
            //in the range of 0.001 to 9.999 and the decimal point
            //is set by the code below
            DPs[0] = LOW;
            DPs[1] = HIGH;
            DPs[2] = HIGH;
            DPs[3] = HIGH;
          }
        }
        else
        {
          //10<=num<100 decimal set as XX.XX
          sepDigs(num/10.0); 
          //num is divided by 10 because sepDigs accepts numbers
          //in the range of 0.001 to 9.999 and the decimal point
          //is set by the code below
          DPs[0] = HIGH;
          DPs[1] = LOW;
          DPs[2] = HIGH;
          DPs[3] = HIGH;
        }
      }
      else
      {
        //100<=num<1000  decimal set as XXX.X
        sepDigs(num/100.0);
        //num is divided by 100 because sepDigs accepts numbers
        //in the range of 0.001 to 9.999 and the decimal point
        //is set by the code below
        DPs[0] = HIGH;
        DPs[1] = HIGH;
        DPs[2] = LOW;
        DPs[3] = HIGH;
      }
    }
    else
    {
      //1000<=num<10000  decimal set as XXXX.
      sepDigs(num/1000.0);
      //num is divided by 1000 because sepDigs accepts numbers
      //in the range of 0.001 to 9.999 and the decimal point
      //is set by the code below
      DPs[0] = HIGH;
      DPs[1] = HIGH;
      DPs[2] = HIGH;
      DPs[3] = LOW;
    }
  }
  else
  {
    //error number too big display will be 4.9.1.4.
    _n[0] = 4;
    _n[1] = 9;
    _n[2] = 1;
    _n[3] = 4;
    DPs[0] = LOW;
    DPs[1] = LOW;
    DPs[2] = LOW;
    DPs[3] = LOW;
  }
  //loop over digits
  for(int i=(0+disp);i<(4+disp);i++)  //this calls lightNumber 
  {                                   //function for each digit
    lightNumber(i,_n[i-disp],DPs[i-disp]);
  }
}

/////////////////////////////////////////////////////////////////
//                        sepDigs Function                     //
/////////////////////////////////////////////////////////////////
//this function accepts a number between 0.001 and 9.999 and   //
//places each digit of the number into a separate index of the //
//_n array. if the number has more than 4 significant digits    //
//additional digits will be truncated.                         //
/////////////////////////////////////////////////////////////////
void sepDigs(double adjnum)
{
  adjnum = adjnum + 0.00001;//this is needed because for some
                            //reason 1.000 will become 0 when
                            //floored, this is a hack
  //this function places each digit of adjnum into _n
  _n[0] = floor(adjnum);
  _n[1] = floor((adjnum-_n[0])/0.1);
  _n[2] = floor((adjnum-_n[0]-(_n[1]*.1))/0.01);
  _n[3] = round((adjnum-_n[0]-(_n[1]*.1)-(_n[2]*.01))/0.001);
}

/////////////////////////////////////////////////////////////////
//                     lightNumber Function                    //
/////////////////////////////////////////////////////////////////
//this function accepts a digit to light, a number to be       //
//displayed, and the state of the decimal point. It sends      //
//commands to the ICM7218 to change the digit and change the   //
//number the digit is displaying.                              //
/////////////////////////////////////////////////////////////////
void lightNumber(int dig, int num, int dec)
{
  changeDig(dig);       //changes active digit to the one sent
  changeNum(num,dec);  //flashes the decimal point
}
    
/////////////////////////////////////////////////////////////////
//                      changeDig Function                     //
/////////////////////////////////////////////////////////////////
//this function changes the active digit to the one that is    //
//sent. it sends the command to the ICM7218 to change the      //
//active digit.
/////////////////////////////////////////////////////////////////
void changeDig(int dig)
{
    digitalWrite(id0, bitRead(dig,0)); //these set the digit
    digitalWrite(id1, bitRead(dig,1)); //coded in binary
    digitalWrite(id2, bitRead(dig,2)); //
    digitalWrite(id3, LOW);                 //Selects the RAM
                                            //bank as BankB
    digitalWrite(mod, HIGH);                //Mode high for
                                            //changing digits
    digitalWrite(id7, LOW);                 //Tells the ICM7218
                                            //that this is a 
                                            //digit change and
                                            //not a new number
    digitalWrite(wrt, HIGH);                //tells the ICM7218
                                            //to accept the new
                                            //active digit
    digitalWrite(wrt, LOW);
}
/////////////////////////////////////////////////////////////////
//                      changeNum Function                     //
/////////////////////////////////////////////////////////////////
//this function changes the number on the currently active     //
//digit. It sends the hex coded command to the ICM7218 to      //
//change the number.                                           //
/////////////////////////////////////////////////////////////////
void changeNum(int num, int decimal)
{
    digitalWrite(mod, LOW);            //mode low for changing
                                       //numbers
    digitalWrite(id0, bitRead(num,0)); //
    digitalWrite(id1, bitRead(num,1)); //These set the number
    digitalWrite(id2, bitRead(num,2)); //coded in binary
    digitalWrite(id3, bitRead(num,3)); //
    digitalWrite(id7, decimal);        //this sets the decimal
    digitalWrite(wrt, HIGH);           //tells the ICM7218 to
                                       //accept the new number
    digitalWrite(wrt, LOW);
}
//************************************************************
//************************************************************


////////////////////////////////////////////
///////////  control parameters ////////////
long controlDelay = 60000; //time between control action in ms
long oldTime = millis();
long elapsed = millis()-oldTime;
////////////////////////////////////////////

//knob variables
int knobPin = 0;
int knobVal = 0;
float knob_mash_temp[16] = {
  130, 131, 132, 133,
  150, 151, 152, 153,
  154, 155, 156, 157,
  158, 165, 170, 80
};
/*
float knob_mash_temp[16] = {
  75, 75, 75, 75,
  75, 75, 75, 75,
 75, 75, 75, 75,
  75, 75, 75, 75,
};*/
float knob_HLT_temp[16] = {
  170, 170, 170, 170,
  170, 170, 170, 170,
  170, 170, 170, 170,
  170, 165, 170, 80
};


//TC reading variables
int nTCs = 2;//number of thermocouples
int SOpin[] = {2,4};//slave out
int CSpin[] = {3,5};//control 
int SCKpin = 6;//clock
int binaryArray[16]; //buffer for the serial communication with tc amps

//heater pins
int HLTPin = A1;
int mashPin = A2;

float T_set[2] = {27.0, 27.0};

int hltIdx = 0;
int mashIdx = 1;

float T[2] = {0.0,0.0};



void setup() {
  // initialize serial communication:
  Serial.begin(9600);

  //initialize heater pins
  pinMode(HLTPin, OUTPUT);
  pinMode(mashPin, OUTPUT);
  //pinMode(7, OUTPUT);
  //pinMode(8, OUTPUT);

  //start out heaters off!
  digitalWrite(HLTPin, LOW);
  digitalWrite(mashPin, LOW);

  //initialize TC pins
  for(int i = 0;i < nTCs;i++){
    pinMode(SOpin[i], INPUT);
    pinMode(CSpin[i], OUTPUT);
  }
  pinMode(SCKpin, OUTPUT);
  
  //************* matt's setup ***************
    //sets all pins as outputs and in a low state
  pinMode(id0, OUTPUT);
  pinMode(id1, OUTPUT);
  pinMode(id2, OUTPUT);
  pinMode(id3, OUTPUT);
  pinMode(id7, OUTPUT);
  pinMode(wrt, OUTPUT);
  pinMode(mod, OUTPUT);
  digitalWrite(id0, LOW);
  digitalWrite(id1, LOW);
  digitalWrite(id2, LOW);
  digitalWrite(id3, LOW);
  digitalWrite(id7, LOW);
  digitalWrite(wrt, LOW);
  digitalWrite(mod, LOW);
  
  pinMode(setpointSwitch, INPUT);
  //*******************************************
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::
//::::::::::: setpoint updating :::::::::::::::::::::
//:::::::::::::::::::::::::::::::::::::::::::::::::::

//get the value of the knob
int queryKnob(int thePin){
  int pinVal, small;
  float num;

  float Ra = 4000.0;
  float Rb = 330.0;

  small = 10;

  pinVal = analogRead(thePin);

  if (pinVal<=small)
    num = 0;
  else
    num = (Ra/Rb)*float(pinVal)/(1023.0-float(pinVal)); 
    
//  Serial.println(pinVal);
  return round(num);
}


//unit conversion for temperature
float F2C(float F){
  return (F-32)*5.0/9.0;
}
float C2F(float C){
  return 9*C/5.0+32;
}

void updateSetpoints(void){
    knobVal = queryKnob(knobPin);
    T_set[mashIdx] = F2C(knob_mash_temp[knobVal]);
    T_set[hltIdx] = F2C(knob_HLT_temp[knobVal]);

}
//:::::::::::::::::::::::::::::::::::::::::::::::::::






/////////////////////////////////////////////////////
/////////////  TC Query Functions///////////////////
////////////////////////////////////////////////////
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
  float T;
  digitalWrite(CS, LOW);
  //delay(1);
  for(int i = 0; i <16; i++){ 
    digitalWrite(SCKpin, HIGH);
    //delay(1);
    binaryArray[i] = digitalRead(SO);
    digitalWrite(SCKpin, LOW);
    //delay(1);


   // Serial.print(binaryArray[i]);
    //Serial.print(" ");

  }
  digitalWrite(CS, HIGH);


  //Serial.print('\n');
  T = cookBinaryArray();
  return T;
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////




void loop() {

  //update temperature

  for(int i = 0;i < nTCs;i++){
    T[i] = .9*T[i]+.1*queryTC(SOpin[i], CSpin[i]);
    
    //Serial.print(T[i]);
    //Serial.print("    ");

    delay(75);
  } 
//  Serial.print('\n');
  updateSetpoints();



  ///////////////////////////////////////////////// 
  ///////////////////control logic/////////////////
  /////////////////////////////////////////////////
  if (digitalRead(setpointSwitch)==showSetpoint)
  {
    dispNumber(C2F(T_set[mashIdx]), mashDisplay);
    delay(10);
    dispNumber(C2F(T_set[hltIdx]), hltDisplay);
    delay(10);
    
  }
  else
  {
    dispNumber(C2F(T[mashIdx]), mashDisplay);
    delay(10);
    dispNumber(C2F(T[hltIdx]), hltDisplay);
    delay(10);
  }
  
  elapsed = millis()-oldTime;
  
 /* Serial.println(elapsed);
  delay(10);	
  Serial.println(controlDelay);
  delay(10);
  Serial.println(elapsed>=controlDelay);
  delay(10);*/
  
  if (elapsed>=controlDelay){
    oldTime = millis();   
    


    /* Serial.print(T[0]);
     Serial.print("  ");
     Serial.print(T_set[0]);
     Serial.print("  ");
     Serial.print(T[1]);
     Serial.print("  ");
     Serial.println(T_set[1]);*/



    //Serial.println("controlling ^^"); 
	//delay(10);	
	
    if (T[hltIdx]<T_set[hltIdx]){
      //Serial.println("writing HLT high");
      digitalWrite(HLTPin, HIGH);
      //digitalWrite(7, HIGH);

    }

    else{
      //Serial.println("writing HLT low");
      digitalWrite(HLTPin, LOW);
      //digitalWrite(7, LOW);
    }

    if (T[mashIdx]<T_set[mashIdx]){
      //Serial.println("writing Mash high");
      digitalWrite(mashPin, HIGH);
      //digitalWrite(8, HIGH);
    }
    else{
      //Serial.println("writing Mash low");
      digitalWrite(mashPin, LOW);
      //digitalWrite(8, LOW);
    }
  }
  ////////////////////////////////////////////////


  //can query the temperatures
  if (Serial.available() > 0) {
    int inByte = Serial.read();

    switch (inByte) {

    
      //report knob position  
    case 'k':
    case 'K':
      Serial.println("knob value");
      Serial.print(knobVal);
      Serial.print("    ");
      Serial.print(analogRead(knobPin));
      Serial.print('\n');
      break;

      //report TC temperature
    case 't':    
    case 'T':    
      for(int i = 0;i < nTCs;i++){
        Serial.print(T[i]);
        Serial.print("  ");
        //Serial.print(T_set[i]);
        //Serial.print("  ");
        delay(10);
      }
      Serial.print('\n');
      break;
      
    case 's':    
    case 'S':    
        Serial.println("mashSet_C mashSet_F HLTSet_C HLTSet_F");
        Serial.print(T_set[mashIdx]);
        Serial.print("  ");
        Serial.print(knob_mash_temp[knobVal]);
        Serial.print("  ");
        Serial.print(T_set[hltIdx]);
        Serial.print("  ");
        Serial.print(knob_HLT_temp[knobVal]);
        Serial.print("  ");
        delay(10);
      
      Serial.print('\n');
      break;
      
    //report pinout
    case 'p':
    case 'P':
      Serial.println("HLTPin    mashPin    SOPin0    SOPin1    hltIdx    mashIdx");
      Serial.print(HLTPin);
      Serial.print("         ");
      Serial.print(mashPin);
      Serial.print("         ");
      Serial.print(SOpin[0]);
      Serial.print("         ");
      Serial.print(SOpin[1]);
      Serial.print("         ");
      Serial.print(hltIdx);
      Serial.print("         ");
      Serial.print(mashIdx);

      Serial.print('\n');
      break;

    default:
      Serial.println("bad command: valid strings are 's', 'S', 'p','P', 'k', 'K', 't' and 'T'");
    }
  }
  delay(10);
}




