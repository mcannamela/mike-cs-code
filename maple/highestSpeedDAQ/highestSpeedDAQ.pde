//interrupt handlers
#include "adc.h"
void sampleHandler(void);
void commHandler(void);

//sampling frequency in Hz, and period in us
//double samplingFrequency = 350000;
uint32 uPeriod = 5;
uint32 commPeriod = 0;

//number of samples, factor by which to oversample
const uint32 NSAMPLES =  5000;

//pin designations
const uint8 LED_PIN =  13;
const int analogInPin = 15; 


int toggle = 0; 

//counters for samples and oversamples
unsigned int cnt = 0;

//allocate an array to hold the data and a variable to hold the current sample
uint16 dataBuffer[NSAMPLES];

unsigned int usbDelay = 0;

void setup() {
   SerialUSB.begin();
    //estimate delay need for communication of data based on baudrate
    usbDelay = uint32((1.0/(115200.0/2.0/16.0/NSAMPLES))*1000);
    if (usbDelay<=1)
    {
      usbDelay = 2;
    }
    
    // Configure pins
    pinMode(analogInPin, INPUT_ANALOG);
    pinMode(LED_PIN, OUTPUT);
    
    //configure the timeer
   // uPeriod = uint32(1e6/samplingFrequency);
   // if (uPeriod==0)
   //   uPeriod  =1;
    
    commPeriod = uint32(uPeriod*NSAMPLES*.9);
    
    adc_set_sample_rate(ADC1, ADC_SMPR_1_5 );
    adc_set_sample_rate(ADC2, ADC_SMPR_1_5 );

    
    Timer2.setChannel1Mode(TIMER_OUTPUTCOMPARE);
    Timer2.setPeriod(uPeriod); // in microseconds
    Timer2.setCompare1(1);      // overflow might be small
    Timer2.attachCompare1Interrupt(sampleHandler);
    
    Timer3.setChannel1Mode(TIMER_OUTPUTCOMPARE);
    Timer3.setPeriod(commPeriod); // in microseconds
    Timer3.setCompare1(1);      // overflow might be small
    Timer3.attachCompare1Interrupt(commHandler);
}

void loop() {
///  SerialUSB.println(ADC_SMPR_1_5);
  /*SerialUSB.print("usbDelay is: ");
  SerialUSB.print(usbDelay);
  SerialUSB.print("\n");*/
  //SerialUSB.println(uPeriod);
  //SerialUSB.println(commPeriod);
}

void sampleHandler(void) {
      dataBuffer[cnt] = analogRead(analogInPin);
      cnt++; 
    }
    
void commHandler(void) {
  Timer2.pause();
  Timer3.pause();
  SerialUSB.print(uPeriod);
  SerialUSB.print(" ");
  for(int i = 0; i<NSAMPLES; i++)
  {
    SerialUSB.print(dataBuffer[i]);
    SerialUSB.print(" ");
    delay(1);
    dataBuffer[i] = 0;
  }
  SerialUSB.print("\n");
  cnt = 0;
  delay(usbDelay);
  
  Timer2.resume();
  Timer3.resume();
          
    }

