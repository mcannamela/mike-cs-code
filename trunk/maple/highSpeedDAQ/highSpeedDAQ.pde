//interrupt handler
void handler(void);

//sampling frequency in Hz, and period in us
double samplingFrequency = 100000;
uint32 uPeriod = 0;

//number of samples, factor by which to oversample
const uint32 NSAMPLES =  5000;
const uint32 NOVERSAMPLE = 1;
float nOversample = NOVERSAMPLE;



uint32 usbDelay = 0;

//pin designations
const uint8 LED_PIN =  13;
const int analogInPin = 15; 


int toggle = 0; 

//counters for samples and oversamples
unsigned int oversampleCount = 0;
unsigned int sampleCount = 0;


//allocate an array to hold the data and a variable to hold the current sample
uint16 dataBuffer[NSAMPLES];
uint16 sample;

//flag tells us whether we are currently writing to the USB or not
bool usbFlag = false;

void setup() {
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
    uPeriod = uint32(1e6/(samplingFrequency*NOVERSAMPLE));
    if (uPeriod==0)
      uPeriod  =1;
    
    Timer2.setChannel1Mode(TIMER_OUTPUTCOMPARE);
    Timer2.setPeriod(uPeriod); // in microseconds
    Timer2.setCompare1(1);      // overflow might be small
    Timer2.attachCompare1Interrupt(handler);
}

void loop() {
  /*SerialUSB.print("usbDelay is: ");
  SerialUSB.print(usbDelay);
  SerialUSB.print("\n");*/
}

void handler(void) {

    if (!usbFlag)
    {
      //sum the oversamples to make one sample; this is safe because adc is 12 bit and we have a 16 bit int
      sample += analogRead(analogInPin);
      
      if( oversampleCount==(NOVERSAMPLE-1))
      {
       // SerialUSB.println(oversampleCount);
        oversampleCount = 0;
        dataBuffer[sampleCount] = sample;
        sample = 0;
        if (sampleCount == (NSAMPLES - 1))
        {
          usbFlag = true;
          SerialUSB.print(uPeriod);
          SerialUSB.print(" ");
          for(int i = 0; i<NSAMPLES; i++)
          {
            SerialUSB.print(dataBuffer[i]);
            SerialUSB.print(" ");
            dataBuffer[i] = 0;
          }
          SerialUSB.print("\n");
          sampleCount = 0;
          delay(usbDelay);
          usbFlag = false;
        }
        else
        {
          sampleCount+=1;
        }
      }
      else
      {
        //toggle ^= 1;
        //digitalWrite(LED_PIN, toggle);
        //SerialUSB.println(oversampleCount);
        oversampleCount+=1;
        
      }
    }

} 


