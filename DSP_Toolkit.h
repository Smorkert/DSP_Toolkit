//
//  DSP_Toolkit.h
//  DSP
//
//  Created by Steven Morkert on 7/16/16.
//  Copyright (c) 2016 Steven Morkert. All rights reserved.
//

#ifndef DSP_DSP_Toolkit_h
#define DSP_DSP_Toolkit_h

//includes
#include <math.h>

/*List of Functions:
 --------------------------------------------------------------------
 1) Moving Average
 2) FIR Filter: float and 16-bit fixed point
 3) IIR Filter: 16-bit fixed point
 4) Float-to-Fixed conversion
 5) LFO generator 
 6) Chorus Effect
 7) 16-bit Compression: Forward and Reverse
 ---------------------------------------------------------------------
 */


/* Moving Average, mAve(order):
 This function takes the moving average of a continuous input
 where N is the "order" of the filter/average
 Type: int16_t
 --------------------------------------------------------------
 unBuffer = [indexer, state,1,...,N];
 1) decalre object: mAve ave(N) where N is order of average
 2) call function: ave.average(input);
 */

struct mAve{
    mAve(uint16_t n){
        unBuffer = new int16_t[n+2]; //buffer
        N = n;
    }
    int16_t *unBuffer;
    uint16_t N;
    int32_t sum;
    //Moving Average Function
    int16_t average(int16_t input){
        //initialize register
        if(unBuffer[1] == 0){
            unBuffer[2+unBuffer[0]] = input;
            unBuffer[0] += 1;
            //take average
            if(unBuffer[0] == N){
                for(uint16_t i=0;i<N; i++){
                    sum += unBuffer[2+i];
                }
                //change state
                unBuffer[1] = 1;
                //reset indexer
                unBuffer[0] = 2;
                return sum/N;
            }
        }
        else{
            //compute current moving average
            sum = ((sum+input)-unBuffer[unBuffer[0]]);
            //shift samples
            unBuffer[unBuffer[0]] = input;
            unBuffer[0] += 1;
            //reset indexer
            if(unBuffer[0] >= N+2){
                unBuffer[0] = 2;
            }
        }
        return sum/N;
    }
};


/* FIR filter, fFIR(input, filter coefficients, order):
 This function implements an FIR filter with user defined filter
 coefficients and order uN
 Type: float
 ---------------------------------------------------------------
 buff = [indexer, output,1,...,N];
 1) declare filter object: fFIR_Filter LPF(N, fFIR_K) where N is
 the number of filter coefficients and fFIR_K are the filter
 coefficients
 2) Call function: output = LPF.fFIR(input);
 */

struct fFIR_Filter{
    //constructor (N, filter coefficients)
    fFIR_Filter(uint8_t n, float *X)
    {
        buff = new float[n+2]; //buffer
        fFIR_k = new float[n]; //filter coefficient array
        //copy in coefficients
        for(int i = 0; i<n; i++){
            fFIR_k[i] = X[i];
        }
        N = n; //store N
    }
    //declare constructor variables
    float *buff;
    float *fFIR_k;
    uint8_t N;
    //FIR function
    float fFIR(float input){
        //clear previous output
        buff[1] = 0;
        //store indexer
        uint16_t inc = buff[0];
        //store current sample in end of buffer
        buff[1+N-inc] = input;
        //convolve
        for(uint16_t i = 0; i<N; i++){
            buff[1]+= buff[1+N-inc]*fFIR_k[i];
            if(inc <= 0){inc = N;}
            inc -= 1;
        }
        buff[0] += 1;
        if(buff[0]>=N){buff[0]=0;}
        return buff[1];
    }
};

/* FIR filter, nFIR(input, filter coefficients, order);
 This function implements an FIR filter with user defined filter
 coefficients and order uN,
 Type: int16_t
 ---------------------------------------------------------------
 buff = [current sample, indexer, output,1,...,N];
 1) declare filter object: nFIR_Filter LPF(N, fFIR_K) where N is
    the number of filter coefficients and fFIR_K are the filter
    coefficients
 2) Call function: output = LPF.nFIR(input);
 */

struct nFIR_Filter
{
    //constructor (N, filter coefficients)
    nFIR_Filter(uint8_t n, float *X)
    {
        x = new float[n]; //floating point coefficient array
        buff = new int16_t[n+2]; //buffer
        nFIR_k = new int16_t[n+1]; //filter coefficient array
        //copy in coefficients
        for(int i = 0; i<n; i++){
            x[i] = X[i];
        }
        //float-to-fixed point conversion
        //check if coefficients are over 1,
        for(uint16_t i = 0; i<n; i++){
            if(x[i] > 1.0 || x[i] < -1.0){
                bits -= 1; //will need to re-scale to prevent overflow
                break;
            }
        }
        //convert to fixed point
        for(uint16_t i = 0; i<n; i++){
            nFIR_k[i] = (int16_t)((x[i]*((1<<bits)/2))+0.5);
        }
        nFIR_k[n] = bits-1; //store bit scale value
        N = n; //store N
    }
    //declare constructor variables
    float *x;
    int16_t *buff;
    int16_t *nFIR_k;
    uint8_t bits = 16;
    uint8_t N;
    //FIR Filter Function
    int16_t nFIR(int16_t input){
        //clear previous output
        buff[1] = 0;
        //store indexer
        uint8_t inc = buff[0];
        //store current sample in end of buffer
        buff[1+N-inc] = input;
        //convolve
        for(uint16_t i = 0; i<N; i++){
            buff[1]+= (int16_t)((buff[1+N-inc]*nFIR_k[i])>>nFIR_k[N]);
            if(inc <= 0){inc = N;}
            inc -= 1;
        }
        buff[0] += 1;
        if(buff[0]>=N){buff[0]=0;}
        return buff[1];
    }
};


/* IIR filter, nIIR(input, filter coefficients, order):
 This function implements an IIR filter with user defined filter
 coefficients and order N given in even order only. Function implements
 cascading biquad sections. Filter order MUST BE EVEN.
 Type: float
 ---------------------------------------------------------------
 buff = [current sample, indexer, feedback indexer, output,a1,a2,b0,b1,b2];
 1) Declare filter object: nIIR_Filter LPF(N, fIIR_K) where N is
    the number of filter coefficients and fFIR_K are the filter
    coefficients in form [n1, n2, n3, d1 ,d2] per biQuad
 2) Call function: output = LPF.nIIR(input);
 */

struct nIIR_Filter{
    //constructor (N, filter coefficients)
    nIIR_Filter(uint8_t n, float *X)
    {
        x = new float[(n/2)*5]; //floating point coefficient array
        IIRfSIZE = ((n/2)*5)+1; //fized point array size
        bScale = IIRfSIZE-1;
        buff = new int16_t[9+(n/2)*3]; //buffer
        nIIR_k = new int16_t[IIRfSIZE]; //filter coefficient array
        //copy in coefficients
        for(int i = 0; i<IIRfSIZE-1; i++){
            x[i] = X[i];
        }
        //float-to-fixed point conversion
        //check if coefficients are over 1,
        for(uint16_t i = 0; i<IIRfSIZE-1; i++){
            if(x[i] > 1.0 || x[i] < -1.0){
                bits -= 1; //will need to re-scale to prevent overflow
                break;
            }
        }
        //convert to fixed point
        for(uint16_t i = 0; i<IIRfSIZE-1; i++){
            nIIR_k[i] = (int16_t)((x[i]*((1<<bits)/2))+0.5);
        }
        
        nIIR_k[bScale] = bits-1; //store bit scale value
        nBq = n/2; //store N
    }
    //declare constructor variables
    float *x;
    int16_t *buff;
    int16_t *nIIR_k;
    uint8_t bits = 16;
    uint8_t IIRfSIZE;
    uint8_t nBq;
    uint8_t bScale;
    uint8_t BQ, BQ2, output;

    int16_t nIIR(int16_t input){
        buff[0] = input;
        for(uint8_t j = 0;j<nBq;j++){
            BQ = 6*j;
            BQ2 = 5*j;
            output = 3+BQ;
            //Clear previous output
            buff[output] = 0;
            //store indexer in new variable
            uint16_t inc = buff[1];
            //store current sample at indexer inc
            buff[((8-inc)+BQ)] = buff[0];
        
            //Single BiQuad Section: convolve "FIR" of N=3
            for(uint8_t i = 0; i<3; i++){
                buff[output] += (int16_t)((buff[(8-inc)+BQ]*nIIR_k[i+BQ2])>>nIIR_k[bScale]);
                if(inc <= 0){inc = 3;}
                inc -= 1;
            }
            //Perform feedback operations
            buff[output] -= (int16_t)((buff[(4+buff[2]+BQ)]*nIIR_k[(3+BQ2)])>>nIIR_k[bScale]);
            buff[output] -= (int16_t)((buff[(4+(1^buff[2])+BQ)]*nIIR_k[(4+BQ2)])>>nIIR_k[bScale]);
            //store output in feedback section of buffer
            buff[(4+(1^buff[2])+BQ)] = buff[output];
            buff[0] = buff[output];
        }
    
        //Increment buffers
        buff[2] ^= 1; //toggle location
        buff[1]+=1;
        if(buff[1]>=3){buff[1] = 0;}
        return buff[0];
    }
};


/* float to fixed variable, float2fixed(input, output, bits, N)
 This function converts an array of floating point values to an
 array of fixed point values. For use with IIR and FIR filters.
 Type: float-int16_t
 ---------------------------------------------------------------
 input: array of float values of length N
 output: blank array to store fixed point values after conversion, length N+1
 bits: bits of system, i.e. 8-16-24-32-64
 N: length of array
 */

float float2fixed(float *x, int16_t *y, uint8_t bits, uint16_t uN){
    //check if coefficients are over 1,
    for(uint16_t i = 0; i<uN; i++){
        if(x[i] > 1.0 || x[i] < -1.0){
            bits -= 1; //will need to re-scale to prevent overflow
            break;
        }
    }
    uint32_t scale = (1<<bits)/2;
    for(uint16_t i = 0; i<uN; i++){
        y[i] = (int16_t)((x[i]*scale)+0.5);
    }
    return 0;
}

/* LFO, LFO(Fmin, Fmax, resolution, wavesize, amplitude, sample rate, type)
 This function configures, generates, and implements a user defined LFO with a given frequency range,
 resolution, sample rate, wavesize, and amplitude. It outputs one period of a
 user defined sinwave and a corresponding phase accumulator table for the desired
 frequency range.
 Type: uint16_t
 --------------------------------------------------------------
 F_min = minimum frequnecy desired
 F_max = maximum frequency desired
 resolution = # of steps between min and max frequency
 wavesize = desired length of one period of waveform (even number, must be smaller than amplitude)
 amplitude = desired amplitude of waveform
 Fs = Sampling rate of system
 phase_inc = phase accumulator table/frequency table
 wavetable = one period of waveform
 type: 0 = sine, 1 = triangle, 2 = saw
 
 1) Declare object LFO tremoloLFO(Fmin, Fmax, Fres, wave, amp, Fs, type);
 2) Configure: tremoloLFO.configure();
 3) Call in ISR: tremoloLFO.read(rate, depth)
 */

struct LFO{
    LFO(float fmin, float fmax, uint16_t fres, uint32_t Wavesize, uint16_t amp,
        float fs, uint8_t Type){
        F_min = fmin;
        F_max = fmax;
        F_res = fres;
        wavesize = Wavesize;
        amplitude = amp;
        Fs = fs;
        phase_inc = new uint32_t[fres+2];
        wavetable = new uint16_t[Wavesize+1];
        type = Type;
    }
    float F_min;
    float F_max;
    uint16_t F_res;
    uint32_t wavesize;
    uint16_t amplitude;
    float Fs;
    uint32_t *phase_inc;
    uint16_t *wavetable;
    uint8_t type;
    
    float configure(){
        //calculate absolute frequency range
        float F_range = F_max - F_min;
        //calculate LFO frequency granularity
        float LFO_res = F_range/(F_res-1);
        //generate frequency/phase increment table
        for(int i = 0; i<F_res; i++){
            phase_inc[i+2] = (uint32_t)((F_min + i*LFO_res)*((wavesize << 21)/Fs));
        }
        //Generate step size for trianlge/saw wave
        float tristep = amplitude/(wavesize>>1);
        float sawstep = amplitude/wavesize;
        //generate wavetable
        switch(type){
            case 0: //sine wave
                for(int i = 0; i<wavesize; i++){
                    wavetable[i+1] = (uint16_t)((sin(i*(2*acos(-1))/wavesize)+1)*(amplitude/2)+0.5);
                }
                break;
            case 1: //triangle wave
                for(int i = 0; i<=wavesize>>1; i++){
                    wavetable[i+1] += (uint16_t) tristep*i+0.5;
                }
                for(int i = 1; i<wavesize>>1; i++){
                    wavetable[((wavesize>>1)+1)+i] = (uint16_t)(wavetable[(wavesize>>1)+1]-tristep*i)+0.5;
                }
                break;
            case 2: //saw wave
                for(int i = 0; i<wavesize; i++){
                    wavetable[i+1] = (uint16_t) sawstep*i+0.5;
                }
                
        }
        //store information
        phase_inc[1] = (wavesize<<21);
        wavetable[0] = (uint16_t)amplitude;
        return 0;
    }
    
    /* LFO, read(rate, depth)
     This function simulates an LFO with a pre-configured frequency and wave table.
     The rate controls which phase incrememnt value to be used and the depth scales
     the amplitude of the waveform.
     --------------------------------------------------------------------------------
     1) Call function: "name".read(rate, depth);
     */
    
    uint32_t read(uint16_t rate, uint16_t depth){
        //increment phase
        phase_inc[0] += phase_inc[rate];
        if(phase_inc[0] > phase_inc[1]){
            phase_inc[0] -= phase_inc[1];
        }
        return (wavetable[(phase_inc[0]>>21)+1]*((depth<<16)/wavetable[0]))>>16;
    }
};


/* Chorus Effect, ChorusFX(input, Buffer, n)
 This function used in pair with 'LFO' implements a chorus effect, where
 n defines the max amount of delay possible, (signal to store in memory).
 type: int16_t
 ------------------------------------------------------------------------
 example:
 
 LFO chorusLFO(.03125, 16, 1024, 1000, 1024, SAMPLE_RATE, 0);
 int16_t delayBuff[n+4]; //n must be <= max amplitude of LFO
 
 int main(){
 //configure LFO
 chorusLFO.configure();
 delayBuff[0] = 4;
 }
 
 ISR{
 //get current LFO value and store in delayBuff[2]
 delayBuff[2] = choursLFO.read(rate,depth);
 output = ChorusFX(input,delayBuff, n);
 }
 
 */

int32_t ChorusFX(int16_t in, int16_t *delBuff, int16_t n){
    //delBuff[N+2] = [pointer, delaypointer, LFO, output, 1,2,...N]
    //get location of delayed sample
    delBuff[1] = delBuff[0] - delBuff[2];
    if(delBuff[1] <= 3){delBuff[1] += n;}
    //add input to delayed sample
    delBuff[3] = (int16_t)((in + delBuff[delBuff[1]])>>1);
    //increment buffer
    delBuff[0] += 1;
    //check overflow
    if(delBuff[0] >= n+4){delBuff[0] = 4;}
    //store result into buffer
    delBuff[delBuff[0]] = delBuff[3];
    return delBuff[3];
}

/* 16-bit Compression, dif_compression16(input, memory, memory map)
 This function receives an input audio sample and compresses the audio stream and stores in
 memory. The amount of memory allocated for the compression of the audio stream must be determined 
 beforehand. For example, to store 1 second sampled at 48kHz, the memory buffer must be 4+N/2 where
 N is equal to the number of samples, 48000 in this case. Only multiples of 16 are allowed for this to
 function correctly. This also contains the option to reverse if R_compress is called.
 memory required = 2*(4+N/2) kBytes
 -----------------------------------------------------------------------------------------------
 mem[state,current,previous,strt,n1,n2..nX]
 sign[bit indexer,current,memory indexer, n1,n2,...nX]
 
 1) Declare object: dif_compression16 memory(N);
 2) Call function: memory.compress(input);
 3) Once stored in memory it can be exapnded by: dif_compression.expand()
*/

struct dif_compression16{
    dif_compression16(uint32_t n){
        mem = new int16_t[4+n/2];
        sign = new uint16_t[3+n/16];
        N = n>>1;
    }
    int16_t *mem;
    uint16_t *sign;
    uint32_t N, I;
    //--------------------------------------------------------------------------------------------
    //Normal 16-bit compression
    
    int32_t compress(int16_t input){
        //reduce size and store
        mem[1] = (input>>4);
        if(I > 3){//first sample will be actual value
            //calculate and store difference
            mem[I] = (mem[1]-mem[2]);
            //configure "sign array"
            if(mem[I] < 0){
                sign[1] |= (1<<sign[0]); //set bit x to 1
            }
            else{
                sign[1] &= ~(1<<sign[0]); //set bit x to 0
            }
            sign[0] += 1; //increment
            //check for overflow and store
            if(sign[0] >=16){
                sign[0] = 0; //reset bit indexer
                sign[sign[2]] = sign[1]; //store current 16bit map at memory index
                sign[1] = 0; //clear current
                sign[2] += 1; //increment memory index
            }
            mem[0] ^= 0x01;
            //ready to combine into 16bits
            if(mem[0] == 1){
                //take absolute values of differential pair: 1
                uint16_t temp = mem[I]>>15;
                mem[I]^=temp;
                mem[I]+= temp&1;
                //take absolute values of differential pair: 2
                temp = mem[I-1]>>15;
                mem[I-1]^=temp;
                mem[I-1]+= temp&1;
                //cast into uint8_t, lsb is older sample
                uint8_t lsb = mem[I-1];
                uint8_t msb = mem[I];
                //combine into 16bit
                mem[I-1] = (msb<<8)|(lsb&0xFF);
                I -= 1; //decrement
            }
        }
        else{
            //first sample, dont calculate difference
            mem[3] = mem[1];
            //toggle state
            mem[0] ^= 0x01;
            I = 3; //reset indexer
            sign[2] = 3; //reset memory indexer
            sign[1] = 0;
            sign[0] = 0;
        }
        I+= 1; //increment
        mem[2] = mem[1]; //store last sample
        if(I >= N+4){
            I = 0;//reset
            sign[sign[2]] = sign[1]; //store current 16bit map at memory index
        }
        return 0;
    }
    
    //---------------------------------------------------------------------------
    //Reverse Memory Compression
    
    int32_t R_compress(int16_t input){
        mem[1] = (input>>4); //reduce size and store
        //initialize
        if(I==0){
            I = N+4; //reset indexer
            sign[2] = 2+((N<<1)>>4); //reset memory indexer
            sign[0] = 14; //reset bit indexer
        }
        else{
            mem[I] = (mem[1]-mem[2]); //calculate and store difference
            //configure "sign array"
            if(mem[I] < 0){
                sign[1] &= ~(1<<sign[0]); //set bit x to 0
            }
            else{
                sign[1] |= (1<<sign[0]); //set bit x to 1
            }
            sign[0] -= 1; //increment
            //check for overflow and store
            if(sign[0] >=16){
                sign[0] = 15; //reset bit indexer
                sign[sign[2]] = sign[1]; //store current 16bit map at memory index
                sign[1] = 0; //clear current
                sign[2] -= 1; //increment memory index
            }
            mem[0] ^= 0x01;
            //ready to combine into 16bits
            if(mem[0] == 1 && I<(N+2)){
                //take absolute values of differential pair: 1
                uint16_t temp = mem[I]>>15;
                mem[I]^=temp;
                mem[I]+= temp&1;
                //take absolute values of differential pair: 2
                temp = mem[I+1]>>15;
                mem[I+1]^=temp;
                mem[I+1]+= temp&1;
                //cast into uint8_t, lsb is older sample
                uint8_t lsb = mem[I];
                uint8_t msb = mem[I+1];
                mem[I+1] = (msb<<8)|(lsb&0xFF); //combine into 16bit
                I += 1; //decrement
            }
        }
        I-=1; //decrement
        mem[2] = mem[1]; //store last sample
        if(I == 3 && mem[0] == 1){
            mem[3] = mem[2];
            I = 0;
            mem[0] ^= 0x01;
            sign[sign[2]] = sign[1];
        }
        return 0;
    }
    
    /* 16-bit Expansion, dif_expansion16(memory, memory map)
     This function takes the memory and memory map buffers generated by dif_compression16()
     and de-compresses them to their original form.
     --------------------------------------------------------------------------------------
     1) Call Function (comp_mem, mem_map);
     */
    
    int32_t expand(){
        //If first sample is established..
        if(I > 3){
            //toggle state now defines whether to extract high or low byte
            if(mem[0] == 1){
                //LSB 1st
                mem[1] = (mem[I])&0xFF;
            }
            else{
                //MSB 2nd
                mem[1] = (uint8_t)(256+(mem[I]>>8));
                I += 1; //increment
            }
            mem[0] ^= 0x01; //toggle
            //Now apply correct sign to differential value
            if(((sign[sign[2]]>>sign[0])&1)==1){
                mem[1] *= -1;
            }
            sign[0] += 1; //increment
            //check for overflow and store
            if(sign[0] >=16){
                sign[0] = 0; //reset bit indexer
                sign[2] += 1; //increment memory index
            }
            //add result to previous
            mem[1] += mem[2];
        }
        else{
            mem[0] ^= 0x01; //toggle
            I = 4; //reset indexer
            sign[2] = 3; //reset memory indexer
            sign[0] = 0; //reset bit indexer
            mem[1] = mem[3];
        }
        mem[2] = mem[1]; //store last sample
        if(I >= N+4){I = 0;}//reset
        return (mem[1]<<4)+16;
    }
};


#endif
