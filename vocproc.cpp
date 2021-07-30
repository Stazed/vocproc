/* 
    LV2 plugin for pitch shifting, pitch correction, vocoding and harmonizing
    singing voice

    (c) Igor Brkic 2010

    Phase vocoder code adapted from:
    http://www.dspdimension.com/admin/pitch-shifting-using-the-ft/

    
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/


#include <lv2plugin.hpp>
#include "vocproc.peg"
#include <fftw3.h>
#include <string.h>     // because of memset
#include <stdlib.h>
#include <math.h>

using namespace LV2;

#define M_2PI 6.2831853071795
#define MAX_FRAME_LENGTH 4096

#define ROUND(a) ((float)((int)(a+0.5)))

class VocProc : public Plugin<VocProc> {
    private:
        float   fSamplingFreq;
        float 	sPitchFactor, sEffect, sOutputGain;
        float   cFormantVoco, cEffect;
        float   powerIn;
        float   sSwitch;
        float   pOffset[2];

        float   cAutoTune;

        float   *workspace, *gInFIFO, *gIn2FIFO, *gOutFIFO, *gOutputAccum;
        float   *window;

        long    fftFrameSize, overlap;

        float   freqOld;

        double  *fftTmpR;
        fftw_complex *fftTmpC;
        fftw_complex *fftOldC;
        fftw_complex *fftCeps;
        fftw_plan fftPlan;

        void    spectralEnvelope(float *env, fftw_complex *fft, uint32_t nframes);

        float   pitchFrequency(fftw_complex *block);
        void    setPitchFactor(float freq);

        void    phaseVocAnalysis(fftw_complex *block, float *gLastPhase, double freqPerBin, double expct, float *gAnaMagn, float *gAnaFreq);
        void    phaseVocSynthesis(fftw_complex *block, float *gSumPhase, float *gSynMagn, float *gSynFreq, double freqPerBin, double expct);



    public:

        VocProc(double rate);
        ~VocProc();

        void run(uint32_t nframes);
};

static int _ = VocProc::register_class("http://hyperglitch.com/dev/VocProc");


// initialization
VocProc::VocProc(double rate) : Plugin<VocProc>(p_n_ports){
    fSamplingFreq = (float)rate;

    sPitchFactor = 1.0;

    sEffect=0.0;
    cEffect=0.0;

    sOutputGain=1.0;

    cFormantVoco=1;

    pOffset[0]=0; pOffset[1]=0;

    sSwitch=0.0;

    cAutoTune=0.0;

    powerIn=0;

    fftFrameSize=2048;  // pitch detection currently doesn't work for 1024
                        // and there is aliasing present for 1024 with formant correction when
                        // large pitch shifting factor is applied
    overlap=4;

    freqOld=0;

    window=(float*)malloc(fftFrameSize*sizeof(float));
    for(int k=0;k<fftFrameSize;k++)
        window[k] = -.5*cos(M_2PI*(float)k/(float)fftFrameSize)+.5;

    gInFIFO=(float*)calloc(fftFrameSize, sizeof(float));
#ifndef NO_VOCODER
    gIn2FIFO=(float*)calloc(fftFrameSize, sizeof(float));
#endif
    gOutFIFO=(float*)calloc(fftFrameSize, sizeof(float));
    gOutputAccum=(float*)calloc(2*fftFrameSize, sizeof(float));

    // FFTW stuff
    fftTmpR=(double*)fftw_malloc(sizeof(double)*fftFrameSize);
    fftTmpC=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*fftFrameSize);
    fftOldC=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*fftFrameSize);
    fftCeps=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*fftFrameSize);

}
VocProc::~VocProc(){

    free(gInFIFO);
#ifndef NO_VOCODER
    free(gIn2FIFO);
#endif
    free(gOutFIFO);
    free(gOutputAccum);

    fftw_free(fftTmpR);
    fftw_free(fftTmpC);
    fftw_free(fftOldC);
    fftw_free(fftCeps);
};


/*************************************************
    here be procezzing
*************************************************/

void VocProc::run(uint32_t nframes) {

    float *input0 = p(p_voice);

#ifndef NO_VOCODER
    float *input1 = p(p_carrier);
#endif

    float *output0 = p(p_output);

    const float conversionTable[]={
        0.500000, 0.529731, 0.561231, 0.594603, 0.629960, 0.667419, 0.707106, 0.749153, 0.793700, 0.840896, 0.890898, 0.943874, 
        1.000000, 1.059463, 1.122462, 1.189207, 1.259921, 1.334840, 1.414214, 1.498307, 1.587401, 1.681793, 1.781797, 1.887749, 2.000000
    };

    float freqScaling=*p(p_pitch_factor);
    if(freqScaling<-12) freqScaling=-12;
    if(freqScaling>12) freqScaling=12;
    sPitchFactor = conversionTable[(int)(freqScaling+0.5)+12];

    if(*p(p_effect)==0) cEffect=0;
    else{
        cEffect=1;
        sEffect=*p(p_effect);
    }
    cFormantVoco=(int)(*p(p_fc_voc_switch)+0.5);

#ifdef NO_VOCODER
    sSwitch=0.0;
#else
    switch((int)(*p(p_fc_voc)+0.5)){
        case 0:
            sSwitch=0.0;
            break;
        case 1:
            sSwitch=1.0;
            break;
    }
#endif

    cAutoTune=ROUND(*p(p_pitch_correction));


    static float gLastPhase[MAX_FRAME_LENGTH/2+1];
    static float gSumPhase[MAX_FRAME_LENGTH/2+1];
    static float gAnaFreq[MAX_FRAME_LENGTH];
    static float gAnaMagn[MAX_FRAME_LENGTH];
    static float gSynFreq[MAX_FRAME_LENGTH];
    static float gSynMagn[MAX_FRAME_LENGTH];

    static long gRover = false, gInit = false;

    double freqPerBin, expct;
    long i,k, index, inFifoLatency, stepSize, fftFrameSize2;

    float *fPointer, *fPointer2, *fPointer3;
    double *dPointer;

    /* set up some handy variables */
    fftFrameSize2 = fftFrameSize/2;

    stepSize = fftFrameSize/overlap;
    freqPerBin = (double)fSamplingFreq/(double)fftFrameSize;
    expct = M_2PI*(double)stepSize/(double)fftFrameSize;

    inFifoLatency = fftFrameSize-stepSize;
    if (gRover == false) gRover = inFifoLatency;

    /* initialize our static arrays */
    if (gInit == false) {
        memset(gLastPhase, 0, (MAX_FRAME_LENGTH/2+1)*sizeof(float));
        memset(gSumPhase, 0, (MAX_FRAME_LENGTH/2+1)*sizeof(float));
        memset(gAnaFreq, 0, MAX_FRAME_LENGTH*sizeof(float));
        memset(gAnaMagn, 0, MAX_FRAME_LENGTH*sizeof(float));
        gInit = true;
    }

    /* main processing loop */
    for (i = 0; i < nframes; i++){

        // As long as we have not yet collected enough data just read in
        gInFIFO[gRover] = input0[i];

#ifndef NO_VOCODER
        gIn2FIFO[gRover] = input1[i];
#endif

        output0[i] = gOutFIFO[gRover-inFifoLatency];
        gRover++;

        // now we have enough data for processing
        if (gRover >= fftFrameSize) {
            gRover = inFifoLatency;

            float tmpPower=0.0;
            dPointer=fftTmpR;
            fPointer=gInFIFO;
            fPointer2=window;
            for (k = 0; k < fftFrameSize;k++) {
                *dPointer=*(fPointer++) * *(fPointer2++);
                tmpPower+= *dPointer * *dPointer;
                dPointer++;
            }

            powerIn=tmpPower/(float)fftFrameSize;

            // do transform
            fftPlan=fftw_plan_dft_r2c_1d(fftFrameSize, fftTmpR, fftTmpC, FFTW_ESTIMATE);
            fftw_execute(fftPlan);
            fftw_destroy_plan(fftPlan);

            // pitch correction
            if(cAutoTune){
                float freq;
                freq=pitchFrequency(fftTmpC);
                setPitchFactor(freq);
            }

            memcpy(fftOldC, fftTmpC, fftFrameSize*sizeof(fftw_complex));

            // pitch shifting with phase vocoder
            phaseVocAnalysis(fftTmpC, gLastPhase, freqPerBin, expct, gAnaMagn, gAnaFreq);
            memset(gSynMagn, 0, fftFrameSize*sizeof(float));
            memset(gSynFreq, 0, fftFrameSize*sizeof(float));
            for (k = 0; k <= fftFrameSize2; k++) {
                index = k*sPitchFactor;
                if (index <= fftFrameSize2) {
                    gSynMagn[index] += gAnaMagn[k];
                    gSynFreq[index] = gAnaFreq[k] * sPitchFactor;
                    if(cEffect)
                        gSynFreq[index] = gSynFreq[index]*sEffect + sEffect*200*(float)rand()/RAND_MAX-100;
                }
            }
            phaseVocSynthesis(fftTmpC, gSumPhase, gSynMagn, gSynFreq, freqPerBin, expct);

            // formant correction + vocoder
            if(cFormantVoco){
                float env1[fftFrameSize2], env2[fftFrameSize2];

#ifndef NO_VOCODER
                if(sSwitch){
                    dPointer=fftTmpR; fPointer=gIn2FIFO; fPointer2=window;
                    for (k = 0; k < fftFrameSize;k++) {
                        *dPointer++=*(fPointer++) * *(fPointer2++);
                    }

                    fftPlan=fftw_plan_dft_r2c_1d(fftFrameSize, fftTmpR, fftTmpC, FFTW_ESTIMATE);
                    fftw_execute(fftPlan);
                    fftw_destroy_plan(fftPlan);
                }
#endif

                spectralEnvelope(env1, fftOldC, fftFrameSize2);
                spectralEnvelope(env2, fftTmpC, fftFrameSize2);

                // modify spectral envelope of spectrum in fftTmpC to look like spectral
                // envelope of spectrum in fftOldC
                float koef;
                fPointer2=env1; fPointer3=env2;
                for(k=0;k<fftFrameSize2;k++){
                    koef = *(fPointer2++) / (*(fPointer3++)+.02)*2;
                    fftTmpC[k][0] *= koef;
                    fftTmpC[k][1] *= koef;
                }
            }

            // do inverse transform
            fftPlan=fftw_plan_dft_c2r_1d(fftFrameSize, fftTmpC, fftTmpR, FFTW_ESTIMATE);
            fftw_execute(fftPlan);
            fftw_destroy_plan(fftPlan);

            fPointer=gOutputAccum; dPointer=fftTmpR; fPointer2=window;
            for(k=0; k < fftFrameSize; k++) {
                *fPointer += 0.7 * *(dPointer++) / (fftFrameSize2*overlap) * sOutputGain * *(fPointer2++);
                fPointer++;
            }

            memcpy(gOutFIFO, gOutputAccum, stepSize*sizeof(float));

            // shift accumulator
            memmove(gOutputAccum, gOutputAccum+stepSize, fftFrameSize*sizeof(float));

            // move input FIFO
            for (k = 0; k < inFifoLatency; k++) gInFIFO[k] = gInFIFO[k+stepSize];

#ifndef NO_VOCODER
            for (k = 0; k < inFifoLatency; k++) gIn2FIFO[k] = gIn2FIFO[k+stepSize];
#endif
        }
    }
}


void VocProc::phaseVocAnalysis(fftw_complex *block, float *gLastPhase, double freqPerBin, double expct, float *gAnaMagn, float *gAnaFreq){

    double real, imag, magn, phase, tmp;
    long qpd, k;

    for (k = 0; k <= fftFrameSize/2; k++) {

        /* de-interlace FFT buffer */
        real = block[k][0];
        imag = block[k][1];

        /* compute magnitude and phase */
        magn = 2.*sqrt(real*real + imag*imag);
        phase = atan2(imag,real);

        /* compute phase difference */
        tmp = phase - gLastPhase[k];
        gLastPhase[k] = phase;

        /* subtract expected phase difference */
        tmp -= (double)k*expct;

        /* map delta phase into +/- Pi interval */
        qpd = tmp/M_PI;
        if (qpd >= 0) qpd += qpd&1;
        else qpd -= qpd&1;
        tmp -= M_PI*(double)qpd;

        /* get deviation from bin frequency from the +/- Pi interval */
        tmp = overlap*tmp/(M_2PI);

        /* compute the k-th partials' true frequency */
        tmp = (double)k*freqPerBin + tmp*freqPerBin;

        /* store magnitude and true frequency in analysis arrays */
        gAnaMagn[k] = magn;
        gAnaFreq[k] = tmp;
    }
}

void VocProc::phaseVocSynthesis(fftw_complex *block, float *gSumPhase, float *gSynMagn, float *gSynFreq, double freqPerBin, double expct){

    int k;
    double magn, tmp, phase;

    for (k = 0; k <= fftFrameSize/2; k++) {
        /* get magnitude and true frequency from synthesis arrays */
        magn = gSynMagn[k];
        tmp = gSynFreq[k];

        /* subtract bin mid frequency */
        tmp -= (double)k*freqPerBin;

        /* get bin deviation from freq deviation */
        tmp /= freqPerBin;

        /* take overlap into acnframes */
        tmp = M_2PI*tmp/overlap;

        /* add the overlap phase advance back in */
        tmp += (double)k*expct;

        /* accumulate delta phase to get bin phase */
        gSumPhase[k] += tmp;
        phase = gSumPhase[k];

        /* get real and imag part and re-interleave */
        block[k][0] = magn*cos(phase);
        block[k][1] = magn*sin(phase);
    }
}

void VocProc::spectralEnvelope(float *env, fftw_complex *fft, uint32_t nframes){

    int nTaps=20;
    int nTaps2=10;
    float tmp[nframes+nTaps];

    float h[]={ // h=(firls(20, [0 0.02 0.1 1], [1 1 0 0]));
        0.0180, 0.0243, 0.0310, 0.0378, 0.0445, 0.0508, 0.0564, 0.0611,
        0.0646, 0.0667, 0.0675, 0.0667, 0.0646, 0.0611, 0.0564, 0.0508,
        0.0445, 0.0378, 0.0310, 0.0243, 0.0180
    };

    // |H(w)|
    memset(tmp,  0, (nframes+nTaps)*sizeof(float));
    for(int k=0;k<nframes;k++)
        tmp[k]=sqrt(fft[k][0]*fft[k][0]+fft[k][1]*fft[k][1]);

    memset(env, 0, nframes*sizeof(float));

    // magnitude spectrum filtering
    int i, j;
    float *p_h, *p_z, accum;

    float z[2 * nTaps];
    memset(z, 0, 2*nTaps*sizeof(float));
    int state = 0;
    for (j = 0; j < nframes+nTaps2; j++) {
        z[state] = z[state + nTaps] = tmp[j];
        p_h = h;
        p_z = z + state;
        accum = 0;
        for (i = 0; i < nTaps; i++) accum += *p_h++ * *p_z++;
        if (--state < 0) state += nTaps;
        if(j>=nTaps2) env[j-nTaps2]=accum;
    }
}

float VocProc::pitchFrequency(fftw_complex *block){
    // cepstral method - needs improvement but good for now

    float freq;
    int i;
    double cepst[fftFrameSize/2];
    int ppMin, ppMax;
    double max1;
    float pitch=0;

    //double fftTmpReal[fftFrameSize];
    //fftw_complex fftTmpCplx[fftFrameSize];

    for(i=0;i<fftFrameSize/2;i++){
        fftCeps[i][0]=log(sqrt(pow(block[i][0], 2)+pow(block[i][1], 2))+1e-6)/(float)fftFrameSize;
        fftCeps[i][1]=0;
    }

    fftPlan=fftw_plan_dft_c2r_1d(fftFrameSize, fftCeps, fftTmpR, FFTW_ESTIMATE);
    fftw_execute(fftPlan);
    fftw_destroy_plan(fftPlan);

    // normalize
    for(i=0;i<fftFrameSize/2;i++) cepst[i]=fabs(fftTmpR[i]/(float)fftFrameSize)+1e6;

    ppMax=fftFrameSize/2-2;
    ppMin=fSamplingFreq/1200; // fMax=1200

    // find maximum of HTP cepstrum
    max1=0;
    for(i=ppMin;i<=ppMax;i++){
        if(cepst[i]>max1){
            max1=cepst[i];
            pitch=i;
        }
    }

    // interpolate between two samples
    int idx=pitch;
    int l1;

    if(cepst[idx-1]>cepst[idx+1]) l1=pitch-1;
    else l1=pitch;

    pitch=l1+1/(cepst[l1]/cepst[l1+1]+1);

    freq=fSamplingFreq/pitch;

    return freq;
}

void VocProc::setPitchFactor(float freq){
    // find nearest tone from chosen scale

    float scale[14]; // full chromatic + edges
    int scaleLen;
    float mult;

    // fill scale
    scaleLen=1;
    if(*p(p_c )==1) scale[scaleLen++]=130.813;
    if(*p(p_cc)==1) scale[scaleLen++]=138.591;
    if(*p(p_d )==1) scale[scaleLen++]=146.832;
    if(*p(p_dd)==1) scale[scaleLen++]=155.563;
    if(*p(p_e )==1) scale[scaleLen++]=164.814;
    if(*p(p_f )==1) scale[scaleLen++]=174.614;
    if(*p(p_ff)==1) scale[scaleLen++]=184.997;
    if(*p(p_g )==1) scale[scaleLen++]=195.998;
    if(*p(p_gg)==1) scale[scaleLen++]=207.652;
    if(*p(p_a )==1) scale[scaleLen++]=220.000;
    if(*p(p_aa)==1) scale[scaleLen++]=233.082;
    if(*p(p_b )==1) scale[scaleLen++]=246.942;

    if(scaleLen==1){
        // no tones have been chosen - do nothing
        sPitchFactor=1;
    }
    else{
        // add edges
        scale[0]=scale[scaleLen-1]/2.0;
        scale[scaleLen]=scale[1]*2.0;

        // calculate multiplicator
        if(scale[scaleLen-1]<freq)
            mult=(int)(freq/scale[scaleLen-1])+1.0;
        else if(scale[1]>freq)
            mult=1.0/((int)(scale[1]/freq)+1.0);
        else
            mult=1;


        float *curr=NULL;
        int i;
      
        for(i=1;i<=scaleLen;i++){
            curr=&scale[i];
            if(freq<scale[i]) break;
        }

        // add transposition if needed (and if possible)
        if((i+*p(p_transpose))<=scaleLen && (i+*p(p_transpose))>=0){
            curr+=(int)(*p(p_transpose)+0.5);
        }

        // add hysteresis to reduce oscillation around middle point

        float sign;
        if((freqOld-freq)>0) sign=-1; // pitch is falling
        else sign=1; // pitch is rising

        float thr=(*(curr-1)+*curr)/2 + 0.3*sign*(*curr-*(curr-1));

        // if current frequency is lower than threshold use lower note
        if(freq<thr) curr--;

        // calculate pitch factor and do pitch factor averaging (attack)
        float factor=sPitchFactor;
        factor*=(1+(float)((int)(*p(p_attack)*20.0)));
        factor+=(*curr/freq);
        factor/=((float)((int)(*p(p_attack)*20.0))+2.0);
        
        // calculate offset (in cents)
        float offset=3986*log10(factor);
        if(offset<-100) offset=-100;
        if(offset>100) offset=100;

        if(powerIn<0.001) offset=0;

        // median (this probably could be much nicer written)
        float tmpSort[3];
        tmpSort[0]=pOffset[0];
        tmpSort[1]=pOffset[1];
        tmpSort[2]=offset;

        float tmp;
        if(tmpSort[0]>tmpSort[1]){
            tmp=tmpSort[0];
            tmpSort[0]=tmpSort[1];
            tmpSort[1]=tmp;
        }
        if(tmpSort[0]>tmpSort[2]){
            tmp=tmpSort[0];
            tmpSort[0]=tmpSort[2];
            tmpSort[2]=tmp;
        }
        if(tmpSort[1]>tmpSort[2]){
            tmp=tmpSort[1];
            tmpSort[1]=tmpSort[2];
            tmpSort[2]=tmp;
        }

        *p(p_offset)=tmpSort[1];

        pOffset[0]=pOffset[1];
        pOffset[1]=offset;

        // see if correction is needed at all
        if((fabs(*curr-freq)/freq)>*p(p_threshold)*0.07)
            sPitchFactor=factor;
        else
            sPitchFactor=1;

        // honor limits
        if(sPitchFactor>2.0 || sPitchFactor<0.5) sPitchFactor=1.0;
    }
}


