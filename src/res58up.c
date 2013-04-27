    /* resamplesubs.c - sampling rate conversion subroutines */
// Altered version
#include "resample.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define IBUFFSIZE 4096 /* Input buffer size */
#define BBDOFF 100    /*padding in the BBD array */
#define BBDBUCKS 512  /* Number of buckets/2 */
#define BBDDEL (BBDBUCKS/2)
#define BBDSIZE (BBDDEL+2*BBDOFF) /* Size of the BBD array*/
#define TSIZE 100

#include "smallfilter.h"
#include "largefilter.h"

#include "filterkit.h"
#include "sndlibextra.h"


static int SrcUDold(FLOAT X[], FLOAT Y[], double factor, UWORD *Time,
                 UHWORD Nx, UHWORD Nwing, UHWORD LpScl,
                 HWORD Imp[], HWORD ImpD[], BOOL Interp, FLOAT X2[], int preset, int channel, double rescale, BOOL onboard);
                 


double custom1;
double custom2;
double custom3;
FLOAT inmin=1;
FLOAT inmax=-1;
FLOAT outmin=1;
FLOAT outmax=-1;
int mindiff=BBDSIZE;
int maxdiff=0;
/* CAUTION: Assumes we call this for only one resample job per program run! */
/* return: 0 - notDone */
/*        >0 - index of last sample */
static int
readData(int   infd,          /* input file descriptor */
         int   inCount,       /* _total_ number of frames in input file */
         FLOAT *outPtr1,      /* array receiving left chan samps */
         FLOAT *outPtr2,      /* array receiving right chan samps */
         int   dataArraySize, /* size of these arrays */
         int   nChans,
         int   Xoff)          /* read into input array starting at this index */
{
   int    i, Nsamps, nret;
   static unsigned int framecount;  /* frames previously read */
   static mus_sample_t **ibufs = NULL;

   if (ibufs == NULL) {             /* first time called, so allocate it */
      ibufs = sndlib_allocate_buffers(nChans, dataArraySize);
      if (ibufs == NULL) {
         fprintf(stderr, "readData: Can't allocate input buffers!\n");
         exit(1);
      }
      framecount = 0;               /* init this too */
   }

   Nsamps = dataArraySize - Xoff;   /* Calculate number of samples to get */
   outPtr1 += Xoff;                 /* Start at designated sample number */
   outPtr2 += Xoff;

   nret = mus_file_read(infd, 0, Nsamps - 1, nChans, ibufs);
   if (nret < 0) {
     fprintf(stderr, "readData: Can't read data!\n");
     exit(1);
   }     

   /* NB: sndlib pads ibufs with zeros if it reads past EOF. */
   if (nChans == 1) {
      for (i = 0; i < Nsamps; i++)
         *outPtr1++ = MUS_SAMPLE_TYPE_TO_FLOAT(ibufs[0][i]);
   }
   else {
      for (i = 0; i < Nsamps; i++) {
         *outPtr1++ = MUS_SAMPLE_TYPE_TO_FLOAT(ibufs[0][i]);
         *outPtr2++ = MUS_SAMPLE_TYPE_TO_FLOAT(ibufs[1][i]);
      }
   }

   framecount += Nsamps;

   if (framecount >= (unsigned)inCount)     /* return index of last samp */
      return (((Nsamps - (framecount - inCount)) - 1) + Xoff);
   else
      return 0;
}


#ifdef DEBUG
static int pof = 0;             /* positive overflow count */
static int nof = 0;             /* negative overflow count */
#endif

static INLINE HWORD WordToHword(WORD v, int scl)
{
    HWORD out;
    WORD llsb = (1<<(scl-1));
    v += llsb;          /* round */
    v >>= scl;
    if (v>MAX_HWORD) {
#ifdef DEBUG
        if (pof == 0)
          fprintf(stderr, "*** resample: sound sample overflow\n");
        else if ((pof % 10000) == 0)
          fprintf(stderr, "*** resample: another ten thousand overflows\n");
        pof++;
#endif
        v = MAX_HWORD;
    } else if (v < MIN_HWORD) {
#ifdef DEBUG
        if (nof == 0)
          fprintf(stderr, "*** resample: sound sample (-) overflow\n");
        else if ((nof % 1000) == 0)
          fprintf(stderr, "*** resample: another thousand (-) overflows\n");
        nof++;
#endif
        v = MIN_HWORD;
    }   
    out = (HWORD) v;
    return out;
}
static UWORD next(double T[],UWORD offset,int TLen)
{
  int i;
  //~ for(i=0; i<TLen-6; i+=6) printf ("%10.7f %10.7f %10.7f %10.7f %10.7f %10.7f\n", T[i],T[i+1],T[i+2], T[i+3],T[i+4],T[i+5]);
  int inc=1<<Np;
  UWORD rread=offset;
  i=(offset>>Np);
  while (inc>0 && i<TLen)
  {
    if ( (rread&Pmask) + ( inc / T[i]) < (1<<Np) )
    {
      rread+=( inc / T[i]);
      break;
    }
    inc-=( (1<<Np) - (rread&Pmask) )*T[i];
    rread=((rread>>Np) + 1) <<Np;
    i++;
    i%=TLen;
  }
  return rread-offset;
}

void raw2(int input)
{
    input+=32768;
    input/=2;
    if (input<0)
        printf("errore");
    printf("%c%c",(int)(input/256),input%256);
}
void raw2e(int input)
{
    input+=32768;
    input/=2;
    if (input<0)
    {    printf("errore");
        fprintf(stderr,"errore");
    }
    fprintf(stderr,"%c%c",(int)(input/256),input%256);
}
void raw4(int input){
    raw2((int)(input/pow(2,16)));
    raw2((int)(input%(int)pow(2,16)));
}

void raw(int input)
{
    printf("%c%c",input/256,input%256);
}
void rawe(int input)
{
    fprintf(stderr,"%c%c",input/256,input%256);
}

#define Blgt 5
#define Algt 5
double A[Algt]={1.000000,-2.099811,2.439923,-1.502647,0.483548};
double B[Blgt]={0.010773,0.043093,0.064640,0.043093,0.010773};
FLOAT filtiir(FLOAT x[],int lengthx, FLOAT y[]){
    int n,i,iend;
    double out=0;
    for (n=0;n<lengthx;n++){
        out=0;
        iend=Blgt;
        for(i=0;i<iend;i++){
            out+=B[i]*x[n-i];
        }
        iend=Algt;
        for (i=1;i<iend;i++){
            out-=A[i]*y[n-i];
        }
        y[n]=out;
    }
}

double fslow=.83;/*frequency of the slow waveform*/
double ffast=6.16;/*frequency of the fast waveform*/

double phs=0;
double phf=0;

UHWORD Xoff;

void reinitLfoGen(UWORD Time){
  double t=(Time/32768.0 - Xoff)/44100.0;
  phs+= Xoff/44100 + t*fslow*2*M_PI;
  phs=fmod(phs,2*M_PI);
  phf+= (t*ffast)*2*M_PI;
  phf=fmod(phf,2*M_PI);
}

double lfoGen(int preset, int channel, UWORD uTime){
    double Time=uTime/32768.0;
    static double pifract=2*M_PI/44100.0;
    static double ph1=2*M_PI/3;/*phase of second modulator*/
    static double ph2=2*M_PI*2/3;/*phase of third modulator*/
    double samp=1;/*gain reduction of slow waveform*/
    double famp=1*68.0/180.0; /*gain reduction of fast waveform*/
    switch (preset){
    case 0 : /*slow*/
        switch (channel){
        case 0 :
            return samp*sin(pifract*Time*fslow+phs);
        case 1 : 
            return 0;
        case 2 :
            return 0;
        }
    case 1 : /*chorus 1*/
        switch (channel){
        case 0:
            return samp*sin(pifract*Time*fslow+phs);
        case 1:
            return samp*sin(pifract*Time*fslow+phs+ph1);
        case 2:
            return samp*sin(pifract*Time*fslow+phs+ph2);
        }
    case 2 : /*chorus 2*/
        switch (channel){
        case 0 :
            return .5*sin(pifract*Time*fslow+phs) + famp/2*sin(pifract*Time*ffast+phf);
        case 1 :
            return .5*sin(pifract*Time*fslow+phs+ph1) + famp/2*sin(pifract*Time*ffast+phf+ph1);
        case 2 : 
            return .5*sin(pifract*Time*fslow+phs+ph2) + famp/2*sin(pifract*Time*ffast+phf+ph2);
        }
    case 3: /*fast*/
        switch (channel){
        case 0:
            return sin(pifract*ffast*Time+phf);
        case 1: 
            return 0;
        case 2:
            return 0;
        }
    }
    fprintf(stdout,"unknown preset number in lfoGen");
    exit(-1);
}

/* Sampling rate conversion subroutine */
FLOAT bbd[3][BBDSIZE]={0};
double bbdT[TSIZE];

static double wash(double v){
  return  v*v*v*v*v*v*1.2722+
          v*v*v*v*v*0.62791+
          v*v*v*v*-1.0484+
          v*v*v*-0.24642+
          v*v*0.17149+
          v*0.92521+
          -0.047669;  
  return  v*v*v*v*v*v*1.4722+
          v*v*v*v*v*0.82791+
          v*v*v*v*-1.2484+
          v*v*v*-0.44642+
          v*v*0.27149+
          v*0.92521+
          -0.047669;
}

static int SrcUD(FLOAT X[], FLOAT Y[], double factor, UWORD *Time,
                 UHWORD Nx, UHWORD Nwing, UHWORD LpScl,
                 HWORD Imp[], HWORD ImpD[], BOOL Interp, FLOAT X2[], int preset, int channels, double rescale, BOOL onboard)
{ 
  static BOOL init =0;
  int bbdp=BBDOFF, bbdsize=BBDSIZE, bbddel=BBDDEL;/*debug*/
  //~ T2b=.7;
  //~ T2=T2b;
  double z0,z1,w;//used for linear interpolation
  static double fs1=44100;
  double fs2=100000;
  double T2f=fs1/fs2;//step through bbd samples, floating point
  UWORD T1=1<<Np; //step through input samples, fixed to 1 sample
  UWORD T2i=(UWORD)(T2f*(1<<Np)+.5);//step through bbd samples, fixed point
  UWORD T2[3]={T2i,T2i,T2i};
  double width=custom1; //percentage of predelay
  double lfo_width=T2f*width/100;
  static int bbdpointer[3]={0};
  int ch,read;
  double lfo=0;
  //time variables
  static UWORD oldBBDevent[3];//initialized below
  static UWORD nextBBDevent[3];//initialized below
  static UWORD nextTevent[3];//initialized below
  if (!init){//first time initialization
     for (ch=0;ch<channels;ch++){
      oldBBDevent[ch]+=*Time-T2[0];
      nextBBDevent[ch]+=*Time;
      nextTevent[ch]=*Time;
    }
    init=1;
  } 
  else{ //reinizialization
    for (ch=0;ch<channels;ch++){
      oldBBDevent[ch]+=*Time-nextTevent[ch];
      nextBBDevent[ch]+=*Time-nextTevent[ch];
      nextTevent[ch]=*Time;
    }
  }
  UWORD startTime=*Time;
  UWORD newTime;
  FLOAT *Ystart = Y;

  while ( *Time-startTime >>Np < Nx ){
    newTime=*Time+(1<<Np);
    for (ch=0;ch<channels;ch++){
      if (*Time==nextBBDevent[ch]){
        lfo=lfoGen(preset,ch,*Time);
        T2[ch]=(UWORD)( (T2f + lfo_width*lfo)*(1<<Np) +.5);
        //~ printf("%.3f write bbd\n",*Time/32768.0 );
        bbdpointer[ch]=(bbdpointer[ch]+1)%bbdsize; //increment in the circular buffer
        z0=X[*Time>>Np];z1=X[1+(*Time>>Np)];w=(*Time&Pmask)/32768.0;
        bbd[ch][bbdpointer[ch]] = z0+(z1-z0)*w;//linear interpolation
        //~ bbd[ch][bbdpointer]=1;
        oldBBDevent[ch]=*Time;
        nextBBDevent[ch]=*Time + T2[ch];
      }
      newTime=newTime<nextBBDevent[ch]?newTime:nextBBDevent[ch];
      if (*Time==nextTevent[ch]){
        //~ printf("%.3f read bbd\n",*Time/32768.0);
        read=(bbdpointer[ch]-bbddel+bbdsize-1)%bbdsize;
        z0=bbd[ch][ read ];
        z1=bbd[ch][ (read+1)%bbdsize];
        w=(*Time-oldBBDevent[ch])/(double)(T2[ch]);
        *Y+=(z0+(z1-z0)*w)/rescale; //linear interpolation
        //~ *Y=lfo; //linear interpolation
        nextTevent[ch]=*Time+(1<<Np); //integer
        if (ch==channels-1) Y++;
      }
      newTime=newTime<nextTevent[ch]?newTime:nextTevent[ch];
    }
    *Time=newTime;/*min(nextBBDevent[],nextTevent[])*/
  }
  reinitLfoGen(*Time);
  return (Y - Ystart);        /* Return the number of output samples */
}



static int SrcUDold2(FLOAT X[], FLOAT Y[], double factor, UWORD *Time,
                 UHWORD Nx, UHWORD Nwing, UHWORD LpScl,
                 HWORD Imp[], HWORD ImpD[], BOOL Interp, FLOAT X2[], int preset, int channels, double rescale, BOOL onboard)
{ 
  static BOOL init =0;
  int bbdp=BBDOFF, bbdsize=BBDSIZE, bbddel=BBDDEL;/*debug*/
  double z0,z1,w;//used for linear interpolation
  static double fs1=44100;
  double fs2=63500;
  double T2f=fs1/fs2;//step through bbd samples, floating point
  UWORD T1=1<<Np; //step through input samples, fixed to 1 sample
  UWORD T2i=(UWORD)(T2f*(1<<Np)+.5);//step through bbd samples, fixed point
  static UWORD T2[3];
  static double T2lfo;
  double width=custom1/100; //percentage of predelay
  static int bbdpointer[3]={0};
  int ch,rread;
  static double lfo=0;
  double v,scl;
  static short comp=(1<<sizeof(short)*8-1)-1;/*maximum short value*/
      UWORD static t;
  static c=-1;

  //time variables
  static UWORD oldBBDevent[3];//initialized below
  static UWORD nextBBDevent[3];//initialized below
  static UWORD nextTevent[3];//initialized below
  double dhw; UWORD dhbw;
  
  if (!init){//first time initialization
     for (ch=0;ch<channels;ch++){
      oldBBDevent[ch]+=*Time-T2[0];
      nextBBDevent[ch]+=*Time;
      nextTevent[ch]=*Time;
      T2[ch]=T2i;
    }
    init=1;
  } 
  else{ //reinizialization
    for (ch=0;ch<channels;ch++){
      oldBBDevent[ch]+=*Time-nextTevent[ch];
      nextBBDevent[ch]+=*Time-nextTevent[ch];
      nextTevent[ch]=*Time;
    }
  }
  int Xp;
  UWORD startTime=*Time;
  UWORD newTime,wi;
  FLOAT *Ystart = Y;
T2lfo=T2f;
  double lfo_width=.3;
  while ( *Time-startTime >>Np < Nx ){
    newTime=*Time+(1<<Np);
    //~ Xp = *Time>>Np;     /* index of X[] to current input sample */
    for (ch=0;ch<channels;ch++){
      if (*Time==nextBBDevent[ch]){
        lfo=lfoGen(preset,ch,*Time);
         T2[ch]=(UWORD)( (T2f + lfo_width*lfo)*(1<<Np) +.5);
        //~ T2lfo= T2f/(width*lfo+1);
        //~ outmin=MIN(outmin,T2lfo);
        //~ outmax=MAX(outmax,T2lfo);
        //~ dhw = MIN(Npc, Npc/T2lfo);  /* Filter sampling period min(Fs,Fs')*/ // Npc > Npc/T2lfo Fs'/Fs
        //~ dhw = Npc;  /* Filter sampling period min(Fs,Fs')*/
        //~ dhbw = dhw*(1<<Na) + 0.5;     /* Fixed-point representation*/
        //~ T2[ch]=(UWORD)( (T2lfo)*(1<<Np) +.5);
        
        //~ scl=MIN(1,1/T2lfo);

        bbdpointer[ch]=(bbdpointer[ch]+1)%bbdsize; //increment in the circular buffer
        z0=X[*Time>>Np];z1=X[1+(*Time>>Np)];w=(*Time&Pmask)/32768.0;
        
        //~ v = FilterUD(Imp, ImpD, Nwing, Interp, Xp, (HWORD)(*Time&Pmask),-1, dhbw, X,comp,0);  /* Perform left-wing inner product */
        //~ v += FilterUD(Imp, ImpD, Nwing, Interp, Xp+1,(HWORD)((((*Time)^Pmask)+1)&Pmask),1, dhbw, X,comp,0);  /* Perform right-wing inner product */
        bbd[ch][bbdpointer[ch]]= z0+(z1-z0)*w;//linear interpolation
        oldBBDevent[ch]=*Time;
        nextBBDevent[ch]=*Time + T2[ch];
      }
      newTime=newTime<nextBBDevent[ch]?newTime:nextBBDevent[ch];
      if (*Time==nextTevent[ch]){
        //~ dhw = MAX(Npc, T2lfo*Npc);  /* Filter sampling period min(Fs,Fs')*/
        //~ dhw = MIN(Npc,T2lfo*Npc);  /* Filter sampling period min(Fs,Fs')*/ // Npc < T2lfo*Npc Fs/Fs'
        //~ dhw = Npc*T2lfo;  /* Filter sampling period min(Fs,Fs')*/
        //~ dhbw = dhw*(1<<Na) + 0.5;     /* Fixed-point representation*/
        
        //~ printf("%.3f read bbd\n",*Time/32768.0);
        rread=(bbdpointer[ch]-bbddel+bbdsize-1)%bbdsize;
        z0=bbd[ch][ rread ];z1=bbd[ch][ (rread+1)%bbdsize];
        w=(*Time-oldBBDevent[ch])/(double)(T2[ch]);
        //~ wi=(UWORD)(w*32767+.5);
        //~ scl=MIN(1,T2lfo);
        
        //~ v=(z0+(z1-z0)*w); //linear interpolation
        //~ v = FilterUD(Imp, ImpD, Nwing, Interp, rread, (HWORD)(wi&Pmask),-1, dhbw, bbd[ch],BBDSIZE,BBDSIZE);  /* Perform left-wing inner product */
        //~ v += FilterUD(Imp, ImpD, Nwing, Interp, (rread+1)%bbdsize,(HWORD)((((wi)^Pmask)+1)&Pmask),1, dhbw, bbd[ch],BBDSIZE,BBDSIZE);  /* Perform right-wing inner product */
                                      c++;
        *Y+=(z0+(z1-z0)*w)/rescale;
        //~ *Y+=bbd[ch][rread];
        //~ *Y=(1/T2lfo -2)*.5;
        //~ *Y=lfo;
        nextTevent[ch]=*Time+(1<<Np); //integer
        if (ch==channels-1) Y++;
      }
      newTime=newTime<nextTevent[ch]?newTime:nextTevent[ch];
    }
    *Time=newTime;/*min(nextBBDevent[],nextTevent[])*/
  }
  reinitLfoGen(*Time);
  t=*Time;
  return (Y - Ystart);        /* Return the number of output samples */
}

static int err_ret(char *s)
{
    fprintf(stderr,"resample: %s \n\n",s); /* Display error message  */
    return -1;
}


static int resampleWithFilter(  /* number of output samples returned */
    double factor,              /* factor = outSampleRate/inSampleRate */
    int infd,                   /* input and output file descriptors */
    int outfd,
    int inCount,                /* number of input samples to convert */
    int outCount,               /* number of output samples to compute */
    int nChans,                 /* number of sound channels (1 or 2) */
    BOOL interpFilt,            /* TRUE means interpolate filter coeffs */
    HWORD Imp[], HWORD ImpD[],
    UHWORD LpScl, UHWORD Nmult, UHWORD Nwing)
{
    static UWORD Time;          /* Current time/pos in input sample */
    UHWORD Xp, Xread;
    //~ int OBUFFSIZE = (int)(((double)IBUFFSIZE)*factor+2.0)*10;
    int OBUFFSIZE=IBUFFSIZE;
    FLOAT X1[IBUFFSIZE];
    FLOAT Y1[IBUFFSIZE+Blgt]={0}; /* I/O buffers */
    FLOAT X2[IBUFFSIZE];/* I/O buffers */
    FLOAT static Y2[IBUFFSIZE+Algt]={0}; /* I/O buffers */
    UHWORD Nout, Nx;
    int i, Ycount, last;
    
    mus_sample_t **obufs = sndlib_allocate_buffers(1, OBUFFSIZE+Algt);
    if (obufs == NULL)
        return err_ret("Can't allocate output buffers");

    /* Account for increased filter gain when using factors less than 1 */
    //~ if (factor < 1)
      //~ LpScl = LpScl*factor + 0.5;

    /* Calc reach of LP filter wing & give some creeping room */
    Xoff = ((Nmult+1)/2.0) * MAX(1.0,1.0/factor) + 20;

    if (IBUFFSIZE < 2*Xoff)      /* Check input buffer size */
      return err_ret("IBUFFSIZE (or factor) is too small");

    Nx = IBUFFSIZE - 2*Xoff;     /* # of samples available to each iteration of the src */
    
    last = 0;                   /* Have not read last input sample yet */
    Ycount = 0;                 /* Current sample and length of output file */
    Xp = Xoff;                  /* Current "now"-sample pointer for input */
    Xread = Xoff;               /* Position in input array to read into */
    for(i=0;i<3;i++)
        Time = (Xoff<<Np);          /* Current-time pointer for converter */
    
    for (i=0; i<Xoff; X1[i++]=0); /* Need Xoff zeros at begining of sample */
    for (i=0; i<Xoff; X2[i++]=0); /* Need Xoff zeros at begining of sample */

    do {
        if (!last)              /* If haven't read last sample yet */
        {
            last = readData(infd, inCount, X1, X2, IBUFFSIZE, 
                            nChans, (int)Xread);
            if (last && (last-Xoff<Nx)) { /* If last sample has been read... */
                Nx = last-Xoff+1; /* ...calc last sample affected by filter */
                if (Nx <= 0)
                  break;
            }
        }
        /* Resample stuff in input buffer */
        //~0: "slow": only the slow sinewave is used as CV for the chorus. Only the output of the first stage is used
        //~ and mixed at a 1:1 rate with the "direct" signal
        //~1: "chor1": only the slow sinewave is used as CV for the chorus, though this time all of the three stages
        //~ are used, and their outputs are mixed together a a 1:1:1 rate
        //~2: "chor2": the CV waveform is the sum of the slow and fast sinewaves, all three stages are used and their
        //~ outputs are mixed together at a 1:1:1 rate
        //~3: "fast": the CV waveform consists of the fast sine only.
        //~ Only the output of the first stage is used and mixed at a 1:1 rate with the "direct" signal(the schematics
        //~ read "vibrato" for this preset, though it is in fact a chorus)
        //~ Preset 5 -default- is the dry channel only
        int channels,direct,preset,channel,dryonly=0,onboard=1;
        double rescale;
        preset=custom3;
        
        switch (preset){
        case 0 :/*slow*/
        case 3 :/*fast*/
            channels=1;
            direct=1;
            rescale=2;
            break;
        case 1:/*chor1*/
        case 2:/*chor2*/
            channels=3;
            direct=0;
            rescale=3;
            break;
        case 4:/*slow wet only*/
            preset=0;
            channels=1;
            direct=0;
            rescale=1;/*should be 1. Debug*/
            break;
        case 5:/*fast wet only*/
            preset=3;
            channels=1;
            direct=0;
            rescale=1;/*should be 1. Debug*/
            break;
        case 6: /*dry only*/
            channels=1;
            direct=1;
            rescale=1;/*should be 1. Debug*/
            dryonly=1;
            break;
        case 7:  /*external lfo wet only*/
            channels=1;
            direct=0;
            rescale=1;
            onboard=0;
            break;
        case 8:/*external lfo mix*/
            channels=1;
            direct=1;
            rescale=2;
            onboard=0;
            break;
        }
        
        memset(&Y1[Blgt], 0, (IBUFFSIZE)*sizeof(Y1[0]));
        //~ for (i=1;i<=channels;channel++){
        
            Nout=SrcUD(X1,Y1+Blgt,factor,&Time,Nx,Nwing,LpScl,Imp,ImpD,interpFilt,X2,preset,channels,rescale,onboard);
        //~ }
        
          if (direct && !dryonly){ /*sums direct signal to the mix*/
            for (i=0;i<IBUFFSIZE;i++){
                Y1[i+Blgt]+=X1[i]/rescale;
            }
        }
        if (direct && dryonly){ /*overwrites the mix with the  direct signal. it's a debug hack*/
            for (i=0;i<IBUFFSIZE;i++){
                Y1[i+Blgt]=X1[i]/rescale;
            }
        }
        

        Xp=Time>>15;
        Time=(Time&Pmask) + (Xoff<<Np);
        //~ if  (Time>>Np!= 17 || Nout!=4062 || Xp!=4079)
            //~ {printf("\nTime is: %u,%05u",Time>>Np,Time&Pmask);printf("\nNout is: %u",Nout );printf("\nXp is : %u\n\n\n",Xp);fflush(stdout);}
        for (i=0; i<IBUFFSIZE-Xp+Xoff; i++) { /* Copy part of input signal */
            X1[i] = X1[i+Xp-Xoff]; /* that must be re-used */
            X2[i] = X2[i+Xp-Xoff];
        }
        if (last) {             /* If near end of sample... */
            last -= Xp;         /* ...keep track were it ends */
            if (!last)          /* Lengthen input by 1 sample if... */
              last++;           /* ...needed to keep flag TRUE */
        }
        Xread = i;              /* Pos in input buff to read new data into */
        Xp = Xoff;
        
        Ycount += Nout;
        if (Ycount>outCount) {
            Nout -= (Ycount-outCount);
            Ycount = outCount;
        }

        if (Nout > OBUFFSIZE) /* Check to see if output buff overflowed */
          return err_ret("Output array overflow");
        
        //~ filtiir(&Y1[Blgt],Nout,&Y2[Algt]);/*applies iir filter*/
        for (i = 0; i < Nout; i++) Y2[i+Algt]=Y1[i+Blgt];/*doesn't apply filter, just copies the samples*/
        
      
        
        for(i=0;i<Algt;i++){
            Y2[i]=Y2[i+Nout];/*copies outputs of the filter that must be reused*/
        }
        for(i=0;i<Blgt;i++){
            Y1[i]=Y1[i+Nout];/*copies inputs of the filter that must be reused*/
        }
        int warn=1;
        for (i = 0; i < Nout; i++) {
            obufs[0][i] = FLOAT_TO_MUS_SAMPLE_TYPE(Y2[i+Algt]);
             if (warn==1 && Y2[i]>=1 || Y2[i]<-1){
                //~ fprintf(stderr, "***** overflow ********");
                warn=0;
                //~ fprintf(stderr," %f %d\n",Y1[i],obufs[0][i] );
            }
        }
        /* NB: errors reported within sndlib */
        mus_file_write(outfd, 0, Nout - 1, 1, obufs);

        printf(".");
        fflush(stdout);

    } while (Ycount<outCount&&Nout!=0); /* Continue until done */

    return(Ycount);             /* Return # of samples in output file */
}

int resample(                   /* number of output samples returned */
    double factor,              /* factor = Sndout/Sndin */
    int    infd,                /* input and output file descriptors */
    int    outfd,
    int inCount,                /* number of input samples to convert */
    int outCount,               /* number of output samples to compute */
    int nChans,                 /* number of sound channels (1 or 2) */
    BOOL interpFilt,            /* TRUE means interpolate filter coeffs */
    int fastMode,               /* 0 = highest quality, slowest speed */
    BOOL largeFilter,           /* TRUE means use 65-tap FIR filter */
    char *filterFile,           /* NULL for internal filter, else filename */
    double variable1,                /* custom variable*/
    double variable2,                /* custom variable*/
    double variable3                /* custom variable*/
)
{   
    custom1=variable1;
    custom2=variable2;
    custom3=variable3;
    UHWORD LpScl;               /* Unity-gain scale factor */
    UHWORD Nwing;               /* Filter table size */
    UHWORD Nmult;               /* Filter length for up-conversions */
    HWORD *Imp=0;               /* Filter coefficients */
    HWORD *ImpD=0;              /* ImpD[n] = Imp[n+1]-Imp[n] */

#ifdef DEBUG
    /* Check for illegal constants */
    if (Np >= 16)
      return err_ret("Error: Np>=16");
    if (Nb+Nhg+NLpScl >= 32)
      return err_ret("Error: Nb+Nhg+NLpScl>=32");
    if (Nh+Nb > 32)
      return err_ret("Error: Nh+Nb>32");
#endif
    
    /* Set defaults */

    if (filterFile != NULL && *filterFile != '\0') {
        if (readFilter(filterFile, &Imp, &ImpD, &LpScl, &Nmult, &Nwing))
          return err_ret("could not find filter file, "
               "or syntax error in contents of filter file");
    } else if (largeFilter) {
        Nmult = LARGE_FILTER_NMULT;
        Imp = LARGE_FILTER_IMP;         /* Impulse response */
        ImpD = LARGE_FILTER_IMPD;       /* Impulse response deltas */
        LpScl = LARGE_FILTER_SCALE;     /* Unity-gain scale factor */
        Nwing = LARGE_FILTER_NWING;     /* Filter table length */
    } else {
        Nmult = SMALL_FILTER_NMULT;
        Imp = SMALL_FILTER_IMP;         /* Impulse response */
        ImpD = SMALL_FILTER_IMPD;       /* Impulse response deltas */
        LpScl = SMALL_FILTER_SCALE;     /* Unity-gain scale factor */
        Nwing = SMALL_FILTER_NWING;     /* Filter table length */
    }
#if DEBUG
    fprintf(stderr,"Attenuating resampler scale factor by 0.95 "
            "to reduce probability of clipping\n");
#endif
    LpScl *= 0.95;
    fprintf(stderr,"di=[" );
    int cane= resampleWithFilter(factor,infd,outfd,inCount,outCount,nChans, 
                              interpFilt, Imp, ImpD, LpScl, Nmult, Nwing);
    fprintf(stderr,"];");
    //~ printf("\n# inmin: %f, inmax: %f\n",inmin,inmax);
    printf("\n# outmin: %f, outmax: %f\n",1/outmin,1/outmax);
    //~ printf("\n# mindiff: %d, maxdiff: %d\n",mindiff,maxdiff);
    //~ printf("\n# %d",MUS_SAMPLE_BITS);
    return cane;
}
void BBDread(UWORD *BBDTime,FLOAT *Y[],UWORD endTime,double Fs2,UWORD *rread,HWORD Imp[],HWORD ImpD[],
                UHWORD Nwing,BOOL Interp,UWORD dhbr,FLOAT bbd[],double Sclr){
    while(*BBDTime<=endTime)
    {   
        UWORD inc;
        FLOAT v;
        *BBDTime+=(1<<Np);
        inc=Fs2 * (1<<Np)*.999808 +.5;/*fac .999808*/
        //~ inc=Fs2 * (1<<Np)+.5;/*fac .999808*/
        //~ inc=next(bbdT,(rread)%(TSIZE<<Np),TSIZE);
        *rread+=inc;
        *rread%=(BBDSIZE<<Np);
        v = FilterUD(Imp, ImpD, Nwing, Interp, (*rread)>>Np, (HWORD)( ((int)(*rread)) &Pmask),-1, dhbr, bbd, BBDSIZE,BBDSIZE);  /* Perform left-wing inner product */
        v += FilterUD(Imp, ImpD, Nwing, Interp, ( (*rread)>>Np)+1,(HWORD)((( ((int)(*rread))^Pmask)+1)&Pmask),1, dhbr, bbd, BBDSIZE,BBDSIZE);  /* Perform right-wing inner product */
        **Y += v*Sclr;   /* deposit output */
        (*Y)++;
                                                outmin= outmin>v ?  v : outmin;
                                                outmax=outmax<v ? v : outmax;
    }
}

static int SrcUDold(FLOAT X[], FLOAT Y[], double factor, UWORD *Time,
                 UHWORD Nx, UHWORD Nwing, UHWORD LpScl,
                 HWORD Imp[], HWORD ImpD[], BOOL Interp, FLOAT X2[], int preset, int channel, double rescale, BOOL onboard)
{   
    FLOAT Xp, *Ystart;
    FLOAT v,out;
    double dhw, dhr;                  /* Step through filter impulse response */
    double dtr,dtw;                  /* Step through input signal */
    UWORD endTime, startTime=*Time;/* When Time reaches EndTime, return to user */
    UWORD dhbw,dhbr, dtbw;             /* Fixed-point versions of Dh,Dt */
    double lfo;
    static double fullTime[3]={0};
    double fullTimeMod=44100/fslow;/*where 1/fslow is the least common multiple of 1/fslow and 1/ffast */
    static int i[3]={0};
    static int diff=0;
    int bbdp=BBDOFF, bbdsize=BBDSIZE, bbddel=BBDDEL;/*debug*/
    static short comp=(1<<sizeof(short)*8-1)-1;
    Ystart = Y;
    double Fs1=44100;
    double modWidth=custom1;/*percentage of preDelay*/
    double Fs2;//=BBDDEL/(preDelay/1000);
    Fs2=100000;
    Fs2/=Fs1;
    factor=Fs2;
    double width=modWidth*factor/100.0;
    double Ts2=1/Fs2;

    dhw = MIN(Npc, factor*Npc);  /* Filter sampling period */
    dhbw = dhw*(1<<Na) + 0.5;     /* Fixed-point representation */
    dhr = MIN(Npc, Npc/factor);  /* Filter sampling period */
    dhbr = dhr*(1<<Na) + 0.5;     /* Fixed-point representation */
    FLOAT Sclw=0.739516430206288;
    FLOAT Sclr=0.414253639218221/rescale;
    UWORD BBDTime0=0;
    UWORD inc;
    WORD offset=*Time-BBDTime0;
    static int wwrite;
    static UWORD rread[3]={0};
    while ( *Time-startTime >>Np < Nx )
    {
        Xp = *Time>>Np;     /* index of X[] to current input sample */
        v = FilterUD(Imp, ImpD, Nwing, Interp, Xp, (HWORD)(*Time&Pmask),-1, dhbw, X,comp,0);  /* Perform left-wing inner product */
        v += FilterUD(Imp, ImpD, Nwing, Interp, Xp+1,(HWORD)((((*Time)^Pmask)+1)&Pmask),1, dhbw, X,comp,0);  /* Perform right-wing inner product */
        wwrite=(i[channel]+++BBDDEL)%BBDSIZE;
        out=v*Sclw;
        out=0.059303154394009*pow(out,3)+0.035460624997318*pow(out,2)+0.910946476824995*out+-0.005710256216322; /*waveshaping function*/
        bbd[channel][wwrite]=out;
                                                    inmin= inmin>out ?  out : inmin;
                                                    inmax=inmax<out ? out : inmax;
        
        if (onboard) lfo=lfoGen(preset,channel,32768*fullTime[channel]);/*onboard lfo*/
        else lfo=X2[*Time>>Np];/*external lfo*/
        Fs2=factor+width*lfo;/*linear modulation*/
        Ts2=1/Fs2;
        //~ fprintf(stderr,"%f\n",Ts2);
        
        dtw=Ts2;
        dtbw = dtw*(1<<Np) + 0.5;
        *Time += dtbw;          /* Move to next sample by time increment */
        fullTime[channel]+=dtw;
        fullTime[channel]=fmod(fullTime[channel],fullTimeMod);
        
        BBDread(&BBDTime0, &Y, *Time-offset, Fs2, &rread[channel], Imp, ImpD, Nwing, Interp, dhbr, bbd[channel],Sclr);
                                                    diff=(wwrite-(rread[channel]>>Np)+bbdsize)%BBDSIZE;
                                                    //~ fprintf(stderr,"%d\n",diff);
                                                    mindiff=mindiff>diff ? diff : mindiff;
                                                    maxdiff=maxdiff<diff ? diff : maxdiff;
    }
    return (Y - Ystart);        /* Return the number of output samples */
}

