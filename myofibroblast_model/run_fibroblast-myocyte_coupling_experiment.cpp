/*
Atrial myocyte-Myofibroblast coupling experiment
- Myofibroblast: MacCannell model based on Ashihara et al. (2012)
- Atrial cell: Courtemanche model
*/

#include <stdio.h>
#include "CRN_MacCannell_Fib.h"
#include <math.h>

#define MAX_TIME 20000
#define DT 0.05
#define PRINT_INTERVAL 100
#define REST_INTERVAL 10000.0
#define NUM_STIM 20

double time;

int Nf;
double stim_duration;   // millisecond (in membrane)
double stim_end;   // millisecond (in membrane)
double stim_period;   // millisecond (in membrane)
double stim_start;   // millisecond (in membrane)
double print_start, print_end; // ms
bool rest;

int isStim(double time)
{
    if(!rest) return 0;
    
	if ((time >= stim_start) && (time <= stim_end) && (time-stim_start-floor((time-stim_start)/stim_period)*stim_period <= stim_duration))
		return 1;
	else
		return 0;
}

int main()
{
	int i = PRINT_INTERVAL;
    
	stim_duration = 3.0;   // millisecond (in membrane)
	stim_period = 1000.0;   // millisecond (in membrane)
	stim_start = 50.0;   // millisecond (in membrane)
    stim_end = stim_start + stim_period*(NUM_STIM-1) + 50.0;   // millisecond (in membrane)
    print_start = stim_start + stim_period*(NUM_STIM-1) - 5.0;  // ms
    print_end = print_start + 400.0;
	time = 0.0;
    rest = false;
    
	CRN_Fib_init();
    G_gap_ = 8.0; // nS
    Nf = 2;
    
    FILE *dataout = fopen("vmrecord[G8.0][N2].txt","wt");
    
	while(!rest || time <= MAX_TIME)
	{
		CRN_Fib_compute(DT, isStim(time), Nf);
        
		if (i == PRINT_INTERVAL)
		{
            // No stimulation for REST_INTERVAL after start.
            if(!rest && time >= REST_INTERVAL)
            {
                rest = true;
                time = 0.0;
            }
            
            // Print time, Vm of myocyte, Vm of fibroblast.
            if(rest && time >= print_start && time < print_end)
                fprintf(dataout, "%lf %lf %lf\n", time-print_start, Y[14], YY[0]);
            
            /*
             fprintf(dataout, "%lf", time);
             fprintf(dataout, " %lf", Y[14]);
             for(int j=0; j<_NB_OF_STATE_VARIABLES_Fibroblast_; j++)
             fprintf(dataout, " %lf", YY[j]);
             fprintf(dataout, "\n");
             */
            
			i -= PRINT_INTERVAL;
		}
		
		i++;
		time += DT;
	}
    
	fclose(dataout);
    
	return 0;
}