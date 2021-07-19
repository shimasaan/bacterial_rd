#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/time.h>
#include "mt19937-64.h"
// #include <random>

//==============SYSTEM SIZE=====================//
#define DeltaT 0.05                             // control
#define TIME_E 100                             // culture duration in the initial phase
#define TIME 500                                // culture duration in the second phase
#define initial_cellnum 10                      // number of cells at the beginning of preculture
#define initial_cellmax 110000                  // number of cells at the end of precultre
#define initial_cellnum2 10                     // number of cells at the beginning of observation
#define shuffle_cellamount 500                   // number of cells at the onset of starvation
#define max_cellamount 300000                   // cells are not added when the population exceeds this value


//===========DISTRIBUTION PARAMETERS============//
#define initiation_vol 0.55                      // 0.9, if completely consistent with the M9 experiment (Wallden et al. 2016) (FREE)
#define sigma_Ithre_coef 0.1                     // std/mean of C period initiation volume (FREE) 0.2

#define E_mean_alp 0.029                            // mean growth rate
#define sigma_septum 0.065                         // std of the septum position
#define replication_sigma 0.05                 // std of the replication threshold in the 1st phase


//===========chromosome replication speed=======//
#define E_mean_vol 4.9
#define E_mean_rep_speed (1.0/(1.3*std::pow(E_mean_alp, -0.84)+42)/E_mean_vol)


//==============STATIC VARIABLES================//
static std::vector<double> l, I, Iinit, Ithre, alp, Xthre;     //cell volume, virtual volume, volume at initiation, threshold of initiator, growth rate, replication threshold
static std::vector<double> t2, t4, t8, t16;     //initiation times of CD period
static std::vector<double> nC;                  //number of oriC: 1, 2, 4, 8 or 16
static double t;                                //external time
static int cell_count;                          //cell number
static std::vector<double> divtime, divsize;    //scatter data of division sizes in the stationary

//tamporal containers for preculture
static double cl[initial_cellnum2], cI[initial_cellnum2], cIinit[initial_cellnum2], cIthre[initial_cellnum2], calp[initial_cellnum2], cXthre[initial_cellnum2];
static double ct2[initial_cellnum2], ct4[initial_cellnum2], ct8[initial_cellnum2], ct16[initial_cellnum2]; 
static double cnC[initial_cellnum2];


//==============FUNCTIONS=======================//

int randomt(int x){ //random number
	struct timeval tv;
	gettimeofday(&tv, NULL);
	init_genrand64(tv.tv_sec + tv.tv_usec);
    return((int)x*genrand64_real2());
}

double randomd(){ //random number
	struct timeval tv;
	gettimeofday(&tv, NULL);
	init_genrand64(tv.tv_sec + tv.tv_usec);
    return(genrand64_real2());
}

double gauss(double mu, double sigma){// (make gaussian random number)
    if(sigma==0){return mu;}

 	struct timeval tv;
 	gettimeofday(&tv, NULL);
 	init_genrand64(tv.tv_sec + tv.tv_usec);

    double z=sqrt( -2.0*log(genrand64_real3()) ) * sin( 2.0*M_PI*genrand64_real3() );
    return mu + sigma*z;
}

void make_initial(){// make random initial condition
    l.reserve(initial_cellmax); // cell volume
    l.resize(initial_cellnum,0);
    I.reserve(initial_cellmax); // virtual total cell volume (amount of the initiator)
    I.resize(initial_cellnum,0);
    Iinit.reserve(initial_cellmax); // amount of the initiator at initiation
    Iinit.resize(initial_cellnum,0);
	Ithre.reserve(initial_cellmax); // replication initiation cell volume
    Ithre.resize(initial_cellnum,0);
    alp.reserve(initial_cellmax); // elongation rate
	alp.resize(initial_cellnum,0);
    Xthre.reserve(initial_cellmax); // replication size
	Xthre.resize(initial_cellnum,0);
    nC.reserve(initial_cellmax); // number of oriC
    nC.resize(initial_cellnum,0);
    t2.reserve(initial_cellmax); // oldest initiation time
	t2.resize(initial_cellnum,-1000);
    t4.reserve(initial_cellmax); // second oldest initiation time
    t4.resize(initial_cellnum,-1000);
    t8.reserve(initial_cellmax); // 3rd oldest initiation time
    t8.resize(initial_cellnum,-1000);
    t16.reserve(initial_cellmax); // 4th oldest initiation time
    t16.resize(initial_cellnum,-1000);

    divtime.reserve(initial_cellmax);
    divsize.reserve(initial_cellmax);

    for(int i=0; i<initial_cellnum; i++){
        l[i] = 0.07+0.06*randomd();
        I[i] = l[i];
        Iinit[i] = 0;
        nC[i] = 1;
        Ithre[i] = gauss(nC[i], sigma_Ithre_coef*nC[i]);
        alp[i] = E_mean_alp;
        Xthre[i] = gauss(1, replication_sigma);
    }

    t = 0;
    cell_count = initial_cellnum; //start from a single cell
}

void cell_division(int label){
    double sep_pos = gauss(1, sigma_septum); //get the septum position
    l.push_back(l[label]*(sep_pos)/2);
    l[label] = l[label]*(2-sep_pos)/2;

    I.push_back(I[label]);
    Iinit.push_back(Iinit[label]);
    nC.push_back(nC[label]/2);
    nC[label] = nC[label]/2;
    Ithre.push_back(Ithre[label]);
    Ithre[label] = Ithre[label];
    alp.push_back(1);
    alp[label] = 1;
    Xthre.push_back(gauss(1, replication_sigma));
    Xthre[label] = gauss(1, replication_sigma);

    if(nC[label]==1){ // if there is no proceeding replication, set tau as infinity
        t2[label] = -1000;
        t2.push_back(-1000);
        t4.push_back(-1000);
        t8.push_back(-1000);
        t16.push_back(-1000);
    }
    else{ // if there IS proceeding replication, reset devision time and initiation times
        t2[label] = t4[label];
        t2.push_back(t4[label]);
        t4[label] = t8[label];
        t4.push_back(t8[label]);
        t8[label] = t16[label];
        t8.push_back(t16[label]);
        t16.push_back(-1000);
    }

    cell_count++;
}

void cell_division_nopush(int label){
    double sep_pos = gauss(1, sigma_septum); //get the septum position
    l[label] = l[label]*(2-sep_pos)/2;

    nC[label] = nC[label]/2;
    Ithre[label] = Ithre[label];
    Xthre[label] = gauss(1, replication_sigma);

    if(nC[label]==1){ // if there is no proceeding replication, set tau as infinity
        t2[label] = -1000;
    }
    else{ // if there IS proceeding replication, reset devision time and initiation times
        t2[label] = t4[label];
        t4[label] = t8[label];
        t8[label] = t16[label];
        t16[label] = -1000;
    }
}

void cell_growth(int label){
    l[label] = l[label]*std::exp(alp[label]*DeltaT);
    I[label] = I[label]*std::exp(alp[label]*DeltaT);
    if(I[label] > Ithre[label]*initiation_vol + Iinit[label]){ // initiation of new chromosome replication
        nC[label] = 2*nC[label];
        I[label] = l[label]; // volume at initiation
        int n = nC[label];
        switch(n){
            case 2:
                t2[label] = 0;
                break;
            case 4:
                t4[label] = 0;
                break;
            case 8:
                t8[label] = 0;
                break;
            case 16:
                t16[label] = 0;
                break;
        }
        Ithre[label] = gauss(nC[label], sigma_Ithre_coef*nC[label]); // reset next initiation volume
        Iinit[label] = l[label];
    }
}

void pre_exponential_phase(){
    while(cell_count < initial_cellmax){
        t = t + DeltaT; //update the external clock
        int cell_count_s = cell_count; //save for the next for loop (avoid changing the conditional formula)

        for(int i=0; i<cell_count_s; i++){
            t2[i] = t2[i]+DeltaT*E_mean_rep_speed*l[i];
            t4[i] = t4[i]+DeltaT*E_mean_rep_speed*l[i];
            t8[i] = t8[i]+DeltaT*E_mean_rep_speed*l[i];
            t16[i] = t16[i]+DeltaT*E_mean_rep_speed*l[i];

            if(t2[i] > Xthre[i]){ //cell division
                cell_division(i);
            }
            else{ //cell growth
                alp[i] = E_mean_alp;
                cell_growth(i);
            }
        }
    }
    
    for(int i=0; i<initial_cellnum2; i++){
        int label = randomt(initial_cellmax);
        cl[i] = l[label];
        while(cl[i] > .52*E_mean_vol){
            label = randomt(initial_cellmax);
            cl[i] = l[label];
        }
        cI[i] = I[label];
        cIinit[i] = Iinit[label];
        cnC[i] = nC[label];
        cIthre[i] = Ithre[label];
        calp[i] = alp[label];
        cXthre[i] = Xthre[label];

        ct2[i] = t2[label];
        ct4[i] = t4[label];
        ct8[i] = t8[label];
        ct16[i] = t16[label];
    }

    l.reserve(max_cellamount); // cell volume
    l.resize(initial_cellnum2,0);
    I.reserve(max_cellamount); // virtual total cell volume (amount of the initiator)
    I.resize(initial_cellnum2,0);
    Iinit.reserve(max_cellamount); // amount of the initiator at initiation
    Iinit.resize(initial_cellnum2,0);
	Ithre.reserve(max_cellamount); // replication initiation cell volume
    Ithre.resize(initial_cellnum2,0);
    alp.reserve(max_cellamount); // elongation rate
	alp.resize(initial_cellnum2,0);
    nC.reserve(max_cellamount); // number of oriC
    nC.resize(initial_cellnum2,0);
    Xthre.reserve(max_cellamount); // replication size
    Xthre.resize(initial_cellnum2,0);
    t2.reserve(max_cellamount); // oldest initiation time
	t2.resize(initial_cellnum2,-1000);
    t4.reserve(max_cellamount); // second oldest initiation time
    t4.resize(initial_cellnum2,-1000);
    t8.reserve(max_cellamount); // 3rd oldest initiation time
    t8.resize(initial_cellnum2,-1000);
    t16.reserve(max_cellamount); // 4th oldest initiation time
    t16.resize(initial_cellnum2,-1000);

    t = 0;
    cell_count = initial_cellnum2; //start from a few cells
    for(int i=0; i<initial_cellnum2; i++){
        l[i] = cl[i];
        I[i] = cI[i];
        Iinit[i] = cIinit[i];
        nC[i] = cnC[i];
        Ithre[i] = cIthre[i];
        alp[i] = calp[i];
        Xthre[i] = cXthre[i];
        t2[i] = ct2[i];
        t4[i] = ct4[i];
        t8[i] = ct8[i];
        t16[i] = ct16[i];
    }

    std::cout << "Preculture end." << std::endl;
}

void exponential_phase(){
    double sep_pos, checker1;
    checker1 = 1.5;

    std::ofstream ofs1("./E_time.csv"); // time points in exponential phase
    std::ofstream ofs2("./E_volume_lists.csv"); // list of cell volume in the last period of Exponential phase

    while(cell_count < shuffle_cellamount){ //increase the population
        int cell_count_s = cell_count; //save for the next for loop (avoid changing the conditional formula)
        for(int i=0; i<cell_count_s; i++){
            t2[i] = t2[i]+DeltaT*E_mean_rep_speed*l[i];
            t4[i] = t4[i]+DeltaT*E_mean_rep_speed*l[i];
            t8[i] = t8[i]+DeltaT*E_mean_rep_speed*l[i];
            t16[i] = t16[i]+DeltaT*E_mean_rep_speed*l[i];

            if(t2[i] > Xthre[i]){ //cell division
                cell_division(i);
            }
            else{ //cell growth
                alp[i] = E_mean_alp;
                cell_growth(i);
            }
        }
    }
    while(t < 15000){ // shuffle the cell cycle across individuals
        t = t + DeltaT; //update the external clock
        for(int i=0; i<cell_count; i++){
            t2[i] = t2[i]+DeltaT*E_mean_rep_speed*l[i];
            t4[i] = t4[i]+DeltaT*E_mean_rep_speed*l[i];
            t8[i] = t8[i]+DeltaT*E_mean_rep_speed*l[i];
            t16[i] = t16[i]+DeltaT*E_mean_rep_speed*l[i];

            if(t2[i] > Xthre[i]){ //cell division
                cell_division_nopush(i);
            }
            else{ //cell growth
                alp[i] = E_mean_alp;
                cell_growth(i);
            }
        }
    }

    t = -TIME_E;
    while(t < 0){
        t = t + DeltaT; //update the external clock
        for(int i=0; i<cell_count; i++){
            t2[i] = t2[i]+DeltaT*E_mean_rep_speed*l[i];
            t4[i] = t4[i]+DeltaT*E_mean_rep_speed*l[i];
            t8[i] = t8[i]+DeltaT*E_mean_rep_speed*l[i];
            t16[i] = t16[i]+DeltaT*E_mean_rep_speed*l[i];

            if(t2[i] > Xthre[i]){ //cell division
                divtime.push_back(t);
                divsize.push_back(l[i]);
                cell_division(i);
            }
            else{ //cell growth
                alp[i] = E_mean_alp;
                cell_growth(i);
            }

            //========= export data =========//
            if(t+TIME_E > checker1){
                ofs1 << t << std::flush;
                if(t < 0){ofs1 << "," << std::flush;}

                for(int i=0; i<cell_count; i++){
                    ofs2 << l[i] << std::flush;
                    if(i != cell_count-1){ofs2 << "," << std::flush;}
                }
                ofs2 << std::endl;

                checker1 += 1.5;
            }
        }
    }

    double sum = 0;
    for(int m=0; m<l.size(); m++){sum += l[m];}
    std::cout << "Mean: " << sum/l.size() << std::endl;
    double max_l = *std::max_element(l.begin(),l.end());
    std::cout << "Max: " << max_l << std::endl;
    std::cout << "Num: " << l.size() << std::endl;
}

void stationary_phase(double kc, double tau1){
    double checker2 = 1.5;
    double meanalp, rep_speed;
    std::ofstream ofs3("./S_time.csv"); // time points in exponential phase
    std::ofstream ofs4("./S_volume_lists.csv"); // cell volume lists in stationary phase

    char fnamed[50];
    sprintf(fnamed,"./cell_div.csv");
    std::ofstream ofsd(fnamed);  

    while(t < TIME){
        t = t + DeltaT; //update the external clock

        meanalp = E_mean_alp*(1-0.84)/(std::exp(0.011*t)-0.84);

        for(int i=0; i<cell_count; i++){
            alp[i] = meanalp;
            rep_speed = E_mean_rep_speed*(kc+1)/(kc*std::exp(t/tau1)+1)*l[i];

            // replication proceeds
            t2[i] = t2[i]+DeltaT*rep_speed;
            t4[i] = t4[i]+DeltaT*rep_speed;
            t8[i] = t8[i]+DeltaT*rep_speed;
            t16[i] = t16[i]+DeltaT*rep_speed;

            if(t2[i]>Xthre[i]){ // cell division
                divtime.push_back(t);
                divsize.push_back(l[i]);            
                cell_division(i);
            }
            else{ //cell growth
                cell_growth(i);
            }
        }
        
        //========= export cell volume =========//
        if(t > checker2){
            ofs3 << t << std::flush;
            if(t < TIME){ofs3 << "," << std::flush;}

            for(int i=0; i<cell_count; i++){
                ofs4 << l[i] << std::flush;
                if(i != cell_count-1){ofs4 << "," << std::flush;}
            }
            ofs4 << "\n" << std::flush;

            checker2 += 1.5;
        }
    }

    for(int m=0; m<divtime.size(); m++){
        ofsd << divtime[m] << std::flush;
        if(m != divtime.size()-1){ofsd << "," << std::flush;}
    }
    ofsd << std::endl;
    for(int m=0; m<divsize.size(); m++){
        ofsd << divsize[m] << std::flush;
        if(m != divsize.size()-1){ofsd << "," << std::flush;}
    }
}

//====================MAIN======================//

int main(int argc, char *argv[]){
	struct timeval s, e;
	gettimeofday(&s, NULL);

    for(int j=0; j<1; j++){
        double kc = 0.01;
        double tau1 = 40;

        make_initial(); //make initial condition
        pre_exponential_phase(); //make uniform cell samples from steady exponential phase

        exponential_phase();
        stationary_phase(kc, tau1);

        // std::cout << j << std::endl;
    }

    return 0;
}