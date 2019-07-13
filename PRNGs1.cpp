//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%% Mir Mohammad Ebrahimi  %%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// You can compile this code using the following command in shell:
// > g++ -o PRNGs1.out PRNGs1.cpp -Ofast
// , then execute it in shell by
// > ./PRNGs.out

/* These code is based on the random number generator function

   double ran2(long *idum);
   
which is implemented in the Numerical recipes in C, chapter 7 (ran2)

``Long period (> 2 × 10^{18}) random number generator of L. Ecuyer with Bays-Durham shuffle
and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
the endpoint values). 

***!!! Call with idum a negative integer to initialize; !!!*** thereafter, do not alter
idum between successive deviates in a sequence. RNMX should approximate the largest floating
value that is less than 1."

Visit www.nr.com for the licence.*/

#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>

using namespace std;

//------------------------------------------------------------------------------

/* following routine is based on the random number generator

   double ran2(long *idum);
   
which is implemented in the Numerical recipes in C, chapter 7 (ran2)

``Long period (> 2 × 10^{18}) random number generator of L. Ecuyer with Bays-Durham shuffle
and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
the endpoint values). 

***!!! Call with idum a negative integer to initialize; !!!*** thereafter, do not alter
idum between successive deviates in a sequence. RNMX should approximate the largest floating
value that is less than 1."

Visit www.nr.com for the licence.*/

// This is a internal, 32 bit random number generator with uniform Distribution in range [0..1)
/* note #undef's at end of file */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum) {
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    double temp;

    if (*idum <= 0) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

//------------------------------------------------------------------------------

// My interface to the Numerical recipes PRNG 

// Random number generator seed
long iseed = -36;

// You can initializing the random number generator seed by using the following routine.
void Randomize() {
  iseed = -time(NULL);  
}

// Return a random number with uniform distribution in the range of [0..1), where 1 is excluded.
inline double Random() { return ran2(&iseed); }

// Return an integer random number with uniform distribution between 0 to N-1 (both boundary are included).
inline int Random(int N) { return int(ran2(&iseed)*N); }

//------------------------------------------------------------------------------

// Main routine
int main (int argc, char *argv[]) {
    cout << "Hello" << endl;
    
    cout << "time(NULL): " << time(NULL) << endl;

    cout << "iseed: " << iseed << endl;
    
    Randomize();
    
    cout << "iseed: " << iseed << endl;
    
    cout << "random sequence from ran2()\n";
    for (int i=0; i<10; i++)
        cout << Random() << endl;
    
    // Following section of code finds the distribution of ran2()
    // Init
    int N_bin = 15;                         // Number of bins
    int p[N_bin];                           // bins
    for (int i=0; i<N_bin; i++)
        p[i] = 0;                           // clear bins
    
    // Run
    for (int i=0; i<1000000; i++) {
        int i_bin = (int) (Random()*N_bin);
        ++p[i_bin];
    }
    
    // Done
    cout << "p[i]\n";
    for (int i=0; i<N_bin; i++)
        cout << p[i] << endl;

    int sum = 0;
    for (int i=0; i<N_bin; i++)
        sum += p[i];
    cout << "sum: " << sum << endl;

    cout << "normalized\n";
    double r = 1.*N_bin/sum;
    //r = 1./sum;
    for (int i=0; i<N_bin; i++)
        cout << p[i] * r << endl;
    // End of section

    // Following section of code plots the result by the gnuplot (an opensource software in Linux)
    ofstream prng("prng.txt", ios::out | ios::trunc);
    prng << "set ylabel \"y\"\n"
         << "set xlabel \"x\"\n" 
         << "set autoscale\n"
         << "set yrange [0:1.2]\n"
         << "set key off" << endl;

    prng << "plot '-' using 1:2 title \"My Plot\" with lp" << endl;
    
    for(int i=0; i<N_bin; i++)
        prng << (i+.5) / N_bin << '\t' << p[i] * r << endl;
    
    prng << "EOF\n"
         << "pause -1 \"Hit any key to continue\"" << endl;

    prng.close();

    system("gnuplot prng.txt");   // This C++ function executes the gnuplot in main OS
    
    cout << "Finish" << endl;

    return 0;
}
