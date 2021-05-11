//
//  Janez Brest, CEC 2020
//


#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>    // std::ostringstream
#include <cstdlib>    // drand48,   system
#include <cmath>
#include <string>
#include <limits>
#include <ctime>
using namespace std;

// void cec20_test_func(double *, double *,int,int,int);   // CEC 2020

void cec21_bias_shift_rot_func(double *, double *,int,int,int);
void cec21_bias_shift_func(double *, double *,int,int,int);
void cec21_bias_rot_func(double *, double *,int,int,int);
void cec21_shift_rot_func(double *, double *,int,int,int);
void cec21_rot_func(double *, double *,int,int,int);
void cec21_shift_func(double *, double *,int,int,int);
void cec21_bias_func(double *, double *,int,int,int);
void cec21_basic_func(double *, double *,int,int,int);
double *OShift,*M,*y,*z,*x_bound;
int ini_flag=0,n_flag,func_flag,*SS;


// CEC 2020 FUNCTIONS
const int nFun=10;      // number of functions
const int nRuns=30;      // number of runs (30 in CEC 2021)
double REZ[8][nFun][nRuns];   // data for Tables II-V   4...D=10,20

const int rawNum=16; // number of rawdata values
double rawData[8][nFun][nRuns][rawNum];
const string DIRNAME="Results";              // directory for files of raw results

const int myPrec=11;    // setprecision(myPrec) for cout
long mySeed;            //
const double eps=1e-12;
const double terminateErrorValue=1e-8;  // --> ZERO

const int maxNP = 1000;   // max. number for population size
const int maxD = 20;       // max. number for problem/function dimension
//const unsigned long maxFES = (unsigned long)1e12;  // max. number of functions evaluations (very high number)

double P[maxNP][maxD];    // population NP x D; NP is population size, i.e., number of individuals

const int printVectorD=10;  // if D is smaller the print solution vectors
long nReset;                // reset counter
long sReset;                // reset counter

const long EMPTY=-1;         // Table2: If the given level of accuracy is not reached, leave the corresponding cell empty.

/* --------- jDE constants ----------- */
const double Finit  = 0.5;  // F   INITIAL FACTOR VALUE
const double CRinit = 0.9;  // CR  INITIAL FACTOR VALUE

double Fl = 0.1;         //
const double Fu = 1.1;   //

double CRl = 0.0;        //
double CRu = 1.0;  //

const double tao1 = 0.1;   // probability to adjust F
const double tao2 = 0.1;   // probability to adjust CR
const double myEqs = 0.25; // for reset populations  CEC2019:0.25

// [0, 1)
double mRandom() {
  //  return rand()/(RAND_MAX + 1.0);
  return drand48();
}

// [0, val) intValue
unsigned int mRandomInt(int upper) {
  return (int)(upper * mRandom());
}

void initialization(double P[maxNP][maxD], const int NP, const int D, const double Xmin, const double Xmax) {
   for (int i = 0; i<NP; i++)
      for (int j = 0; j<D; j++) {
         P[i][j] = Xmin + mRandom()*(Xmax - Xmin);
      }
}

void printIndividual(const double X[maxD], const int D) {
   cout << " [";
   for (int j = 0; j<D; j++) {
      //if(j>2) { cout << "..."; break; }   // PRINT ONLY THREE
      cout << X[j];
      if (j != D-1) cout << ", ";
   }
   cout << "] ";
}

void printPOP(const double P[maxNP][maxD], const int NP, const int D, const double cost[maxNP]) {  // print POP E
   for (int i = 0; i<NP; i++) {
      cout << i << " ";
      if (D <= printVectorD) printIndividual(P[i], D);  // smaller print else skip
      cout << " value= "<< cost[i] << endl;
   }
}

double Dist(const double A[maxD], const double B[maxD], const int D) {
   double dist=0.0;
   for (int j = 0; j<D; j++)
      dist += (A[j]-B[j])* (A[j]-B[j]);
   // without sqrt(dist)
   return dist;
}

int crowding(const double P[maxNP][maxD], const double U[maxD], const int NP, const int D) {
   double dist = Dist(P[0],U,D);
   double min_dist= dist;
   int min_ind = 0;
   for (int i = 1; i<NP; i++) {
      dist = Dist(P[i],U,D);
      if (dist < min_dist) {
         min_dist = dist;
         min_ind = i;
      }
   }
   return min_ind;
}

// count how many individuals have similar fitness function value as the best one
int stEnakih(const double cost[], const int NP, const double cBest) {
   int eqs=0;  // equals
   for (int i = 0; i<NP; i++) {
      if(fabs(cost[i]-cBest) < eps)
         eqs++;
   }
   return eqs;
}

// Are there too many individuals that are equal (very close based on fitness) to the best one
bool prevecEnakih(const double cost[], const int NP, const double cBest) {
   int eqs=0;  // equals
   for (int i = 0; i<NP; i++) {
      if(fabs(cost[i]-cBest) < eps)
         eqs++;
   }
   //if(eqs>myEqs*NP) return true;
   if(eqs>myEqs*NP && eqs > 2) return true;
   else return false;
}

void swap(double &a, double &b) {
   double t=a; a=b; b=t;
}


// print Tables I (rawData tables) CEC 2020
void genTablesRawData(const double rawData[8][nFun][nRuns][rawNum], const int D) {
   for (int trans = 0; trans < 8; trans++) { // CEC2020 dimensions 5,10,15,20
     string trans_name;
     switch (trans) {
       case 0: trans_name="(000)"; break;
       case 1: trans_name="(001)"; break;
       case 2: trans_name="(010)"; break;
       case 3: trans_name="(011)"; break;
       case 4: trans_name="(100)"; break;
       case 5: trans_name="(101)"; break;
       case 6: trans_name="(110)"; break;
       case 7: trans_name="(111)"; break;

     }
      for (int iFunk=0; iFunk<nFun; iFunk++) {
         ostringstream ss;
         ss << DIRNAME<<"/j2020_"<<trans_name<<"_"<<iFunk+1<<"_"<<D<<".txt";  // for example, j2020_1_5.txt, ...
         string s = ss.str();
         ofstream outFile(s.c_str());

         outFile << fixed << setprecision(8);
         // Transpose raws and runs
         for (int k=0; k<rawNum; k++) {  // for each raw
            for (int iRun=0; iRun<nRuns; iRun++) {  // for each RUN
               outFile << rawData[trans][iFunk][iRun][k]<<" ";
            }
            outFile << endl;
         }
         outFile.close();
      } // iFunk
   } // trans
}

// print Table2-Table5 (CEC 2020
void genTableREZ(const double REZ[8][nFun][nRuns], const int D) {
   for (int trans = 0; trans < 8; trans++) { // CEC2021 dimensions 10, 20
      string trans_name;
      switch (trans) {
       case 0: trans_name="(000)"; break;
       case 1: trans_name="(001)"; break;
       case 2: trans_name="(010)"; break;
       case 3: trans_name="(011)"; break;
       case 4: trans_name="(100)"; break;
       case 5: trans_name="(101)"; break;
       case 6: trans_name="(110)"; break;
       case 7: trans_name="(111)"; break;

     }
      ostringstream ss;
      ss <<"Tables"<<"/table_"<<trans_name<<"_"<<D<<".tex";  // table5D.tex, ...
      string s = ss.str();
      ofstream outFile(s.c_str());
      outFile << fixed << setprecision(4);
      outFile << "\\begin{table}[ht]" << endl
              << "\\centering" << endl
              << "\\caption{ Results for "<<D<<"D. }"<<endl
              << "\\label{tab:"<<D<<"D}" << endl
              << "\\begin{tabular}{crrrrr}" << endl
              << "\\hline" << endl
              << "{Func.} & Best & Worst & Median & Mean & Std \\\\" << endl
              << "\\hline" << endl;

      for (int iFunk=0; iFunk<nFun; iFunk++) {
         double dataOne[nRuns];
         for (int iRun=0; iRun<nRuns; iRun++) {  // for each RUN
            dataOne[iRun] = REZ[trans][iFunk][iRun];
         }
         // sort
         for (int i=0; i<nRuns; i++)
            for (int j=0; j<nRuns-1; j++)
               if(dataOne[j]>dataOne[j+1])
                  swap(dataOne[j], dataOne[j+1]);

         // statistic
         double best=dataOne[0];
         double worst=dataOne[nRuns-1];
         double median=dataOne[nRuns/2];
         double mean, std;

         double sum=0.0;
         for (int i=0; i<nRuns; i++) sum += dataOne[i];
         mean=sum/nRuns;

         sum=0;
         for (int i=0; i<nRuns; i++)
            sum += (dataOne[i]-mean)*(dataOne[i]-mean);
         sum *= (1.0/(nRuns-1));
         std = sqrt(sum);
         outFile << iFunk+1 << " & " << best << " & " << worst << " & " << median << " & " <<  mean << " & " <<  std << " \\\\" << endl;
      } // iFunk
      outFile << "\\hline"<<endl;
      outFile << "\\end{tabular}"<<endl;
      outFile << "\\end{table}"<<endl;
   }
}

/**
// print table TIMES (Table IV -- algorithm complexity)
void genTableTime(double T2[4][nRuns]) {
   // T0
   clock_t c_start = clock();
   double x,y;
   for (int i=1; i <= 1000000; i++) {
      x= 0.55 + (double) i;
      //x=x + x; x=x./2; x=x*x; x=sqrt(x); x=log(x); x=exp(x); y=x/x;
      x=x + x; x=x/2; x=x*x; x=sqrt(x); x=log(x); x=exp(x); y=x/x;
   }
   clock_t c_end = clock();
   double T0 = (c_end-c_start) / (double)CLOCKS_PER_SEC;  // time [s]
   cout << "T0="<<T0<<" [s]"<<endl;  // here 0  since compiler optimized code
   // measure with -O0  ::  Time T0: 0.066323
   //   T0=0.066323;

   // T1  200000 evals of F7
   double T1[4];  // times for each dimensions
   double c; //, sum=0;
   double vec[maxD];
   for (int D=5; D<=20; D=D+5) {
      c_start = clock();
      for (int i=1; i <= 200000; i++) {
         for (int j=0; j<D; j++) vec[j] = mRandom();
         cec20_test_func(vec, &c, D, 1, 7);  // evaluation results is returned in c
         //sum+=c; //
      }
      c_end = clock();
      T1[D/5-1] = (c_end-c_start) / (double)CLOCKS_PER_SEC;  // time [s]
   }
   //cout << "   (sum="<<sum<<")"<<endl;

   // T2 ... calc mean of all runs
   double sum;
   double meanT2[4];
   for (int D=5; D<=20; D=D+5) {
      sum = 0;
      for (int iRun=0; iRun<nRuns; iRun++) {
        sum += T2[D/5-1][iRun];
      }
      meanT2[D/5-1] =  sum / nRuns;   // MEAN VALUE of T2
   }
   ofstream outFile ("tableAlgComp.tex");
   outFile << setprecision(4);
   outFile << "\\begin{table}[ht]" << endl
           << "\\centering" << endl
           << "\\caption{ Computation complexity. }"<<endl
           << "\\label{tab:AlgComp}" << endl
           << "\\begin{tabular}{ccccc}" << endl
           << "\\hline" << endl
           << "  & T0 & T1 & {\\it T2} & {\\it (T2 - T1)/T0} ???? \\\\" << endl
           << "\\hline" << endl;

   for (int D=5; D<=20; D=D+5) {
      if(D==20) continue;
      outFile << "$D=$" << D << " & " << T0 << " & " << T1[D/5-1] << " & " << meanT2[D/5-1] << " & " << (meanT2[D/5-1]-T1[D/5-1])/T0 <<" \\\\" << endl;
   }
   outFile << "\\hline" << endl;
   outFile << "\\hline"<<endl;
   outFile << "\\end{tabular}"<<endl;
   outFile << "\\end{table}"<<endl;
   outFile.close();
}
**/


//====================================================================================
int main() {

   stringstream ss;
   ss << "mkdir -p " << DIRNAME;
   string s = ss.str();
   int tmpRez;
   tmpRez=system(s.c_str());  // This should work on LINUX

   double T2[4][nRuns];  // we measure for all runs for each dimenstion

   cout << "******************************** " << endl;
   cout << "*** Algorithm j2020 CEC 2020 *** " << endl;
   cout << "******************************** " << endl;
   int NP, bNP,sNP; // two populations, big and small
   unsigned long maxFES = 0; // depends on dimension

   // cout << "seedInit for random generator (long value): "; cin >> mySeed;
   // mySeed=1;
   mySeed=time(NULL);  //
   cout << "seed: "<<mySeed<<endl;
   srand48(mySeed);    // seeding random generator drand48

   // CEC 2020 FUNCTIONS

   //int D;                // DIMENSION
   double cost[maxNP];                 // vector of population costs (fitnesses)
   int indBest = 0;                    // index of best individual in current population
   double globalBEST[maxD];            // copy of GLOBAL BEST (we have restarts)
   double bestCOST = std::numeric_limits<double>::max(); // cost of globalBEST
   double parF[maxNP], parCR[maxNP];   // each individual has own F and CR
   double U[maxD];                     // trial vector

   // lower and upper bounds (the same for all components) CEC2020
   const double Xmin = -100.0;
   const double Xmax = +100.0;

   unsigned long it;    // iteration counter
   unsigned long age;   // how many times the best cost does not improve
   unsigned long cCopy; // counter for copy

   for (int D = 10; D <= 20; D=D+10) { // CEC2020 dimensions 10, 20

      bNP=D*7;     // cec2019: bNP=1000
      sNP=D;       //D*1;
      NP=bNP+sNP;   // both population size together
      if (NP < 5 || NP > maxNP) { cout << "Error maxNP! Increase maxNP"<<endl; exit(1); }

      switch (D) {                   // CEC 2021
         case 10: maxFES =  200000; break;
         case 20: maxFES = 1000000; break;
      }
     for (int trans = 0; trans < 8; trans++) {  // for each TRANFORMATION
        for (int iFunk = 0; iFunk < nFun; iFunk++) {  // for each FUNCTION
           double optimum = 0;
           if (trans > 3) {
             switch (iFunk+1) {            // CEC 2021
                case 1: optimum =  100; break;
                case 2: optimum = 1100; break;
                case 3: optimum =  700; break;
                case 4: optimum = 1900; break;
                case 5: optimum = 1700; break;
                case 6: optimum = 1600; break;
                case 7: optimum = 2100; break;
                case 8: optimum = 2200; break;
                case 9: optimum = 2400; break;
                case 10: optimum= 2500; break;
             }
           }
           cout << "Function number=F"<<iFunk+1<<endl;
           cout << "Transformation T="<<trans<<endl;
           cout << "Dimension D="<<D<<endl;
           cout << "Search range: Xmin="<<Xmin<<" Xmax="<<Xmax<<endl;
           cout << "Population size NP="<<NP<<" (bNP="<<bNP<<", sNP="<<sNP<<") "<<endl;

           for (int iRun=0; iRun<nRuns; iRun++) {  // for each RUN
              clock_t c_start = clock();
              bestCOST = std::numeric_limits<double>::max(); // cost of globalBEST
              nReset=0;
              sReset=0;
              cCopy=0;
              age=0;
              indBest = 0;

              for (int k=0; k<rawNum; k++)
                 rawData[trans][iFunk][iRun][k] = -1; // -1 means no data

              initialization(P,NP,D, Xmin, Xmax);

              for (int i=0; i<NP; i++) {
                 parF[i] = Finit;     // init
                 parCR[i]= CRinit;    // init
              }

              // evaluate initialized population
              for (int i=0; i<NP; i++) {

                 // cec20_test_func(P[i], &cost[i], D, 1, iFunk+1);  // evaluation results is returned in cost[i]
                 switch (trans) {            // CEC 2021
                    case 0: cec21_basic_func(P[i], &cost[i], D, 1, iFunk+1); break;
                    case 1: cec21_rot_func(P[i], &cost[i], D, 1, iFunk+1); break;
                    case 2: cec21_shift_func(P[i], &cost[i], D, 1, iFunk+1); break;
                    case 3: cec21_shift_rot_func(P[i], &cost[i], D, 1, iFunk+1); break;
                    case 4: cec21_bias_func(P[i], &cost[i], D, 1, iFunk+1); break;
                    case 5: cec21_bias_rot_func(P[i], &cost[i], D, 1, iFunk+1); break;
                    case 6: cec21_bias_shift_func(P[i], &cost[i], D, 1, iFunk+1); break;
                    case 7: cec21_bias_shift_rot_func(P[i], &cost[i], D, 1, iFunk+1); break;
                   // cec21_bias_shift_func(pop[i],  &fitness[i], problem_size, 1, function_number);
                   // cec21_bias_rot_func(pop[i],  &fitness[i], problem_size, 1, function_number);
                   // cec21_shift_rot_func(pop[i],  &fitness[i], problem_size, 1, function_number);
                   // cec21_bias_func(pop[i],  &fitness[i], problem_size, 1, function_number);
                   // cec21_shift_func(pop[i],  &fitness[i], problem_size, 1, function_number);
                   // cec21_rot_func(pop[i],  &fitness[i], problem_size, 1, function_number);
                   // cec21_basic_func(pop[i],  &fitness[i], problem_size, 1, function_number);
                 }
                 cost[i] -= optimum; // CEC 2020
                 if(cost[i] < terminateErrorValue) cost[i] = 0.0;

                 if (cost[i] < cost[indBest] || i==0) {  // i == 0  for output only
                    indBest = i;
                    //cout << "it="<<i<<" value="<<cost[indBest]<<" "; if (D <= printVectorD) printIndividual(P[indBest],D); cout << endl;
                 }

                 // CEC 2020
                 for (int k=0; k<16; k++) {
                    if (i+1 == (int)(pow(D,(k/5.0 - 3.0)) * maxFES)) {
                       rawData[trans][iFunk][iRun][k] = cost[indBest];
                       //cout << "it="<<i<<" k="<<k<<" rawData[D/5-1][iFunk][iRun][k]="<< rawData[D/5-1][iFunk][iRun][k] << endl;
                    }
                 }

              }
              bestCOST = cost[indBest];
              for (int j=0; j<D; j++)
                 globalBEST[j] = P[indBest][j];

              // main loop:  maxFES..maximum number of evaluations
              for (it=NP; it < maxFES; it++) {  // in initialization, NP function evaluations
                 int i = it % (2*bNP);
                 int r1, r2, r3;
                 double F, CR;
                 double c;

                 // reinitialization big population
                 if (i==0 && ( prevecEnakih(cost,bNP,cost[indBest]) || age > maxFES/10 ) ) {
                    nReset++;
                    for (int w = 0; w<bNP; w++) {
                       for (int j = 0; j<D; j++) {
                          P[w][j] = Xmin + mRandom()*(Xmax - Xmin);
                       }
                       parF[w] = Finit;     // init
                       parCR[w]= CRinit;    // init
                       cost[w]=std::numeric_limits<double>::max();
                    }
                    age=0;
                    indBest=bNP; for (int w = bNP+1; w<NP; w++) { if(cost[w] < cost[indBest]) indBest=w; }
                 }

                 // reinitialization small pop
                 if (i==bNP && indBest>=bNP && prevecEnakih(&cost[bNP],sNP,cost[indBest])) {
                    sReset++;
                    for (int w = bNP; w<NP; w++) {
                       if(indBest==w) continue;
                       for (int j = 0; j<D; j++) {
                          P[w][j] = Xmin + mRandom()*(Xmax - Xmin);
                       }
                       parF[w] = Finit;     // init
                       parCR[w]= CRinit;    // init
                       cost[w]=std::numeric_limits<double>::max();
                    }
                 }

                 if (i==bNP && indBest < bNP) {   // copy best solution from the big pop in the small pop
                    cCopy++;
                    cost[bNP]=cost[indBest];
                    for (int j = 0; j<D; j++) {
                       P[bNP][j]= P[indBest][j];
                    }
                    indBest=bNP;
                 }

                 bool bigNP;
                 if (i < bNP) {
                    bigNP=true;

                    // Parameters for big pop
                    Fl=0.01;
                    CRl=0.0;
                    CRu=1.0;

                    int mig=0;
                    if (it < maxFES/3)
                       mig = 1;
                    else if (it < 2*maxFES/3)
                       mig = 2;
                    else
                       mig = 3;

                    do {
                       r1 = mRandomInt(bNP+1);
                    } while (r1 == i && r1 == indBest );   // r1 also should differ from indBest
                    do {
                          r2 = mRandomInt(bNP+mig);        // HERE: index bNP is the first element in small pop
                    } while (r2 == i || r2 == r1);
                    do {
                       r3 = mRandomInt(bNP+mig);           // HERE: index bNP is the first element in small pop
                    } while (r3 == i || r3 == r2 || r3 == r1);
                 }
                 else {
                    bigNP=false;
                    // Parameters for small pop
                    Fl=0.17;
                    CRl=0.1;
                    CRu=0.7;

                    i = (i-bNP) % sNP;
                    do {
                       r1 = mRandomInt(sNP);
                    } while (r1 == i );
                    do {
                       r2 = mRandomInt(sNP);
                    } while (r2 == i || r2 == r1);
                    do {
                       r3 = mRandomInt(sNP);
                    } while (r3 == i || r3 == r2 || r3 == r1);
                    r1 += bNP;
                    r2 += bNP;
                    r3 += bNP;
                    i += bNP;
                 }

                 int jrand = mRandomInt(D);

                 // SELF-ADAPTATION OF CONTROL PARAMETERS  jDE
                 if (mRandom()<tao1) {                      // F
                    F = Fl + mRandom() * Fu;
                 }
                 else {
                    F = parF[i];
                 }
                 if (mRandom()<tao2) {                      // CR
                    CR = CRl + mRandom() * CRu;
                 }
                 else {
                    CR = parCR[i];
                 }

                 bool border;                        // CEC 2020
                 if(mRandom() < 0.0) border=true;
                 else border=false;

                 for(int j=0; j<D; j++) {    // mutation and crossover
                    if (mRandom() < CR || j == jrand) {
                       U[j] = P[r1][j] + F*(P[r2][j] - P[r3][j]);  // DE/rand/1/bin (jDEbin)

                       if (border) {
                          // on border
                          if(U[j] < Xmin) { U[j] = Xmin; }
                          if(U[j] > Xmax) { U[j] = Xmax; }
                       } else {
                          // border check  & repair
                          while(U[j] < Xmin) { U[j] += (Xmax-Xmin); }
                          while(U[j] > Xmax) { U[j] -= (Xmax-Xmin); }
  			if(U[j] < Xmin || U[j] > Xmax) { cerr << "ERROR BORDER" << endl; exit(1); }
                       }
                    }
                    else {
                       U[j] = P[i][j];
                    }
                 }

                 // cec20_test_func(U, &c, D, 1, iFunk+1);  // c .. evaluation result
                 switch (trans) {            // CEC 2021
                    case 0: cec21_basic_func(U, &c, D, 1, iFunk+1); break;
                    case 1: cec21_rot_func(U, &c, D, 1, iFunk+1); break;
                    case 2: cec21_shift_func(U, &c, D, 1, iFunk+1); break;
                    case 3: cec21_shift_rot_func(U, &c, D, 1, iFunk+1); break;
                    case 4: cec21_bias_func(U, &c, D, 1, iFunk+1); break;
                    case 5: cec21_bias_rot_func(U, &c, D, 1, iFunk+1); break;
                    case 6: cec21_bias_shift_func(U, &c, D, 1, iFunk+1); break;
                    case 7: cec21_bias_shift_rot_func(U, &c, D, 1, iFunk+1); break;
                 }
                 c -= optimum;   // CEC 2020
                 if(c < terminateErrorValue) c = 0.0;

                 if(i<bNP) age++;

                 if(i<bNP) i=crowding(P, U, bNP, D);


                 // selection
                 if (c < cost[indBest]) {   // best  (min)
                    age=0;  // reset counter
                    // globalBEST
                    if (c < bestCOST) {
                       bestCOST = c;
                       for (int j=0; j<D; j++)
                          globalBEST[j] = U[j];
                    }

                    cost[i] = c;
                    for (int j=0; j<D; j++) P[i][j] = U[j];
                    parF[i]=F;
                    parCR[i]=CR;
                    indBest = i;

                    if(cost[indBest] < terminateErrorValue) {   // ZERO -- FINISH?
                      // fill uncomplete rawData entries
                      for (int k=0; k<16; k++) {
                         if(rawData[trans][iFunk][iRun][k] == -1) rawData[trans][iFunk][iRun][k] = 0;  // if we have solve problem before maxFES, we fill rawData with zeros
                       }
                       break;
                    }
                 }
                 else if (c <= cost[i]) {               // MIN
                    cost[i] = c;
                    for (int j=0; j<D; j++) P[i][j] = U[j];
                    parF[i]=F;
                    parCR[i]=CR;
                 }

                 // CEC 2020
                 for (int k=0; k<16; k++) {
                    if (it+1 == (unsigned int)(pow(D,(k/5.0 - 3.0)) * maxFES)) {
                       rawData[trans][iFunk][iRun][k] = bestCOST;
                    }
                 }

              } // it  (FEs counter)
              cout <<"Fun=F"<<iFunk+1<<" run="<<iRun+1<<" nReset="<<nReset<<" sReset="<<sReset<<" cCopy="<<cCopy<<" FES="<<it<< " NP="<<NP<<" (bNP="<<bNP<<", sNP="<<sNP<<") D="<<D<< "  bestValue="<<setprecision(myPrec+1+10)<<bestCOST<<" "; /*  printIndividual(globalBEST,D); */  cout <<endl;
              REZ[trans][iFunk][iRun] = bestCOST;

              clock_t c_end = clock();
              T2[D/5-1][iRun] = (c_end-c_start) / (double)CLOCKS_PER_SEC;  // time [s]
              double x = 200000.0*T2[D/5-1][iRun] / maxFES;     // expression for 200000 FES
              T2[D/5-1][iRun] = x;
           } // iRun
        } // iFunk
   }
   genTableREZ(REZ, D);
   genTablesRawData(rawData, D);
   // genTableTime(T2);
} // D

   cout << "Done." << endl << endl;
   return 0;
}
