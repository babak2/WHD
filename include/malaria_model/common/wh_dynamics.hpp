#ifndef WH_DYNAMICS_HPP
#define WH_DYNAMICS_HPP

//#pragma once

#include <iostream>
#include <random>

#include "Eigen/Dense"
#include "Eigen-unsupported/Eigen/MatrixFunctions" 

#include "human/human.h"

#include "util/CSVWriter.h"
#include "util/util.h"

#include "rapidjson/document.h"

using namespace std;
using namespace human;

namespace common 
{

   std::random_device rd;

    // Engine
   //std::mt19937 r_engin(rd());
   std::mt19937 r_engin(1);

   // Gaussian/ Normal Distribtuion
   std::normal_distribution<double> randn{0.0, 1.0}; //normally distributed having zero mean and variance one.
   // Uniform Distribtuion
   std::uniform_real_distribution<double> rand{0.0, 1.0}; //uniformly distributed on the interval (0, 1)
   //used "normrnd" "randn" (identified by rdn)

   //P0[0] = 4*pow(10,4); //nb of initial merozoit?
   //std::normal_distribution<> normrnd(P0[0], 20000.0); 
   //std::normal_distribution<> normrnd(40000.0, 20000.0); 
                  
   class WH_Dynamics //Within Host Dynamics
   { 
      private:
        
         //Note that for easier potential debugging for now, the name of parameters here
         //are kept the same as their Matlab original code counterparts. Many of these
         //can be ideally changed to more meaningful names later. 

         const int timepoints=0;
         const int level = 0;
         //int level = 1;
         // Initialize system
         int asteps=timepoints;
         
         //unsigned int II50 = 2*pow(10,9);  //pyrogenic threshold, 2B
         const unsigned int II50 = 0;  //pyrogenic threshold; 2B, = 2*10^9 = 2*pow(10,9); min_max=[0, ~ 2*10^12[

         const double pmf = 0;  //replication rate: parasites multiply pmf times every 48 hours; =11.0; min_max[0, Z+[
         const double gamma=0; //=0.75;  min_max[0, ~ 30[
         const double maxcyt=0; //=0.225;  min_max[0, ~ 10[
         const double kfever=0; //fever// =1.175;  min_max[0, ~ 80+[
         const double kfever_multiplier = 0;
         const double memory=0; // host’s immunological memory //=0.35; min_max: value change seems to have no effect!
         const double decay=0; //=150; min_max: value change seems to have no effect!
         const double switchrate=0; //=0.03;  min_max[0, R+[
         const double immunestack=0; //=2.0; min_max: value change seems to have no effect!
         const double switcheffect = 0; //=22.0;  min_max[0, ~ 50+[
         const double maxinit=0; //=5e4 = 50000;  min_max[0, ~ 5e4]

         const double etac=0; //=5.0; min_max[0, 100+[
         const double etah=0; // Values of Ycapacity,i above 0.4 correspond to hyperimmunity, and Ycapacity,i then grows towards 1 with a time constant of ?hyperimmunity?=?3 days, regardless of antigen concentrations; =3.0;; min_max[0, ~1000+[
         const unsigned long long int x50ab=0; //=5e10 = 50B = 50000000000; min_max[0, 5e15]
         const double x50ab_multiplier=0;

         const double kab=0; //minimum growth rate at low antigen concentration; =0.115;  min_max[0, ~ 1]
         const double kab_multiplier=0;
         const double abprod=0; //=1.0; min_max[0, Z+[
         const double abdecay=0; //=20.0  //min_max: value change seems to have no effect!

         const double killcyt=0;  //innate immune killing maximal rate; =3.15; min_max[0, 100+[
         const double killab=0;   //adaptive immunity maximal kill rate; =2.15; min_max[0, R+[
         const double killab_multiplier=0;

         const double kmsp=0; // rate of development of merozoite-blocking immunity// Per cycle growth rate of MSP antibody response
                              //unknown, and the effects of its variation are explored; =0.0853; min_max[0, ~ 10]
         const double Cmer=0; //Cmerozoite: Maximum blocking of merozoite invasion by MSP antibodies ?!// =0.358;  min_max[0, R+[

         const double kg_divisor=0; //=1/9.0 = 0.1111111;  min_max[0, ~ 5]
         const double kmg_divisor=0; //=1/10.0 = 0.1;  min_max[0, ~ 10+]
         
         //~~~ Gametocytes develop into mature gametocytes after a mean time of kg days which then live for a mean of kmg days
         double kg=0; //=1/9.0 = 0.1111111;  min_max[0, ~ 5]  //is calculated once by constructor  after which will not vary
         double kmg=0; //=1/10.0 = 0.1;  min_max[0, ~ 10+]  //is calculated once by constructor after which will not vary

         const double kgam=0; // proportion that turns gametocyte - sexual commitment;// =0.01;  min_max[0.01, 0.95]
         const double kgametocyte = 0; // sucess rate per gametocyte in a mosquito; =0.0025;//Gametocyte production rate?!// value change only affect "infect": min_max[0, R+[

         //unsigned long int pt = 4*10^9;  //pyrogenic threshold, 4B
         const unsigned long int pt = 0;  // pyrogenic threshold, 4B; =4*pow(10,9)= 4000000000; min_max[0, ~ 4*10^12]

         const double C50g = 0; //Cytokine level Yinnate at which gametocyte inactivation is 0.5?! //= 0.005 ; min_max: value change seems to have no effect, though C50g > 0
         //unsigned long int gam50 = 6*pow(10,6);  //50% gametocyte sucess at infecting mosquitoes, 6M
         const double gam50 = 0;  // 50% gametocyte sucess at infecting mosquitoes, 6M; =6*pow(10,6)= gam50= 6e+06 = 6000000;
                                  //min_max:[0, ~ 6*10^10]: value change seems to only effect "infect"

         //pk
         const double dosingpk = 0;   //=4000.0; dos in g PK; min_max[0, Z+[
         const double dosingpip = 0;  //=960.0; dos in g Pip; min_max[0, ~ 5000]

         // PD parameters
         const double killratedha=0; // max kill rate DHA; =3.25; min_max[0, ~ 30]
         const double killratepip=0; // Kill rate Piperaquine; =1.75; min_max[0, ~ 30]
         const double ce50dha=0;   //inflexion point DHA; =300.0;  min_max[0, Z+[
         const double ce50pip=0;   //inflexion point Pip; =450.0;  min_max[0, Z+[
         const double hdha=0;    //hill coefficent DHA; =3.0; min_max[0, ~ 30]
         const double hpip=0;   //hill coefficient Piperaquine; =3.0;  min_max[0, ~ 30]

         const int npk=0;//=9    //must be currently fixed (=9) as transition matrices dimensions are currently fixed (for now can be removed from config.json file)
         const int npip=0; //=6  //must be currently fixed (=6) as transition matrices dimensions are currently fixed (for now can be removed from config.json file)

         const double P0_initial_val = 0; //nb of initial parasitaemia (merozoit)
         bool with_drug;
         //bool with_drug = true;
         //bool with_drug = false;

         //--------transition matrices:

         //Jpk:
         const double CL = 0; //=3.2500
         const double V2 = 0; //=5.3750

         const double THETA_3_divident = 0; //0.982, THETA(3)
         const double THETA_3_divisor = 0; //24
         double MT = 0; //~  0.04091666666666666369, THETA(3)/ 24  //is calculated once by constructor after which will not vary

         const double NN = 0; //=7.0

         //Jpip:
         const double ACL2 = 0; //0.75
         const double M_WE2 = 0; //54
         const double AGE2 = 0;  //25
         const double WT2 = 0;  //54
         const double THETA2_1 = 0;  //1329.6 , THETA2[1]*24
         const double EM502 = 0;  //0.575
         const double HILL2 = 0;  //5.51
         double CL2 = 0; //~ 1329.599998  //is calculated once by constructor after which will not vary

         const double Q12 = 0; //=7440.0
         const double Q22 = 0; //=2520.0
         const double V22 = 0; //=69840.0
         const double V32 = 0; //=117840.0
         const double V42 = 0; //=741600.0

         const double THETA2_7_divident = 0; // 2.11, THETA2(7)
         const double THETA2_7_divisor = 0; //24
         double MT2 = 0; // ~ 0.0879  //is calculated by constructor once after which will not vary

         const double NN2 = 0; //=2.0

         //----

         typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_dyn_t;
         typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::DontAlign> vector_dyn_t; //vertical vector
 


      public:

         WH_Dynamics(const int timepoints_v, const int level_v, const unsigned int II50_v, const double pmf_v, const double gamma_v,
                     const double maxcyt_v, const double kfever_v, const double kfever_multiplier_v, const double memory_v, const double decay_v, const double switchrate_v,
                     const double immunestack_v, const double switcheffect_v, const double maxinit_v, const double etac_v, const double etah_v, 
                     const unsigned long long int x50ab_v, const double x50ab_multiplier_v, const double kab_v, const double kab_multiplier_v, const double abprod_v, const double abdecay_v, const double killcyt_v, 
                     const double killab_v, const double killab_multiplier_v, const double kmsp_v, const double Cmer_v, const double kg_divisor_v, const double kmg_divisor_v, const double kgam_v,
                     const double kgametocyte_v, const unsigned long int pt_v, const double C50g_v, const double gam50_v, const double dosingpk_v,
                     const double dosingpip_v, const double killratedha_v, const double killratepip_v, const double ce50dha_v, const double ce50pip_v,
                     const double hdha_v, const double hpip_v, const int npk_v, const int npip_v,  const double P0_initial_val_v, bool with_drug_v, 
                     const double CL_v, const double V2_v, const double THETA_3_divident_v, const double THETA_3_divisor_v, const double NN_v, 
                     const double ACL2_v, const double M_WE2_v, const double AGE2_v, const double WT2_v, const double THETA2_1_v, const double EM502_v,  const double HILL2_v,
                     const double Q12_v, const double Q22_v, const double V22_v, const double V32_v, const double V42_v,
                     const double THETA2_7_divident_v, const double THETA2_7_divisor_v, const double NN2_v);

         ~WH_Dynamics();

         void set_with_drug(bool status);
         bool get_with_drug();

         void export_to_csv_total_MDA_vals(const std::string output_file_prefix) const;

         void step(int i, int agent_id, double& tt, int& mda1, int& mda2, int& mda3, 
                          double &kfeverj, double &x50abj, double &kabj, double &killabj,
                          matrix_dyn_t &vpk, matrix_dyn_t &vvpk, matrix_dyn_t &vpip, matrix_dyn_t &vvpip,
                          vector<double>& P0, vector<double>& AI, vector<double> &infect, vector<double> &MG, 
                          vector<double>& G, vector<double> &MSP, vector<double> &cyt, vector<double> &AB,
                          vector<int> &recieved_drug_timestamps, vector<int> &scheduled_drug_timestamps);

         struct miss_drug_record {
            int t;
            int id;
            double p0;  //P0 value
            int rand_caused; //is it caused by the rand condtion, 1 or 0
         };

         std::vector<miss_drug_record> miss_drug_record_vec;

         unsigned long int fever_threshold;


   };

   //----------------------           

   WH_Dynamics::WH_Dynamics(const int timepoints_v, const int level_v, const unsigned int II50_v, const double pmf_v, const double gamma_v,
                            const double maxcyt_v, const double kfever_v, const double kfever_multiplier_v, const double memory_v, const double decay_v, const double switchrate_v,
                            const double immunestack_v, const double switcheffect_v, const double maxinit_v, const double etac_v, const double etah_v, 
                            const unsigned long long int x50ab_v, const double x50ab_multiplier_v, const double kab_v, const double kab_multiplier_v, const double abprod_v, const double abdecay_v, const double killcyt_v, 
                            const double killab_v, const double killab_multiplier_v, const double kmsp_v, const double Cmer_v, const double kg_divisor_v, const double kmg_divisor_v, const double kgam_v,
                            const double kgametocyte_v, const unsigned long int pt_v, const double C50g_v, const double gam50_v, const double dosingpk_v,
                            const double dosingpip_v, const double killratedha_v, const double killratepip_v, const double ce50dha_v, const double ce50pip_v,
                            const double hdha_v, const double hpip_v, const int npk_v, const int npip_v, const double P0_initial_val_v, bool with_drug_v, 
                            const double CL_v,  const double V2_v, const double THETA_3_divident_v,const double THETA_3_divisor_v, const double NN_v, 
                            const double ACL2_v, const double M_WE2_v, const double AGE2_v, const double WT2_v, const double THETA2_1_v, const double EM502_v,  const double HILL2_v,
                            const double Q12_v, const double Q22_v, const double V22_v, const double V32_v, const double V42_v,
                            const double THETA2_7_divident_v, const double THETA2_7_divisor_v, const double NN2_v):

                            timepoints(timepoints_v), level(level_v), II50(II50_v), pmf(pmf_v), gamma(gamma_v),
                            maxcyt(maxcyt_v), kfever(kfever_v), kfever_multiplier(kfever_multiplier_v), memory(memory_v), decay(decay_v), switchrate(switchrate_v),
                            immunestack(immunestack_v), switcheffect(switcheffect_v), maxinit(maxinit_v), etac(etac_v), etah(etah_v), 
                            x50ab(x50ab_v), x50ab_multiplier(x50ab_multiplier_v), kab(kab_v), kab_multiplier(kab_multiplier_v), abprod(abprod_v), abdecay(abdecay_v), killcyt(killcyt_v), 
                            killab(killab_v), killab_multiplier(killab_multiplier_v), kmsp(kmsp_v), Cmer(Cmer_v), kg_divisor(kg_divisor_v), kmg_divisor(kmg_divisor_v), kgam(kgam_v), 
                            kgametocyte(kgametocyte_v), pt(pt_v), C50g(C50g_v), gam50(gam50_v), dosingpk(dosingpk_v),
                            dosingpip(dosingpip_v), killratedha(killratedha_v), killratepip(killratepip_v), ce50dha(ce50dha_v), ce50pip(ce50pip_v),
                            hdha(hdha_v), hpip(hpip_v), npk(npk_v), npip(npip_v), P0_initial_val(P0_initial_val_v), with_drug(with_drug_v), 
                            CL(CL_v), V2(V2_v), THETA_3_divident(THETA_3_divident_v), THETA_3_divisor(THETA_3_divisor_v), NN(NN_v),
                            ACL2(ACL2_v), M_WE2(M_WE2_v), AGE2(AGE2_v), WT2(WT2_v), THETA2_1(THETA2_1_v), EM502(EM502_v), HILL2(HILL2_v),
                            Q12(Q12_v), Q22(Q22_v), V22(V22_v), V32(V32_v), V42(V42_v),
                            THETA2_7_divident(THETA2_7_divident_v), THETA2_7_divisor(THETA2_7_divisor_v), NN2(NN2_v)

   {

      //resize_matrices(timepoints_v, level_v);
      asteps=timepoints;

      std::cout <<"maxinit= " << maxinit << "\n";
      std::cout <<"x50ab= " << x50ab << "\n";
      std::cout <<"II50= " << II50 << "\n";
      std::cout <<"pt= " << pt << "\n";
      std::cout <<"gam50= " << gam50 << "\n";
      std::cout <<"--++--\n"; 

      kg = 1/ kg_divisor;
      kmg = 1/kmg_divisor;
      double TVMT = THETA_3_divident/ THETA_3_divisor;
      MT = TVMT;

      //MF2 = ((AGE2)^HILL2)/(((EM502)^HILL2)+((AGE2)^HILL2));
      double MF2 = pow(AGE2,HILL2)/(pow(EM502,HILL2)+pow(AGE2,HILL2));
      //std::cout << "MF2= " <<MF2<< " \n";

      //TVCL2 = THETA2(1)*MF2*(WT2/M_WE2)^ACL2;
      double TVCL2 = THETA2_1*MF2*pow((WT2/M_WE2),ACL2);
      //std::cout << "TVCL2= " <<TVCL2<< " \n";
      CL2 = TVCL2;

      //double THETA2_7 = 2.11 / 24; // 2.11 / 24 ;  //THETA2(7) / 24
      double TVMT2 = THETA2_7_divident/ THETA2_7_divisor;
      MT2 = TVMT2;

      fever_threshold = pt;

   } 


   WH_Dynamics::~WH_Dynamics()
   {        
   }

   void WH_Dynamics:: set_with_drug(bool status) {
      with_drug = status;
   }

   bool WH_Dynamics:: get_with_drug() {
      return with_drug;
   }

   void WH_Dynamics::step(int i, int agent_id, double &tt, int& mda1, int& mda2, int& mda3, 
                          double &kfeverj, double &x50abj, double &kabj, double &killabj,
                          matrix_dyn_t &vpk, matrix_dyn_t &vvpk, matrix_dyn_t &vpip, matrix_dyn_t &vvpip,
                          vector<double>& P0, vector<double>& AI, vector<double>& infect, vector<double> &MG, 
                          vector<double>& G, vector<double>& MSP, vector<double>& cyt, vector<double>& AB,
                          vector<int> &recieved_drug_timestamps, vector<int> &scheduled_drug_timestamps)
   {
      
      if (i <1) 
      {
            //P0(:,1)=4*10^4;
            //P0[0] = 4*pow(10,4); //nb of initial parasitaemia/ merozoit?  
            P0[0] = P0_initial_val;  //nb of initial parasitaemia /merozoit?  

            //tt=(j-1)*immunestack+1;  //in matlab j=level starts at 1
            tt=(((level+1)-1)*immunestack)+1.0;

            //std::cout << "i= " << i << ", level= "<< level<<", tt=  " << tt <<"\n";

            //P0(j,1)=max(round(normrnd(P0(j,1),20000) ),5*10^2);
            std::normal_distribution<> normrnd(P0[0], 20000.0);  
            //P0[0] = max(round(normrnd(r_engin)), 5*pow(10,2));  //original
            P0[0] = max(round(normrnd(r_engin)), 5*pow(10,3));  //changed to 10^3

            //P0[0] = 10000;  //rdn
            //P0[0] = 2431;  //rdn
            //P0[0] = 2004;  

            //std::cout << "P0[0] = "<< P0[0] <<"\n";
            //std::cout << "normrnd(r_engin):" << normrnd(r_engin)<<"\n";

            //if (P0[j,1] > maxinit)
            if (P0[0] > maxinit)  //j
               //P0(j,1)=maxinit;
               P0[0] = maxinit;

            //std::cout << "after > maxinit , P0[0] = "<< P0[0] <<"\n";

            //~~memory level for each person will depend on history of infection, to guarantee that individuals with 
            //~more past exposures are able to control initial parasitaemia more efficiently upon reinfection

            std::cout << std::fixed;
            std::cout << std::setprecision(10);
            // AI(j,1)=(j>1)* memory;   
            AI[0] = ((level>0)* memory);
            //std::cout << "AI[0] = " <<AI[0]<< " \n";

            /*
            vpk.resize(npk, asteps+1); //DHA/PK
            vvpk.resize(level+1, asteps+1); //[DHA/PK]

            vpip.resize(npip, asteps+1); //PIP/PD
            vvpip.resize(level+1, asteps+1); //[Pip/PD] 

            //vpk=zeros(npk,asteps+1);
            vpk.setZero();
            vvpk.setZero();

            //vpip=zeros(npip,asteps+1);
            vpip.setZero();

            //vvpip=zeros(levels,asteps+1);
            vvpip.setZero();   */

            //kfeverj=randn(1)*kfever*0.1+kfever;
            //kfeverj=randn(r_engin)*kfever*0.1+kfever;
            kfeverj=randn(r_engin)*kfever*kfever_multiplier+kfever;
            //kfeverj=0.5*kfever*0.1+kfever;  //rdn
            //std::cout << "kfeverj = "<< kfeverj <<"\n";


            //x50abj=randn(1)*x50ab*0.5+x50ab;
            //x50abj= randn(r_engin)*x50ab*0.5+x50ab;
            x50abj= randn(r_engin)*x50ab*x50ab_multiplier+x50ab;
            //x50abj= 0.5*x50ab*0.5+x50ab;  //rdn
            //std::cout << "x50abj = "<< x50abj <<"\n";

            //kabj=randn(1)*kab*0.02+kab;
            //kabj= randn(r_engin)*kab*0.02+kab;
            kabj= randn(r_engin)*kab*kab_multiplier+kab;
            //kabj= 0.5*kab*0.02+kab;  //rdn
            //std::cout << "kabj = "<< kabj <<"\n";

            //killabj=randn(1)*killab*0.05+killab;
            //~~~killab : adaptive immunity maximal kill rate
            //killabj= randn(r_engin)*killab*0.05+killab;
            killabj= randn(r_engin)*killab*killab_multiplier+killab;
            //killabj= 0.5*killab*0.05+killab;  //rdn
            //std::cout << "killabj = "<< killabj <<"\n";

      }
      else 
      {

         std::cout << std::fixed;
         std::cout << std::setprecision(10);
         // num=P0(j,i-1)*(1-AI(j,i-1));
         const double num=P0[i-1]*(1-AI[i-1]);

         //std::cout << "^^^^AI[i-1] " <<AI[i-1]<< " \n";
         //std::cout << "^^^^P0[i-1]= " <<P0[i-1] << " \n";
         //std::cout << "^^^num= " <<num << " \n";
         //std::cout << "^^b4^^cyt[i-1]= " <<cyt[i-1] << " \n";

         // cyt(j,i)=cyt(j,i-1)*exp(-1/gamma)+maxcyt*(num/(num+II50));
         cyt[i] = cyt[i-1]*exp(-1/gamma)+maxcyt*(num/(num+II50));

         if ( cyt[i] < 0) {
            std::cout << "cyt[i] < 0" <<"\n";
         }

         //std::cout << "^^aft^cyt[i]= " <<cyt[i] << " \n";

         // fever = kfeverj*cyt(j,i);
         const double fever = kfeverj*cyt[i];
         
         //std::cout << "^^fever= " << fever << " \n";

         // adaptive immunity
         //if (P0(j,i-1)>0)
         if (P0[i-1]>0) {

               //if (AI(j,i-1)<0.4) {
               if (AI[i-1]<0.4) {
                  //AI(j,i)= min(AI(j,i-1)+(1/etac)*(1-AI(j,i-1))*((P0(j,i-1)+x50abj*kabj)/(P0(j,i-1)+x50abj)),1);
                  if (with_drug)    
                     AI[i] = min(AI[i-1]+(1/etac)*(1-AI[i-1])*((P0[i-1]+x50abj*kabj)/(P0[i-1]+x50abj)),1.0);
                  else  AI[i] = AI[i-1]+(1/etac)*(1-AI[i-1])*((P0[i-1]+x50abj*kabj)/(P0[i-1]+x50abj));
                  //std::cout << std::fixed;
                  //std::cout << std::setprecision(10);
                  //std::cout << "AdIm if AI[i]= " << AI[i]  << " \n";

               }
               else {
                  //?etah=3.0; //Values of Ycapacity,i above 0.4 correspond to hyperimmunity, and Ycapacity,i then grows towards 1 with a time constant of ?hyperimmunity?=?3 days, regardless of antigen concentrations
                  //AI(j,i)= AI(j,i-1)+(1/etah)*(1-AI(j,i-1));
                  AI[i] = AI[i-1]+(1/etah)*(1-AI[i-1]);
                  //std::cout << "AdIm else AI[i]= " << AI[i]  << " \n";
               }

         }
         else {
               //if (AI(j,i-1)>memory)
               if (AI[i-1]>memory) {
                  //AI(j,i)= AI(j,i-1)-(1/decay)*(AI(j,i-1)-memory);
                  //~~ When parasites get cleared from circulation, Ycapacity decays towards Ymemory with a time constant τcapdecay,
                  //corresponding to loss of hyperimmunity approximately 120 days after disappearance of antigen
                  AI[i] = AI[i-1]-(1/decay)*(AI[i-1]-memory);
                  //std::cout << "***********AI[i] memory if = " << AI[i]  << " \n";

               }

               else {
                  //AI(j,i)= AI(j,i-1)-(1/decay)*(AI(j,i-1));
                  AI[i] = AI[i-1]-(1/decay)*(AI[i-1]);
                  //std::cout << "***********AI[i] memory else = " << AI[i]  << " \n";

               }

         }
         //if (P0(j,i-1)>0)
         if (P0[i-1]>0) {
               //if (AI(j,i-1)>0.3)
               if (AI[i-1]>0.3) {
                     //~~~~If Ycapacity >0.3, antibody is released, and antibody levels Yantibody rapidly approach Ycapacity
                     //~ with a time constant τabprod = 6 hours, in the presence of the parasite.
                  //AB(j,i)=AB(j,i-1)+(1/abprod)*(AI(j,i-1)-AB(j,i-1));
                  AB[i] = AB[i-1]+(1/abprod)*(AI[i-1]-AB[i-1]);
                  //std::cout << "AB[i] if = " <<AB[i] << " \n";
               }

               else {
                  //AB(j,i)=AB(j,i-1);
                  AB[i] = AB[i-1];
                  //std::cout << "AB[i] else = " <<AB[i] << " \n";

               }
         }
         else {
               //AB(j,i)=AB(j,i-1)-(1/abdecay)*AB(j,i-1);
               AB[i] = AB[i-1]-(1/abdecay)*AB[i-1];
               //std::cout << "AB[i] else2 = " <<AB[i] << " \n";
         }

         double rr =0;
         // MEROZOITES
         //if (i>2) {
         if (i>1) {
               //if ((i%2)==0) {
               if (((i+1)%2)==0) {

                  //~~~Upon schizont rupture, pmf*gam gametocytes are generated, whilst pmf*(1-gam) new merozoites join circulation
                  // #rr=pmf*(1-kgam)*(P0(j,i-2));
                  rr=pmf*(1-kgam)*(P0[i-2]);
                  //std::cout << "rr if= " <<rr << " \n";
                  //std::cout << "P0[i-2]= " <<P0[i-2] << " \n";
               }

               else {
                  rr=0;
                  //std::cout << "rr else= " <<rr << " \n";
               }
               
               //if (AI(j,i-1)>0.4)
               if (AI[i-1]>0.4) {
                  //MSP(j,i)= AI(j,i);
                  MSP[i] = AI[i];
                  //std::cout << "MSP[i] if = " <<MSP[i]<< " \n";
               }

               else {
                  //MSP(j,i)= MSP(j,i-1)+kmsp*(rr)/(rr+x50abj)*(1-MSP(j,i-1));
                  MSP[i] = MSP[i-1]+kmsp*(rr)/(rr+x50abj)*(1-MSP[i-1]);
                  //std::cout << "MSP[i] else = " <<MSP[i]<< " \n";
               }
         }

         if (with_drug) 
         {

            /*double TVCL= 3.2500;
            double TVV2= 5.3750;
            double TVMT= 0.0409;
            double TVV22= 69840.0000;
            double CL = TVCL;
            double V2 = TVV2;
            double MT = TVMT; */

            //double CL = 3.2500;
            //double V2 = 5.3750;
            //double MT = 0.0409;
            //double NN = 7.0; 

            /*
            This following transition matrix (Jpk) is replicated from the 
            original Matlab code. Although the key main values are 
            are read and set from the config.json file, at this moment, 
            it is assumed that transition matrices will have the current 
            format (~semi-hard-coded) which may not be correct. 
            @TODO:
            Verify how transition matrices should be built and make the 
            this transition matrix as flexible as possible
            */

            const double KTR = (NN+1.0)/MT;
            const double K13 = KTR;
            const double K20 = CL/V2;
            const double K34 = KTR;
            const double K45 = KTR;
            const double K56 = KTR;
            const double K67 = KTR;
            const double K78 = KTR;
            const double K89 = KTR;
            const double K92 = KTR;

            matrix_dyn_t Jpk(npk, npk);  //9x9 

            // state transition matrix
            Jpk <<
               -K13, 0, 0, 0, 0, 0, 0, 0, 0,
               0, -K20, 0, 0, 0, 0, 0, 0, K92,
               K13, 0, -K34, 0, 0, 0, 0, 0, 0,
               0, 0, K34, -K45, 0, 0, 0, 0, 0,
               0, 0, 0, K45, -K56, 0, 0, 0, 0,
               0, 0, 0, 0, K56, -K67, 0, 0, 0,
               0, 0, 0, 0, 0, K67, -K78, 0, 0,
               0, 0, 0, 0, 0, 0, K78, -K89, 0,
               0, 0, 0, 0, 0, 0, 0, K89, -K92;
            
            matrix_dyn_t ejpk = Jpk.exp();

            //std::cout << std::fixed;
            //std::cout << std::setprecision(0);

            //std::cout << "ejpk:\n" <<  ejpk << std::endl << std::endl; 

            //std::cout << "vpk.size:\n" <<  vpk.rows()<< " x "  << vpk.cols() <<std::endl << std::endl; 
            //std::cout << i<< ", vpk (before):\n" <<  vpk << std::endl << std::endl; 
            //std::cout << i-1<<", vpk.col(i-1):\n" <<   vpk.col(i-1) << std::endl << std::endl; 
            //std::cout << i << ", vpk.col(i):\n" <<   vpk.col(i) << std::endl << std::endl; 
            //std::cout << i<<", ejpk * vpk.col(i-1):\n" <<    ejpk * vpk.col(i-1) << std::endl << std::endl; 
            
            //  vpk(:,i)=ejpk*vpk(:,i-1);
            vpk.col(i) = ejpk * vpk.col(i-1);

            //std::cout << " ejpk * vpk.col(i-1):\n" <<   ejpk * vpk.col(i-1) << std::endl << std::endl; 
            //std::cout << "vpk (after):\n" <<  vpk << std::endl << std::endl; 
            //std::cout << agent_id <<",  t= " <<i <<"\n";
            //std::cout << " vpk.col(i):\n" <<   vpk.col(i) << std::endl << std::endl; 


            //-----
            /*double CL2 = TVCL2;
            double V22 = TVV22;
            double Q12 = TVQ12;
            double V32 = TVV32;
            double Q22 = TVQ22;
            double V42 = TVV42;
            double MT2 = TVMT2; */

            //double CL2 = 1329.6000;
            //double Q12 = 7440.0000;
            //double Q22 = 2520.0000;
            //double V22 = 69840.0000;
            //double V32 = 117840.0000;
            //double V42 = 741600.0000;
            //double MT2 = 0.0879;
            //double NN2 = 2.0;

            /*
            This following transition matrix (Jpip) is replicated from the 
            original Matlab code. Although the key main values are 
            are read and set from the config.json file, at this moment, 
            it is assumed that transition matrices will have the current 
            format (~semi-hard-coded) which may not be correct. 
            @TODO:
            Verify how transition matrices should be built and make the 
            this transition matrix as flexible as possible
            */

            const double KTR2 = (NN2+1)/MT2;
            const double K15 = KTR2;
            const double K20b = CL2/V22;  //K20b
            const double K23 = Q12/V22;
            const double K32 = Q12/V32;
            const double K24 = Q22/V22;
            const double K42 = Q22/V42;
            const double K56b = KTR2; // K56b 
            const double K62 = KTR2;

            matrix_dyn_t Jpip(npip, npip); //6 x 6

            Jpip <<-K15, 0, 0, 0, 0, 0,
                  0, (-K23-K20b-K24), K32, K42, 0, K62,
                  0, K23, -K32, 0, 0, 0,
                  0, K24, 0, -K42, 0, 0,
                  K15, 0, 0, 0, -K56b, 0,
                  0, 0, 0, 0, K56b, -K62;

            matrix_dyn_t ejpip = Jpip.exp();

            //std::cout << "ejpip:\n" <<  ejpip  << std::endl << std::endl; 
            //std::cout << "vpip (before):\n" <<  vpip << std::endl << std::endl; 

            //std::cout << " vpip.col(i):\n" <<  vpip.col(i) << std::endl << std::endl; 
            //std::cout << " vpip.col(i-1):\n" <<  vpip.col(i-1) << std::endl << std::endl; 
            //std::cout << " ejpip*vpip.col(i-1):\n" << ejpip*vpip.col(i-1) << std::endl << std::endl; 

            //vpip(:,i)=ejpip*vpip(:,i-1); 
            vpip.col(i)=ejpip*vpip.col(i-1); 

            //std::cout << "vpip (after):\n" <<  vpip << std::endl << std::endl; 



            // # treatment
            //# assign different dosing absorption rates per dose.
            if (i==mda1) {
                  //std::cout << "*i==mda1 " <<i<< ", mda1= "<< mda1 << " \n";
                  // vpk(1,i)=dosingpk*F1+vpk(1,i);
                  //vpk(1,i)=dosingpk+vpk(1,i);
                  vpk(0,i)=dosingpk+vpk(0,i);
                  //std::cout << "vpk(0,i) = " << vpk(0,i)<< " \n";

                  //vpip(1,i)=dosingpip*F12*0.577+vpip(1,i);
                  //vpip(1,i)=dosingpip*0.577+vpip(1,i);
                  vpip(0,i)=dosingpip*0.577+vpip(0,i);

                  //std::cout << "vpip(0,i) = " << vpip(0,i)<< " \n";
                 recieved_drug_timestamps.push_back(mda1);

            }

            if (i==mda2) {
                  //std::cout << "i==mda2 " <<i<< ", mda2= "<< mda2 << "\n";

                  //#vpk(1,i)=dosingpk*F1+vpk(1,i);
                  //vpk(1,i)=dosingpk+vpk(1,i);
                  vpk(0,i)=dosingpk+vpk(0,i);
                  //std::cout << "vpk(0,i) = " << vpk(0,i)<< " \n";

                  //vpip(1,i)=dosingpip*1.24*0.577+vpip(1,i); //F12=1
                  vpip(0,i)=dosingpip*1.24*0.577+vpip(0,i);
                  //std::cout << "vpip(0,i) = " << vpip(0,i)<< " \n";
                  recieved_drug_timestamps.push_back(mda2);

            }

            if (i==mda3) {

                  //std::cout << "i==mda3 " <<i<< ", mda3= "<< mda3 <<" \n";
            
                  //#vpk(1,i)=dosingpk*F1+vpk(1,i);
                  //vpk(1,i)=dosingpk+vpk(1,i); //F1=1
                  vpk(0,i)=dosingpk+vpk(0,i);
                  //std::cout << "vpk(0,i) = " << vpk(0,i)<< " \n";

                  //#vpip(1,i)=dosingpip*F12*(1.24^2)*0.577+vpip(1,i); //F12=1
                  vpip(0,i)=dosingpip*pow(1.24,2)*0.577+vpip(0,i);
                  //std::cout << "vpip(0,i) = " << vpip(0,i)<< " \n";
                  recieved_drug_timestamps.push_back(mda3);

            }

           //std::cout << "vpk.row(1): " << vpk.row(1) << std::endl; 

            //vvpk(j,:)=vpk(2,:);
            vvpk.row(level)=vpk.row(1);
           // std::cout << "vvpk(level,i):" <<  vvpk(level,i) << std::endl; 

            // vvpip(j,:)=vpip(2,:);
            vvpip.row(level)=vpip.row(1);
            //std::cout << "vvpip(level,i):" <<  vvpip(level,i) << std::endl; 
         
         } //(with_drug)

         //#Parasite Killing
         //killcyt: innate immune killing maximal rate
         //ky=killcyt*fever/(1+fever);
         const double ky=killcyt*fever/(1+fever);
         
         //std::cout << "ky= " <<ky<< " \n";

         //kz=killabj*AB(j,i);
         const double kz=killabj*AB[i];
         //std::cout << "^^^AB[i]= " <<AB[i]<< " \n";
         //std::cout << "kz= " <<kz<< " \n";

         //double killed = 1-exp(-(ky+kz));
      
      double kda = 0;
      double kdb = 0;

      if (with_drug) 
      {

         //# DRUGS PD
         //~~~~~~~ DHA kill rate Cdha=killratedha/(1+(ce50dha/[DHA])^h);
         //kda(w,j,i)=killratedha/(1+(ce50dha/vvpk(j,i))^hdha);
         //std::cout << "vvpk(level,i):" <<  vvpk(level,i) << std::endl; 
         kda = killratedha/ (1+ pow((ce50dha/ vvpk(level,i)),hdha)) ; //kda: drug effect/ DHA kill rate Cdha  ?

         //if (kda > 0) std::cout << "kda = " << kda<< " \n";

         //kdb(w,j,i)=killratepip/(1+(ce50pip/vvpip(j,i))^hpip);
         kdb =  killratepip/(1+ pow((ce50pip/vvpip(level,i)),hpip));

         //if (kdb> 0) std::cout << "kdb = " << kdb<< " \n";
         //std::cout << "vvpip(level,i) = " << vvpip(level,i)<< " \n";

      }
            
         //# Parasite Dynamics
         //# switching
         //#if P0(j,i-1)>1
         if (P0[i-1]>1) {
               //~~~Ycapacity = Ycapacity *[switchcount / (switcheffect+ switchcount)]
               //~~ as switchcount increases & becomes much larger than switcheffect, antibody capacity becomes insensitive to further switching events
               //#if (rand<switchrate){
               //std::cout << "rand(r_engin) = " <<rand(r_engin)<< " \n";
               if (rand(r_engin) < switchrate) {
               //if (0.02 < switchrate) {  //rdn

                  //#switch occurs
                  tt=tt+1;  //switchcount 
                  //#AI(j,i)=AI(j,i)*((tt/(switcheffect+tt)));
                  AI[i] = AI[i]*(tt/(switcheffect+tt));
                  //std::cout << "AI[i] switch= " <<AI[i]<< " \n";

                  //#AB(j,i)=AB(j,i)*((tt/(switcheffect+tt)));
                  AB[i] = AB[i]*(tt/(switcheffect+tt));
                  //std::cout << "AB[i] switch= " <<AB[i]<< " \n";
               }
         }

         // this was added because AI>1 value causes 'num' value to become negative which
         // in turns makes 'cyt' and eventually 'infect' become negative (<0)
         if (AI[i]> 1) {
             AI[i] = 1;
         }

         //# infected RBCs
         // if i>2 && mod(i,2)==0
         if ((i>1) && (((i+1)%2)==0)) { //#mod(i,2)

               //std::cout << "RBC, MSP[i-2]  if = " <<MSP[i-2] << " \n";

               //#P0(j,i)=P0(j,i-1)+(pmf*(1-kgam)*((P0(j,i-2))))*exp(-MSP(j,i-2)*Cmer*2); # survival of merozoites
               //~~~~Upon schizont rupture, pmf*gam gametocytes are generated, whilst pmf*(1-gam) new merozoites join circulation
               P0[i] = P0[i-1]+(pmf*(1-kgam)*((P0[i-2])))*exp(-MSP[i-2]*Cmer*2); //# survival of merozoites

               //P0(j,i)=(P0(j,i))*exp(-(ky+kz+kda(w,j,i)+kdb(w,j,i)));    

            if (with_drug)             
               P0[i] = P0[i]*exp(-(ky+kz+kda+kdb));  
            else P0[i] = P0[i]*exp(-(ky+kz));  

               //std::cout << "RBC, P0[i] if = " <<P0[i] << " \n";
         }
         else {

               //#P0(j,i)=(P0(j,i-1))*exp(-(ky+kz+kda(w,j,i)+kdb(w,j,i)));
               //std::cout << "RBC, P0[i-1] else = " <<P0[i-1] << " \n";

            if (with_drug)    
               P0[i] = (P0[i-1])*exp(-(ky+kz+kda+kdb));
            else P0[i] = (P0[i-1])*exp(-(ky+kz));
            
            //std::cout << "RBC, P0[i] else = " <<P0[i] << " \n";
         }

         //if (P0(j,i)<1) 
         if (P0[i]<1) {
               //P0[j,i]=0.000000001
               P0[i] = 0.000000001;
               //std::cout << "P0[i] = 0.000000001 " << " \n";
         }

         // # Gametocytes
         // if i>2 && mod(i,2)==0
         if ((i>1) && (((i+1)%2)==0)) {  // #(mod(i,2)
               //#G(j,i)=G(j,i-1)+(pmf*(kgam)*((P0(j,i-2))))*exp(-MSP(j,i-2)*Cmer*2); # survival of merozoites
               //~~~Upon schizont rupture, pmf*gam gametocytes are generated, whilst pmf*(1-gam) new merozoites join circulation
               //std::cout << "G[i-1]= " <<G[i-1]<< " \n";
               //std::cout << "`P0[i-2]= " <<P0[i-2]<< " \n";
               //std::cout << "`MSP[i-2]= " <<MSP[i-2]<< " \n";
               G[i] = G[i-1]+(pmf*(kgam)*((P0[i-2])))*exp(-MSP[i-2]*Cmer*2); //# survival of merozoites

               //#G(j,i)=(G(j,i))*exp(-(kg));
               //~~~ Gametocytes develop into mature gametocytes after a mean time of kg days which then live for a mean of kmg days
               //std::cout << "G[i]= " <<G[i]<< " \n";
               G[i] = G[i]*exp(-(kg));

               //std::cout << "Gametocytes, G[i] if = " <<G[i]<< " \n";
         }
         else {
               //#G(j,i)=(G(j,i-1))*exp(-(kg));
               //~~~ Gametocytes develop into mature gametocytes after a mean time of kg days which then live for a mean of kmg days
               //std::cout << "G[i-1]= " <<G[i-1]<< " \n";
               G[i] = G[i-1]*exp(-(kg));
               //std::cout << "Gametocytes, G[i] else = " <<G[i]<< " \n";
         }

         //if (G[j,i]<1) 
         if (G[i]<1) {

               //G[j,i]=0
               G[i] = 0;
               //std::cout << "Gametocytes, G[i] = 0 " << " \n";
         }

         //#MG(j,i)=MG(j,i-1)*exp(-(kmg))+(G(j,i-1))*(1-exp(-(kg)));
         //~~~ Gametocytes develop into mature gametocytes after a mean time of kg days which then live for a mean of kmg days
         MG[i] = MG[i-1]*exp(-(kmg))+(G[i-1])*(1-exp(-(kg)));
         //std::cout << " MG[i] = " << MG[i] << " \n";

         //if (MG[j,i]<1)
         if (MG[i]<1) {
               //MG[j,i]=0
               MG[i] = 0;
               //std::cout << " MG[i] =  0" << " \n";
         }

      //if (with_drug) 
      //{     
         //#FEVER OUTPUT
         //if (P0(j,i)>=pt) {
         if (P0[i]>=pt) {
               //std::cout << "FEVER OUTPUT  \n";
               //std::cout << "rand(r_engin) = " <<rand(r_engin)<< " \n";
               if (with_drug) {
                  //#if rand<0.33
                  if ( rand(r_engin) <0.33){ 
                  //if ( 0.2 <0.33){   //rdn
                     
                     //if mda1==0 || (i-mda1)>30
                     //if ((mda1==0) || (((i+1)-mda1)>30)) {
                     if ((mda1==0) || ((i-mda1)>30)) {

                        //# 'get drug'
                        //std::cout << "DDDDD get drug " << " \n";

                        mda1=i+2;
                        mda2=mda1+1;
                        mda3=mda2+1;

                        scheduled_drug_timestamps.push_back(mda1);
                        scheduled_drug_timestamps.push_back(mda2);
                        scheduled_drug_timestamps.push_back(mda3);
                     }
                     else {
                        //agent with fever did not recieve drug (becasue if (i-mad1>30 condition)
                        miss_drug_record record;
                        record.t = i;
                        record.id = agent_id;
                        record.p0 = P0[i];
                        record.rand_caused = 0;
                        miss_drug_record_vec.push_back(record);
                     }
                  }
                  else {
                     //agent with fever did not recieve drug (becasue rand<0.33 failed)
                     miss_drug_record record;
                     record.t = i;
                     record.id = agent_id;
                     record.p0 = P0[i];
                     record.rand_caused = 1;
                     miss_drug_record_vec.push_back(record);
                  }
               }
         }
      //}

         //#Infectivity
         //~~~~~Infectivity=(1-exp(-(MatureGam/(gam50))*kgametocyte*(1-Cyt/(Cyt+C50g))));
         //~~~~~mature gametocyte numbers are used to inform the infectivity function
         //~~ right hand side eq: ~ innate immunity interference on gametocyte viability 
         //kgametocyte: sucess rate per gametocyte in a mosquito ?

         //#infect(j,i)=(1-exp(-(MG(j,i)/(gam50))*kgametocyte*(1-cyt(j,i)/(cyt(j,i)+C50g))));
         infect[i] = ( 1-exp(- MG[i]/(gam50) *kgametocyte*  ( 1-cyt[i]  /    (cyt[i]+ C50g)  )  ) )  ;

         //std::cout << "infect[i] = " << infect[i] << " \n";
         //std::cout << "MG[i] = " << MG[i] << " \n";
         //std::cout << "cyt[i] = " << cyt[i] << " \n";

         //if (infect(j,i)<0) {
         if (infect[i]<0) {
               //stop("error:: stop");
               std::cout << "ERROR: infect[i]<0, stop" <<"\n";
               exit(0);
         }

      }
   }

   //==========================  WHD_Manager : 

   /*
    Within-Host Dynamics Manager takes care of managing the within-host dynamics, including 
    printing and recording the WHD related values. 
   */
   class WHD_Manager {

      private:

         const rapidjson::Document config;
         const std::string output_directory;
         const std::string output_prefix;
         const int sum_total_population;
         Human_Manager* human_mgr;
         const int num_time_steps = 0;
         
         //P0 metrics:
         const double p0_patent_indicator;
         const double p0_peak_divisor;
         const double infection_min_threshold;


         //const int whd_level;

         // constant output filenames: 
         const std::string whd_param_vals_export_file_name = "humans_whd_param_vals_daily";
         const std::string whd_param_avg_vals_export_file_name = "humans_whd_param_avg_vals_daily";
         const std::string p0_patents_export_file_name = "humans_P0_patents.csv";
         const std::string p0_peakpar_export_file_name = "humans_P0_peakpar.csv";
         const std::string p0_vals_export_file_name = "humans_P0_vals";

         const std::string mda_vals_export_file_name = "humans_mda_vals.csv";
         const std::string p0_patents_for_no_drug_recievers_export_file_name = "humans_P0_patents_for_no_drug_recievers.csv";

         const std::string p0_duration_export_file_name = "humans_P0_dur";

         const std::string recieved_drug_export_file_name = "humans_recieved_drug";
         const std::string scheduled_drug_export_file_name = "humans_scheduled_drug";
         const std::string agents_missed_drug_export_file_name = "humans_missed_drug";


      public: 
         //vector<HumanAgent> list_of_human_agents;
         WH_Dynamics whd;

      public:

         WHD_Manager(const std::string& configuration_file_name,
                     const std::string output_directory_name,
                     const std::string& output_prefix_name, 
                     const int initial_population, human::Human_Manager* h_mgr);

         ~WHD_Manager();

         void step_within_host_dynamics(int tick_time); 

         void print_WH_vals_for_all_humans() const;
         void print_human_WH_val(const int index_human) const;

         void print_whd_param_vals_for_all_humans() const;
         void print_whd_param_vals_for_human(const int index_human) const;

         int  calculate_P0_patent_for_human(HumanAgent agent) const;
         void print_P0_patent_for_all_humans() const;
         double  calculate_P0_peakpar_for_human(HumanAgent agent) const;

         void print_P0_duration_for_all_humans() const;
         int  calculate_P0_duration_for_human(HumanAgent agent) const;


         //to be called at the end of simulation: 
         void export_to_csv_whd_param_vals_for_all_humans(const std::string output_file_prefix) const;
         void export_to_csv_whd_param_avg_vals_of_all_humans(const std::string output_file_prefix) const;
         void export_to_csv_P0_patent_for_all_humans(const std::string output_file_prefix) const;
         void export_to_csv_P0_peakpar_for_all_humans(const std::string output_file_prefix) const;
         void export_to_csv_P0_vals_for_all_humans(const std::string output_file_prefix) const;

         void export_to_csv_total_MDA_vals(const std::string output_file_prefix) const;
         void export_to_cvs_P0_patent_for_no_drug_receivers(const std::string output_file_prefix) const;

         void export_to_csv_P0_duration_for_all_humans(const std::string output_file_prefix) const;

         void export_to_csv_recieved_drug_timestamps_for_all_humans(const std::string output_file_prefix) const;
         void export_to_csv_scheduled_drug_timestamps_for_all_humans(const std::string output_file_prefix) const;

         void export_to_csv_agents_info_missed_drug(const std::string output_file_prefix) const;

   };

   WHD_Manager::WHD_Manager(const std::string& configuration_file_name,
                            const std::string output_directory_name,
                            const std::string& output_prefix_name,
                            const int initial_population, human::Human_Manager* h_mgr) :

         config(util::get_json_from_file(configuration_file_name)),
         output_directory(output_directory_name),
         output_prefix(output_prefix_name.empty() ? util::get_output_prefix() : output_prefix_name),
         sum_total_population (initial_population),
         human_mgr(h_mgr),
         num_time_steps(config["simulation"]["total_steps"].GetInt()),
         p0_patent_indicator(config["whd"]["P0_related_metrics"]["p0_patent_indicator"].GetDouble()),
         p0_peak_divisor(config["whd"]["P0_related_metrics"]["p0_peak_divisor"].GetDouble()),
         infection_min_threshold(config["whd"]["P0_related_metrics"]["infection_min_threshold"].GetDouble()),


         whd(config["simulation"]["total_steps"].GetInt(), 
             config["whd"]["whd_level"].GetInt(),
             config["whd"]["II50"].GetInt(),
             config["whd"]["pmf"].GetDouble(),
             config["whd"]["gamma"].GetDouble(),
             config["whd"]["maxcyt"].GetDouble(),
             config["whd"]["kfever"].GetDouble(),
             config["whd"]["kfever_multiplier"].GetDouble(),
             config["whd"]["memory"].GetDouble(),
             config["whd"]["decay"].GetDouble(),
             config["whd"]["switchrate"].GetDouble(),
             config["whd"]["immunestack"].GetDouble(),
             config["whd"]["switcheffect"].GetDouble(),
             config["whd"]["maxinit"].GetDouble(),
             config["whd"]["etac"].GetDouble(),
             config["whd"]["etah"].GetDouble(),
             config["whd"]["x50ab"].GetUint64(),
             config["whd"]["x50ab_multiplier"].GetDouble(),
             config["whd"]["kab"].GetDouble(), 
             config["whd"]["kab_multiplier"].GetDouble(),       
             config["whd"]["abprod"].GetDouble(),
             config["whd"]["abdecay"].GetDouble(),
             config["whd"]["killcyt"].GetDouble(),
             config["whd"]["killab"].GetDouble(),
             config["whd"]["killab_multiplier"].GetDouble(),
             config["whd"]["kmsp"].GetDouble(),
             config["whd"]["Cmer"].GetDouble(),
             config["whd"]["kg_divisor"].GetDouble(),
             config["whd"]["kmg_divisor"].GetDouble(),
             config["whd"]["kgam"].GetDouble(),
             config["whd"]["kgametocyte"].GetDouble(),
             config["whd"]["pt"].GetUint64(),
             config["whd"]["C50g"].GetDouble(),
             config["whd"]["gam50"].GetDouble(),
             config["whd"]["dosingpk"].GetDouble(),
             config["whd"]["dosingpip"].GetDouble(),
             config["whd"]["killratedha"].GetDouble(),
             config["whd"]["killratepip"].GetDouble(),
             config["whd"]["ce50dha"].GetDouble(),
             config["whd"]["ce50pip"].GetDouble(),
             config["whd"]["hdha"].GetDouble(),
             config["whd"]["hpip"].GetDouble(),
             config["whd"]["npk"].GetInt(),
             config["whd"]["npip"].GetInt(),
             config["whd"]["p0_initial_val"].GetDouble(),
             config["whd"]["with_drug"].GetBool(),
             config["whd"]["transition_matrix"]["Jpk"]["CL"].GetDouble(),
             config["whd"]["transition_matrix"]["Jpk"]["V2"].GetDouble(),
             config["whd"]["transition_matrix"]["Jpk"]["THETA_3_divident"].GetDouble(),
             config["whd"]["transition_matrix"]["Jpk"]["THETA_3_divisor"].GetDouble(),
             config["whd"]["transition_matrix"]["Jpk"]["NN"].GetDouble(),
             config["whd"]["transition_matrix"]["Jpip"]["ACL2"].GetDouble(),
             config["whd"]["transition_matrix"]["Jpip"]["M_WE2"].GetDouble(),
             config["whd"]["transition_matrix"]["Jpip"]["AGE2"].GetDouble(),
             config["whd"]["transition_matrix"]["Jpip"]["WT2"].GetDouble(),
             config["whd"]["transition_matrix"]["Jpip"]["THETA2_1"].GetDouble(),
             config["whd"]["transition_matrix"]["Jpip"]["EM502"].GetDouble(),
             config["whd"]["transition_matrix"]["Jpip"]["HILL2"].GetDouble(),
             config["whd"]["transition_matrix"]["Jpip"]["Q12"].GetDouble(),
             config["whd"]["transition_matrix"]["Jpip"]["Q22"].GetDouble(),
             config["whd"]["transition_matrix"]["Jpip"]["V22"].GetDouble(),
             config["whd"]["transition_matrix"]["Jpip"]["V32"].GetDouble(),
             config["whd"]["transition_matrix"]["Jpip"]["V42"].GetDouble(),
             config["whd"]["transition_matrix"]["Jpip"]["THETA2_7_divident"].GetDouble(),
             config["whd"]["transition_matrix"]["Jpip"]["THETA2_7_divisor"].GetDouble(),
             config["whd"]["transition_matrix"]["Jpip"]["NN2"].GetDouble()
             
         )
      
      {

         std::cout << "simulate for length of: " << num_time_steps 
         << " time steps, with number of agents: " << sum_total_population << "\n";
         std::cout << "\n";

         /*for (int i=1; i<=sum_total_population; i++) {
            list_of_human_agents.push_back(HumanAgent(i, num_time_steps));
         } */

      }


   WHD_Manager::~WHD_Manager() {}

   void WHD_Manager::step_within_host_dynamics(int timetick)
   {

      for (int hh=0; hh<sum_total_population; hh++) {

         //HumanAgent &agent = list_of_human_agents.at(i);
         HumanAgent &agent = human_mgr->list_of_human_agents.at(hh);

         //std::cout << "=============== t= "<< timetick<< ", step for agent # " << agent.id<<" : \n";

         whd.step(timetick, agent.id, agent.tt, agent.mda1, agent.mda2, agent.mda3, 
                     agent.kfeverj, agent.x50abj, agent.kabj, agent.killabj,
                     agent.vpk, agent.vvpk, agent.vpip, agent.vvpip,
                     agent.P0, agent.AI, agent.infect, agent.MG,
                     agent.G, agent.MSP, agent.cyt, agent.AB,
                     agent.recieved_drug_timestamps, agent.scheduled_drug_timestamps);
		}
   }

   //-------------
   void WHD_Manager::print_WH_vals_for_all_humans() const{
      //std::cout << " HumanManager::print_all(), total_population= " <<total_population<<"\n";
      for (int i=0; i<sum_total_population; i++) {
			//HumanAgent h_agent = list_of_human_agents.at(i);
         HumanAgent h_agent = human_mgr->list_of_human_agents.at(i);

         h_agent.print_human_WH_params_values();
		}
   }

   void WHD_Manager::print_human_WH_val(const int human_index) const{
			//HumanAgent h_agent = list_of_human_agents.at(human_index);
         HumanAgent h_agent = human_mgr->list_of_human_agents.at(human_index);
         h_agent.print_human_WH_params_values();
   }


  
   void WHD_Manager::export_to_csv_whd_param_vals_for_all_humans(const std::string output_file_prefix) const
   {

      int level_n =  config["whd"]["whd_level"].GetInt();
      std::string level_s = std::to_string(level_n);

      std::cout << "Export human(s) whd params values to csv file: " << output_file_prefix+ whd_param_vals_export_file_name + "_L"+ level_s+".csv" << "\n";

      util::CSVWriter csv_writer(output_file_prefix+ whd_param_vals_export_file_name + "_L"+ level_s+".csv");

      if (sum_total_population <= 1)
      {
         std::vector<std::string> header = {"AI", "P0", "infect", "MG", "G", "MSP", "cyt", "AB"};
         csv_writer.addHeader(header);

         std::vector<std::string> vals;

         HumanAgent agent = human_mgr->list_of_human_agents.at(0);

         for (int t = 0; t < num_time_steps; t++)
         {

            vals.push_back(to_string(agent.AI[t]));
            vals.push_back(to_string(agent.P0[t]));
            vals.push_back(to_string(agent.infect[t]));
            vals.push_back(to_string(agent.MG[t]));
            vals.push_back(to_string(agent.G[t]));
            vals.push_back(to_string(agent.MSP[t]));
            vals.push_back(to_string(agent.cyt[t]));
            vals.push_back(to_string(agent.AB[t]));

            csv_writer.addRow(vals);
            vals.clear();
         }
      }

      else
      {
         std::vector<std::string> header;
         std::vector<std::string> headerNames = {"AI", "P0", "infect", "MG", "G", "MSP", "cyt", "AB"};

         for (int hh = 0; hh < sum_total_population; hh++)
         {
            HumanAgent agent = human_mgr->list_of_human_agents.at(hh);
            for (unsigned int i = 0; i < headerNames.size(); i++)
            {
               // cout << headerNames.at(i) <<", ";
               header.push_back(headerNames.at(i) + "_" + std::to_string(agent.getID()));
            }
         }

         /*for (int aa=1; aa<=total_population; aa++) {
            for (unsigned int i=0; i<headerNames.size(); i++) {
               //cout << headerNames.at(i) <<", ";
               header.push_back(headerNames.at(i)+"_"+std::to_string(aa));
            }
         } */

         csv_writer.addHeader(header);

         std::vector<std::string> AI_vals;
         std::vector<std::string> P0_vals;
         std::vector<std::string> infect_vals;
         std::vector<std::string> MG_vals;
         std::vector<std::string> G_vals;
         std::vector<std::string> MSP_vals;
         std::vector<std::string> cyt_vals;
         std::vector<std::string> AB_vals;

         vector<vector<std::string>> agent_vals;
         vector<vector<vector<std::string>>> vals;

         for (int i = 0; i < sum_total_population; i++)
         {

            HumanAgent agent = human_mgr->list_of_human_agents.at(i);

            for (int t = 0; t < num_time_steps; t++)
            {
               AI_vals.push_back(to_string(agent.AI[t]));
               P0_vals.push_back(to_string(agent.P0[t]));
               infect_vals.push_back(to_string(agent.infect[t]));
               MG_vals.push_back(to_string(agent.MG[t]));
               G_vals.push_back(to_string(agent.G[t]));
               MSP_vals.push_back(to_string(agent.MSP[t]));
               cyt_vals.push_back(to_string(agent.cyt[t]));
               AB_vals.push_back(to_string(agent.AB[t]));
            }

            agent_vals.push_back(AI_vals);
            agent_vals.push_back(P0_vals);
            agent_vals.push_back(infect_vals);
            agent_vals.push_back(MG_vals);
            agent_vals.push_back(G_vals);
            agent_vals.push_back(MSP_vals);
            agent_vals.push_back(cyt_vals);
            agent_vals.push_back(AB_vals);

            AI_vals.clear();
            P0_vals.clear();
            infect_vals.clear();
            MG_vals.clear();
            G_vals.clear();
            MSP_vals.clear();
            cyt_vals.clear();
            AB_vals.clear();

            vals.push_back(agent_vals);
            agent_vals.clear();
         }

         /*for (unsigned int i=0; i<vals.size(); i++) { //number of agents
            std::cout << "agent" << i << ": \n ";
            for (unsigned int j=0; j<vals.at(i).size(); j++) { //each param (AI, P0, ...)
               std::cout << "\n";
               for (unsigned int k=0; k<vals.at(i).at(j).size(); k++) {  //timesteps
                  std::cout << vals[i][j][k] << " ";
               }
               std::cout << "\n";
            }
         } */

         // write in CSV:
         std::vector<std::string> line_vals;

         for (int k = 0; k < num_time_steps; k++)
         { // timesteps
            for (unsigned int i = 0; i < vals.size(); i++)
            { // number of agents
               for (unsigned int j = 0; j < vals.at(i).size(); j++)
               { // each param (AI, P0, ...)
                  // std::cout << vals[i][j][k] << " ";
                  line_vals.push_back(vals[i][j][k]);
               }
            }
            // std::cout << "\n \n";
            csv_writer.addRow(line_vals);
            line_vals.clear();
         }
      }
   }

   void WHD_Manager::export_to_csv_whd_param_avg_vals_of_all_humans(const std::string output_file_prefix) const
   {

      int level_n =  config["whd"]["whd_level"].GetInt();
      std::string level_s = std::to_string(level_n);

      std::cout << "Export human(s) whd param AVG values to csv file: " << output_file_prefix+ whd_param_avg_vals_export_file_name + "_L" + level_s + ".csv" << "\n";

      util::CSVWriter csv_writer(output_file_prefix + whd_param_avg_vals_export_file_name + "_L" + level_s + ".csv" );

      if (sum_total_population <= 1)
      {
         std::vector<std::string> header = {"AI_avg", "P0_avg", "infect_avg", "MG_avg", "G_avg", "MSP_avg", "cyt_avg", "AB_avg"};
         csv_writer.addHeader(header);

         std::vector<std::string> vals;

         HumanAgent agent = human_mgr->list_of_human_agents.at(0);

         for (int t = 0; t < num_time_steps; t++)
         {

            vals.push_back(to_string(agent.AI[t]));
            vals.push_back(to_string(agent.P0[t]));
            vals.push_back(to_string(agent.infect[t]));
            vals.push_back(to_string(agent.MG[t]));
            vals.push_back(to_string(agent.G[t]));
            vals.push_back(to_string(agent.MSP[t]));
            vals.push_back(to_string(agent.cyt[t]));
            vals.push_back(to_string(agent.AB[t]));

            csv_writer.addRow(vals);
            vals.clear();
         }
      }

      else
      {
         std::vector<std::string> header;
         std::vector<std::string> headerNames = {"AI_avg", "P0_avg", "infect_avg", "MG_avg", "G_avg", "MSP_avg", "cyt_avg", "AB_avg"};

         for (unsigned int i = 0; i < headerNames.size(); i++)
         {

            header.push_back(headerNames.at(i));
         }

         csv_writer.addHeader(header);

         vector<std::string> agents_avg_vals;

         double AI_avg = 0;
         double P0_avg = 0;
         double infect_avg = 0;
         double MG_avg = 0;
         double G_avg = 0;
         double MSP_avg = 0;
         double cyt_avg = 0;
         double AB_avg = 0;

         for (int t = 0; t < num_time_steps; t++)
         {
            AI_avg = 0;
            P0_avg = 0;
            infect_avg = 0;
            MG_avg = 0;
            G_avg = 0;
            MSP_avg = 0;
            cyt_avg = 0;
            AB_avg = 0;

            for (int i = 0; i < sum_total_population; i++)
            {

               HumanAgent agent = human_mgr->list_of_human_agents.at(i);

               AI_avg = AI_avg + agent.AI[t];
               P0_avg = P0_avg + agent.P0[t];
               infect_avg = infect_avg + agent.infect[t];
               MG_avg = MG_avg + agent.MG[t];
               G_avg = G_avg + agent.G[t];
               MSP_avg = MSP_avg + agent.MSP[t];
               cyt_avg = cyt_avg + agent.cyt[t];
               AB_avg = AB_avg + agent.AB[t];
            }

            agents_avg_vals.push_back(to_string(AI_avg / sum_total_population));
            agents_avg_vals.push_back(to_string(P0_avg / sum_total_population));
            agents_avg_vals.push_back(to_string(infect_avg / sum_total_population));
            agents_avg_vals.push_back(to_string(MG_avg / sum_total_population));
            agents_avg_vals.push_back(to_string(G_avg / sum_total_population));
            agents_avg_vals.push_back(to_string(MSP_avg / sum_total_population));
            agents_avg_vals.push_back(to_string(cyt_avg / sum_total_population));
            agents_avg_vals.push_back(to_string(AB_avg / sum_total_population));

            csv_writer.addRow(agents_avg_vals);

            agents_avg_vals.clear();
         }
      }
   }

   void WHD_Manager::export_to_csv_P0_vals_for_all_humans(const std::string output_file_prefix) const
   {

      int level_n =  config["whd"]["whd_level"].GetInt();
      std::string level_s = std::to_string(level_n);

      std::cout << "Export human(s) P0-values to csv file: " << output_file_prefix+ p0_vals_export_file_name + "_L" + level_s + ".csv" << "\n";

      util::CSVWriter csv_writer(output_file_prefix+ p0_vals_export_file_name + "_L" + level_s + ".csv" );
      if (sum_total_population <= 1)
      {
         std::vector<std::string> vals;

         HumanAgent agent = human_mgr->list_of_human_agents.at(0);

         for (int t = 0; t < num_time_steps; t++)
         {
            vals.push_back(to_string(agent.P0[t]));
            csv_writer.addRow(vals);
            vals.clear();
         }
      }

      else
      {

         std::vector<std::string> P0_vals;
         vector<vector<std::string>> agent_vals;
         vector<vector<vector<std::string>>> vals;

         for (int i = 0; i < sum_total_population; i++)
         {

            HumanAgent agent = human_mgr->list_of_human_agents.at(i);

            for (int t = 0; t < num_time_steps; t++)
            {
               P0_vals.push_back(to_string(agent.P0[t]));
            }

            agent_vals.push_back(P0_vals);
            P0_vals.clear();

            vals.push_back(agent_vals);
            agent_vals.clear();
         }

         // write in CSV:
         std::vector<std::string> line_vals;

         for (int k = 0; k < num_time_steps; k++)
         { // timesteps
            for (unsigned int i = 0; i < vals.size(); i++)
            { // number of agents
               for (unsigned int j = 0; j < vals.at(i).size(); j++)
               { // each param (AI, P0, ...)
                  line_vals.push_back(vals[i][j][k]);
               }
            }
            csv_writer.addRow(line_vals);
            line_vals.clear();
         }
      }
   }


   /*
   A simple method to calculate the infection duration (days where parasite is above infection_min_threshold)
   */
   int WHD_Manager::calculate_P0_duration_for_human(HumanAgent agent) const
   {

      int dur = 0;

      for (int t = 0; t < num_time_steps; t++)
      {

         if (agent.P0[t] > infection_min_threshold)
         {
            dur++;
         }
         else {
            break;  //if we want the count of days stops as soon as parasites falls bellow min val
         }
      }

      return dur;
   }


   /*
   A simple method to calculate the parasitemia patent (sum/total of days where P0 is higher than a specific value for a given agent)
   Each iteration/ day when agent P0 value is above 'p0_patent_indicator', a day is added to the parasitemia patent.
   */
   int WHD_Manager::calculate_P0_patent_for_human(HumanAgent agent) const
   {

      int patent = 0;

      for (int t = 0; t < num_time_steps; t++)
      {

         if (agent.P0[t] > p0_patent_indicator)
         {
            patent++;
         }
      }
      return patent;
   }

   void WHD_Manager::print_P0_patent_for_all_humans() const
   {

      for (int i = 0; i < sum_total_population; i++)
      {
         HumanAgent agent = human_mgr->list_of_human_agents.at(i);
         int patent = calculate_P0_patent_for_human(agent);
         std::cout << "P0 patent for agent #" << agent.id << " = " << patent << "\n";
      }
   } 

   void WHD_Manager::print_P0_duration_for_all_humans() const
   {

      for (int i = 0; i < sum_total_population; i++)
      {
         HumanAgent agent = human_mgr->list_of_human_agents.at(i);
         int dur = calculate_P0_duration_for_human(agent);
         std::cout << "P0 duration for agent #" << agent.id << " = " << dur << "\n";
         std::cout << " -------\n";
      }
   } 

   void WHD_Manager::export_to_csv_P0_patent_for_all_humans(const std::string output_file_prefix) const
   {

      std::cout << "Export human(s) P0 patent values to csv file: " << output_file_prefix + p0_patents_export_file_name << "\n";

      util::CSVWriter csv_writer(output_file_prefix + p0_patents_export_file_name);

      std::vector<std::string> header = {"agentID", "P0_patent"};
      csv_writer.addHeader(header);

      std::vector<std::string> vals;

      for (int i = 0; i < sum_total_population; i++)
      {

         HumanAgent agent = human_mgr->list_of_human_agents.at(i);

         int patent = 0;

         for (int t = 0; t < num_time_steps; t++)
         {
            patent = calculate_P0_patent_for_human(agent);
         }

         vals.push_back(to_string(agent.id));
         vals.push_back(to_string(patent));
         csv_writer.addRow(vals);
         vals.clear();
      }
   } 

   /*
   A simple method to first find the peak of parasitemia (max of P0 among all days)
   This max/peak value is then divided by a specific value ('p0_peak_divisor') and returned.
   */
   double WHD_Manager::calculate_P0_peakpar_for_human(HumanAgent agent) const
   {
      double peak_value = -1;

      for (int t = 0; t < num_time_steps; t++)
      {
         if (agent.P0[t] > peak_value)
         {
            peak_value = agent.P0[t];
         }
      }

      return peak_value / p0_peak_divisor;
   }

   void WHD_Manager::export_to_csv_P0_peakpar_for_all_humans(const std::string output_file_prefix) const
   {

      std::cout << "Export human(s) P0 peakpar values to csv file: " << output_file_prefix + p0_peakpar_export_file_name << "\n";

      util::CSVWriter csv_writer(output_file_prefix + p0_peakpar_export_file_name);

      std::vector<std::string> header = {"agentID", "P0_peakpar"};
      csv_writer.addHeader(header);

      std::vector<std::string> vals;

      for (int i = 0; i < sum_total_population; i++)
      {

         HumanAgent agent = human_mgr->list_of_human_agents.at(i);

         double peakpar = 0;

         for (int t = 0; t < num_time_steps; t++)
         {
            peakpar = calculate_P0_peakpar_for_human(agent);
         }

         vals.push_back(to_string(agent.id));
         vals.push_back(to_string(peakpar));
         csv_writer.addRow(vals);
         vals.clear();
      }
   }


   void WHD_Manager::export_to_csv_P0_duration_for_all_humans(const std::string output_file_prefix) const
   {

      int level_n =  config["whd"]["whd_level"].GetInt();
      std::string level_s = std::to_string(level_n);
      std::cout << "level_s= " << level_s<<"\n";

      std::cout << "Export human(s) P0 duration values to csv file: " << output_file_prefix + p0_duration_export_file_name + "_L"+level_s +".csv" << "\n";
   
      util::CSVWriter csv_writer(output_file_prefix + p0_duration_export_file_name + "_L"+ level_s +".csv" );

      std::vector<std::string> header = {"agentID", "dur", "treatment_count", "clinical_count"};
      csv_writer.addHeader(header);

      std::vector<std::string> vals;

      for (int i = 0; i < sum_total_population; i++)
      {

         HumanAgent agent = human_mgr->list_of_human_agents.at(i);

         double dur = 0;
         int clinical_count =0;

         for (int t = 0; t < num_time_steps; t++)
         {
            dur = calculate_P0_duration_for_human(agent);
            
            if (agent.P0[t] >= whd.fever_threshold) {
               clinical_count++;
            }
         }

         vals.push_back(to_string(agent.id));
         vals.push_back(to_string(dur));

         int treatment_count = 0;
         int nb_of_drug_recieved = agent.recieved_drug_timestamps.size();

         if (nb_of_drug_recieved > 0) {
            treatment_count = nb_of_drug_recieved/3;
         }

         vals.push_back(to_string(treatment_count));
         vals.push_back(to_string(clinical_count));

         csv_writer.addRow(vals);
         vals.clear();
      }
   }

    void WHD_Manager::export_to_csv_recieved_drug_timestamps_for_all_humans(const std::string output_file_prefix) const{

      int level_n =  config["whd"]["whd_level"].GetInt();
      std::string level_s = std::to_string(level_n);
      std::cout << "level_s= " << level_s<<"\n";

        std::string end_file_name = ".csv";

        std::string file_name = this->recieved_drug_export_file_name + "_L"+ level_s + end_file_name; 

        std::cout << "Export recieved drug timestamps to csv file: " << output_file_prefix+file_name <<"\n";
        util::CSVWriter  csv_writer(output_file_prefix + file_name);

        std::vector<std::string> header;
        //std::vector<std::string> headerNames = {"infect_times", "bite_times"};
        std::vector<std::string> headerNames = {"a"};

        //for (int hh=0; hh<sum_total_population; hh++) {
        for (int hh=1; hh<=sum_total_population; hh++) {

            for (unsigned int i=0; i<headerNames.size(); i++) {
                header.push_back(headerNames.at(i)+"_"+std::to_string(hh));  //human agents do not have ID yet?!
            }
        }

        csv_writer.addHeader(header);

        int max_nb_of_recieved_drugs = -1;
        for (int i=0; i< sum_total_population; i++) {

            HumanAgent agent = human_mgr->list_of_human_agents.at(i);

            int nb_of_recieved_drugs = agent.recieved_drug_timestamps.size(); 

            if (nb_of_recieved_drugs > max_nb_of_recieved_drugs ) {
                max_nb_of_recieved_drugs = agent.recieved_drug_timestamps.size();
            }
        }

        std::vector<std::string> recieved_drug_time_vals;
        std::vector<std::vector<std::string>> agent_vals;
        std::vector<std::vector<std::vector<std::string>>> vals;

        for (int i=0; i< sum_total_population; i++) {
                     
            HumanAgent agent = human_mgr->list_of_human_agents.at(i);

            int nb_of_recieved_drugs = agent.recieved_drug_timestamps.size();

            for (int k=0; k< nb_of_recieved_drugs; k++) {

                recieved_drug_time_vals.push_back(std::to_string(agent.recieved_drug_timestamps[k]));
            }

            for (int j=0; j< (max_nb_of_recieved_drugs - nb_of_recieved_drugs); j++) {
                recieved_drug_time_vals.push_back("");
            }

            agent_vals.push_back(recieved_drug_time_vals);
            recieved_drug_time_vals.clear();

            vals.push_back(agent_vals);
            agent_vals.clear();
         }

        //write in CSV:
        //std::cout << "write in CSV: \n";
        std::vector<std::string> line_vals;
        for (int k=0; k<max_nb_of_recieved_drugs; k++) {  //timesteps
            for (unsigned int i=0; i<vals.size(); i++) { //number of agents
                for (unsigned int j=0; j<vals.at(i).size(); j++) { //each param (AI, P0, ...)
                    //std::cout << vals[i][j][k] << " ";
                    line_vals.push_back(vals[i][j][k]);
                }				
            }
            //std::cout << "\n \n";
            csv_writer.addRow(line_vals);
            line_vals.clear();
         }          
    } 

   void WHD_Manager::export_to_csv_scheduled_drug_timestamps_for_all_humans(const std::string output_file_prefix) const{

      int level_n =  config["whd"]["whd_level"].GetInt();
      std::string level_s = std::to_string(level_n);
      std::cout << "level_s= " << level_s<<"\n";

        std::string end_file_name = ".csv";

        std::string file_name = this->scheduled_drug_export_file_name + "_L"+ level_s + end_file_name; 

        std::cout << "Export scheduled drug timestamps to csv file: " << output_file_prefix+file_name <<"\n";
        util::CSVWriter  csv_writer(output_file_prefix + file_name);

        std::vector<std::string> header;
        //std::vector<std::string> headerNames = {"infect_times", "bite_times"};
        std::vector<std::string> headerNames = {"a"};

        //for (int hh=0; hh<sum_total_population; hh++) {
         for (int hh=1; hh<=sum_total_population; hh++) {

            for (unsigned int i=0; i<headerNames.size(); i++) {
                header.push_back(headerNames.at(i)+"_"+std::to_string(hh));  //human agents do not have ID yet?!
            }
        }

        csv_writer.addHeader(header);

        int max_nb_of_scheduled_drugs = -1;
        for (int i=0; i< sum_total_population; i++) {

            HumanAgent agent = human_mgr->list_of_human_agents.at(i);

            int nb_of_scheduled_drugs = agent.scheduled_drug_timestamps.size(); 

            if (nb_of_scheduled_drugs > max_nb_of_scheduled_drugs ) {
                max_nb_of_scheduled_drugs = agent.scheduled_drug_timestamps.size();
            }
        }

        std::vector<std::string> scheduled_drug_time_vals;
        std::vector<std::vector<std::string>> agent_vals;
        std::vector<std::vector<std::vector<std::string>>> vals;

        for (int i=0; i< sum_total_population; i++) {
                     
            HumanAgent agent = human_mgr->list_of_human_agents.at(i);

            int nb_of_scheduled_drugs = agent.recieved_drug_timestamps.size();

            for (int k=0; k< nb_of_scheduled_drugs; k++) {

                scheduled_drug_time_vals.push_back(std::to_string(agent.scheduled_drug_timestamps[k]));
            }

            for (int j=0; j< (max_nb_of_scheduled_drugs - nb_of_scheduled_drugs); j++) {
                scheduled_drug_time_vals.push_back("");
            }

            agent_vals.push_back(scheduled_drug_time_vals);
            scheduled_drug_time_vals.clear();

            vals.push_back(agent_vals);
            agent_vals.clear();
        }

        //write in CSV:
        //std::cout << "write in CSV: \n";
        std::vector<std::string> line_vals;
        for (int k=0; k<max_nb_of_scheduled_drugs; k++) {  //timesteps
            for (unsigned int i=0; i<vals.size(); i++) { //number of agents
                for (unsigned int j=0; j<vals.at(i).size(); j++) { //each param (AI, P0, ...)
                    //std::cout << vals[i][j][k] << " ";
                    line_vals.push_back(vals[i][j][k]);
                }				
            }
            //std::cout << "\n \n";
            csv_writer.addRow(line_vals);
            line_vals.clear();
        }          
   } 

   void WHD_Manager::export_to_csv_agents_info_missed_drug(const std::string output_file_prefix) const
   {

      int level_n =  config["whd"]["whd_level"].GetInt();
      std::string level_s = std::to_string(level_n);
      std::cout << "level_s= " << level_s<<"\n";

      std::cout << "Export humans info missed drug to csv file: " << output_file_prefix + agents_missed_drug_export_file_name + "_L"+level_s +".csv" << "\n";
   
      util::CSVWriter csv_writer(output_file_prefix + agents_missed_drug_export_file_name + "_L"+ level_s +".csv" );

      std::vector<std::string> header = {"agentID", "time", "p0", "rand_caused"};
      csv_writer.addHeader(header);

      std::vector<std::string> vals;

      int nb_of_agents_missed_drug = this->whd.miss_drug_record_vec.size();

      for (int i = 0; i < nb_of_agents_missed_drug; i++)
      {

         vals.push_back(to_string(whd.miss_drug_record_vec[i].id));
         vals.push_back(to_string(whd.miss_drug_record_vec[i].t));
         vals.push_back(to_string(whd.miss_drug_record_vec[i].p0));
         vals.push_back(to_string(whd.miss_drug_record_vec[i].rand_caused));

         csv_writer.addRow(vals);
         vals.clear();
      }
   }

   
   //++++++++++++++++
   
   /*
   Export the number of doses each agent has received during 3 mass drug administration (MDA) campaign
   */
   void WHD_Manager::export_to_csv_total_MDA_vals(const std::string output_file_prefix) const
   {

      std::cout << "Export human(s) total MDA values to csv file: " << output_file_prefix+ mda_vals_export_file_name << "\n";

      util::CSVWriter csv_writer(output_file_prefix+ mda_vals_export_file_name );

      std::vector<std::string> header = {"mda_1", "mda_2", "mda_3"};
      csv_writer.addHeader(header);

      vector<std::string> agenst_mda_vals;

      for (int i = 0; i < sum_total_population; i++)
      {

         //HumanAgent agent = list_of_human_agents.at(i);
         HumanAgent agent = human_mgr->list_of_human_agents.at(i);

         agenst_mda_vals.push_back(to_string(agent.mda1));
         agenst_mda_vals.push_back(to_string(agent.mda2));
         agenst_mda_vals.push_back(to_string(agent.mda3));

         csv_writer.addRow(agenst_mda_vals);

         agenst_mda_vals.clear();
      }
   }

   
   void  WHD_Manager::export_to_cvs_P0_patent_for_no_drug_receivers(const std::string output_file_prefix) const
   {

      std::cout << "Export human(s) P0 patent values for those who did not recieve any drug during 3 MDAs to csv file: " << output_file_prefix + p0_patents_for_no_drug_recievers_export_file_name<< "\n";

      util::CSVWriter csv_writer(output_file_prefix + p0_patents_for_no_drug_recievers_export_file_name);

      std::vector<std::string> header = {"agentID", "P0_patent"};
      csv_writer.addHeader(header);

      std::vector<std::string> vals;

      for (int i = 0; i < sum_total_population; i++)
      {

         HumanAgent agent = human_mgr->list_of_human_agents.at(i);

         if ((agent.mda1 == 0) && (agent.mda2 == 0) && (agent.mda3 == 0))
         {

            int patent = 0;

            for (int t = 0; t < num_time_steps; t++)
            {
               patent = calculate_P0_patent_for_human(agent);
            }

            vals.push_back(to_string(agent.id));
            vals.push_back(to_string(patent));
            csv_writer.addRow(vals);
            vals.clear();
         }
      }
   }








} //common 

#endif