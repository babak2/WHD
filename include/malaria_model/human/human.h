#ifndef HUMAN_H
#define HUMAN_H

#include "Eigen/Dense"
#include "Eigen-unsupported/Eigen/MatrixFunctions" 
//#include <cassert>
//#include <vector>

using namespace std;

namespace human {

/* 
   HumanAgent class is a temporary test class representing a simplified version of the base Human class
   in the base malaria model. In other words, the human agents in the base malaria model should be used instead
   of this class in order to use the functionaliites of WH_Dynamics classes. 
   
   To do so the humans (agents) in the base malaria_model should implement the following data strcutures 
   (i.e. AI, PO, MG, ...) and initialize them correspondingl (like here in HumanAgent).

   Note that similar to the most human parameters in the base malaria_model, these data strucutres were 
   initially of type double*, but have been recently changed to vector<double>

   */
   class HumanAgent {

      private: 

         int total_steps = 0;  //simulation time
		 typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_dyn_t;
         typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::DontAlign> vector_dyn_t; //vertical vector


      public:

         int id;


		 /*these vectors are initialized to 0 at the begining but their values will be changed
		 (i.e. they are not fixed) since in each timestep they will be replaced/updated by
		  new calculated values. 
		 */
		
		 double tt=0;
		 int mda1=0;  
         int mda2=0;  
         int mda3=0; 
		
		 double kfeverj=0; 
         double x50abj=0;
         double kabj= 0;  
         double killabj= 0; 

		 matrix_dyn_t vpk;
         matrix_dyn_t vvpk;

         matrix_dyn_t vpip; 
         matrix_dyn_t vvpip;
         
    	 std::vector<double> AI;   // Y[capacity] //Antibody capacity ?
         std::vector<double> P0;   // Total Parasitaemia ~ Par ~ gg
         std::vector<double> infect;  //Infectivity ~ Infectv ~ ii
         std::vector<double> MG;   //gams, Mature Gametocytes ~ MGG ~ ggmd
         std::vector<double> G;    //gametocyte?
         std::vector<double> MSP;  // merozoite proteins?
         std::vector<double> cyt;  //Innate immunity? //cytokine-mediated? 
         std::vector<double> AB;   // Y[antibody] // Antibody?

		std::vector<int> scheduled_drug_timestamps; //keep track of each time human is scheduled for drug (each time drug is scheduled)
        std::vector<int> recieved_drug_timestamps; //keep track of each time human recived drug (each time drug is recieved)

		HumanAgent(int h_id, int ts, int level, int npk, int npip) 
        {
            id= h_id;
            total_steps = ts;


            AI.resize(total_steps, 0.0);
            P0.resize(total_steps, 0.0);
            infect.resize(total_steps, 0.0);
            MG.resize(total_steps, 0.0);
            G.resize(total_steps, 0.0);
            MSP.resize(total_steps, 0.0);
            cyt.resize(total_steps, 0.0);
            AB.resize(total_steps, 0.0);


			vpk.resize(npk, total_steps+1); //DHA/PK
            vvpk.resize(level+1, total_steps+1); //[DHA/PK]

            vpip.resize(npip, total_steps+1); //PIP/PD
            vvpip.resize(level+1, total_steps+1); //[Pip/PD] 

            //vpk=zeros(npk,asteps+1);
            vpk.setZero();
            vvpk.setZero();

            //vpip=zeros(npip,asteps+1);
            vpip.setZero();

            //vvpip=zeros(levels,asteps+1);
            vvpip.setZero();  

		}

        int getID(){
            return id;
        }

		void print_human_WH_param_value(vector<double> vec) const{

           std::cout << "Human agent ID= " << id <<"\n";

		   std::cout << std::fixed;
			std::cout << std::setprecision(4);
			for (int i = 0; i < total_steps; i++) {
				std::cout << vec[i] << " ";
			}
			std::cout << " \n ";
		}

		void print_human_WH_params_values() const{

			//{"AI", "P0", "infect", "MG", "G", "MSP", "cyt", "AB"};
		    std::cout << std::fixed;
			std::cout << std::setprecision(4);
		
            std::cout << "Human agent ID= " << id <<"\n";

			std::cout << "AI= \n";
			for (int i = 0; i < total_steps; i++) 
				std::cout << AI[i] << " ";
			std::cout << " \n ";

			std::cout << "P0= \n";
			for (int i = 0; i < total_steps; i++) 
				std::cout << P0[i] << " ";
			std::cout << " \n ";

			std::cout << "infect= \n";
			for (int i = 0; i < total_steps; i++) 
				std::cout << infect[i] << " ";
			std::cout << " \n ";

			std::cout << "MG= \n";
			for (int i = 0; i < total_steps; i++) 
				std::cout << MG[i] << " ";
			std::cout << " \n ";

			std::cout << "G= \n";
			for (int i = 0; i < total_steps; i++) 
				std::cout << G[i] << " ";
			std::cout << " \n ";

			std::cout << "MSP= \n";
			for (int i = 0; i < total_steps; i++) 
				std::cout << MSP[i] << " ";
			std::cout << " \n ";

			std::cout << "cyt= \n";
			for (int i = 0; i < total_steps; i++) 
				std::cout << cyt[i] << " ";
			std::cout << " \n ";

			std::cout << "AB= \n";
			for (int i = 0; i < total_steps; i++) 
				std::cout << AB[i] << " ";
			std::cout << " \n ";			

            std::cout << "--------- \n ";			
		}
	};

    //================================
	 class Human_Manager {

		public: 
		   const int sum_total_population;
           const int num_time_steps = 0;
		   vector<HumanAgent> list_of_human_agents;

		
        private:
		   const int whd_level;
           const int whd_npk;
           const int whd_npip;

		public:

			Human_Manager(int initial_pop, int time_step, int level, int npk, int npip): 
			                sum_total_population(initial_pop), num_time_steps(time_step), 
			                whd_level(level), whd_npk(npk), whd_npip(npip) 
			{

				for (int i=1; i<=sum_total_population; i++) {
				//for (int i=0; i<sum_total_population; i++) {

					list_of_human_agents.push_back(HumanAgent(i, num_time_steps, whd_level, whd_npk, whd_npip));
				}
			}

            ~Human_Manager() {};
	 };

}

#endif