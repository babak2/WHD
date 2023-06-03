#include "third_party/catch2/catch.hpp"
#include <iomanip>
#include "common/wh_dynamics.hpp"
#include "util/CSVWriter.h"
#include "Eigen/Dense"
#include "Eigen-unsupported/Eigen/MatrixFunctions" 

using namespace std;

/*
This test uses WH_Dynamics class directly (and not through HumanManager) for simulating 
within host parameters for a single (human) agent for a period of 20 time steps. 
*/


TEST_CASE( "001.1: Within Host Dynamic Test", "[whd:test]" ) {

   class Agent 
   {

	private: 
		int total_steps = 0;  //simulation time
		typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_dyn_t;
        //typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::DontAlign> vector_dyn_t; //vertical vector

    public:

	     int id;
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
		std::vector<double> AB;   // Y[antibody] // Antibody ?

		Agent(int ts) {
			
			total_steps = ts;
			
			AI.resize(total_steps, 0.0);
			P0.resize(total_steps, 0.0);
			infect.resize(total_steps, 0.0);
			MG.resize(total_steps, 0.0);
			G.resize(total_steps, 0.0);
			MSP.resize(total_steps, 0.0);
			cyt.resize(total_steps, 0.0);
			AB.resize(total_steps, 0.0);

		}

		void print_array(double* arr) {

		    std::cout << std::fixed;
			std::cout << std::setprecision(4);
			for (int i = 0; i < total_steps; i++) {
				std::cout << arr[i] << " ";
			}
			std::cout << " \n ";
		}

		void print_all_arrays() {

		    std::cout << std::fixed;
			std::cout << std::setprecision(4);
		
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
			
		}

	};
	
	//------------------------------------

	int total_steps = 20;  //simulation length

	Agent agent1 =  Agent(total_steps);
	
	common::WH_Dynamics  whd = common::WH_Dynamics
	                        (
		                        total_steps, //timepoints
	                            0, //level
								2000000000, //II50, 
								11.0, //pmf, 
								0.75, // gamma,
                                0.225, // maxcyt, 
								1.175, //kfever, 
								0.1, // kfever_multiplier
								0.35, //memory, 
								150, //decay, 
								0.03, //switchrate,
                                2.0, //immunestack, 
								22.0, //switcheffect, 
								50000, // maxinit, 
								5.0, //etac, 
								3.0, //etah, 
                                50000000000, //x50ab, 
								0.5,    //x50ab_multiplier
								0.115, //kab, 
								0.02, // kab_multiplier
								1.0, //abprod, 
								20.0, // abdecay, 
								3.15, // killcyt, 
                                2.15, //killab, 
								0.05, //killab_multiplier
								0.0853, // kmsp,
								0.358, // Cmer, 
								9,    //kg_divisor
                                10,   //kmg_divisor
					            0.01, // kgam,
                                0.0025, //kgametocyte, 
								4000000000, //pt, 
								0.005, //C50g, 
								6000000, //gam50, 
								4000.0, //dosingpk,
                                960.0, //dosingpip, 
								3.25, //killratedha, 
								1.75, //killratepip, 
								300.0, //ce50dha, 
								450.0, //ce50pip,
                                3.0, //hdha, 
								3.0, //hpip, 
								9, //npk, 
								6, //npip, 
								4000, //p0_initial_val
                                true, //with_drug
								3.2500, //CL,  
								5.3750, //V2, 
								0.982, //THETA_3_divident
                                24,    //THETA_3_divisor
                                7.0, //NN, 
								0.75, //ACL2
                                54, // M_WE2
                                25, //AGE2
                                54, //WT2
                                1329.6, //THETA2_1
                                0.575, //EM502
                                5.51, //HILL2
								7440.0, //Q12, 
								2520.0, //Q22, 
								69840.0, //V22, 
								117840.0, //V32, 
								741600.0,//V42,
							    2.11, //THETA2_7_divident
                                24,  //THETA2_7_divisor
								2.0 //NN2
	                        );

	cout << "simul time/ total steps: "<< total_steps<< "\n"; 

	for (int t=0; t<total_steps; t++) {

		whd.step(t, agent1.tt, agent1.mda1, agent1.mda2, agent1.mda3, 
                     agent1.kfeverj, agent1.x50abj, agent1.kabj, agent1.killabj,
                     agent1.vpk, agent1.vvpk, agent1.vpip, agent1.vvpip,
                     agent1.P0, agent1.AI, agent1.infect, agent1.MG,
                     agent1.G, agent1.MSP, agent1.cyt, agent1.AB);
	}

	cout << " Agent1:\n"; 
    agent1.print_all_arrays();

	cout << "---------------------\n"; 

	std::cout << " PASS\n";

}


//==========================================================

/*
This test uses WH_Dynamics class directly (and not through HumanManager) for simulating 
within host parameters for a 10 (human) agents for a period of 360 time steps. 
It also calculate the averages of these simulations among 10 agents and write the outcome in a CSV file.
*/

TEST_CASE( "001.2: Within Host Dynamic Test with Multiple Agents", "[whd_multiple:test]" ) {

    cout << "Test Within Host Dynamics with 10 agents:\n"; 


    class Agent {

		private: 
			int total_steps = 0;  //simulation time
			typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_dyn_t;
			//typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::DontAlign> vector_dyn_t; //vertical vector

        public:

			int id;
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
			std::vector<double> AB;   // Y[antibody] // Antibody ?

			Agent(int ts) {
				total_steps = ts;

				AI.resize(total_steps, 0.0);
				P0.resize(total_steps, 0.0);
				infect.resize(total_steps, 0.0);
				MG.resize(total_steps, 0.0);
				G.resize(total_steps, 0.0);
				MSP.resize(total_steps, 0.0);
				cyt.resize(total_steps, 0.0);
				AB.resize(total_steps, 0.0);				
			}

			void print_array(double* arr) {

				std::cout << std::fixed;
				std::cout << std::setprecision(4);
				for (int i = 0; i < total_steps; i++) {
					std::cout << arr[i] << " ";
				}
				std::cout << " \n ";
			}

			void calculate_averages(int num_agents, Agent agent1, Agent agent2,  Agent agent3,  Agent agent4,  Agent agent5,
									Agent agent6,  Agent agent7,  Agent agent8,  Agent agent9,  Agent agent10) {
				
				//std::string filename = "all_avg_values.csv";
				util::CSVWriter  csv_writer("all_avg_values.csv"); //write all together
				//util::CSVWriter  csv_writer(filename); //write all together
				
				/*****if you want to write them separatly then: 
				- comment-out the parameter value you want to record and comment-in 'addHeader' line below.
				- Also comment-in the other parameter values (except the desired one) in 'avg_vals.push_back' instrunctions further below */
			
				//util::CSVWriter  csv_writer("AI_avg_cpp.csv");
				//util::CSVWriter  csv_writer("P0_avg_cpp.csv");
				//util::CSVWriter  csv_writer("infect_avg_cpp.csv");
				//util::CSVWriter  csv_writer("MG_avg_cpp.csv");
				//util::CSVWriter  csv_writer("G_avg_cpp.csv");
				//util::CSVWriter  csv_writer("MSP_avg_cpp.csv");
				//util::CSVWriter  csv_writer("cyt_avg_cpp.csv");
				//util::CSVWriter  csv_writer("AB_avg_cpp.csv");

				//std::vector<std::string> header = {"AI", "P0", "infect", "MG", "G", "MSP", "cyt", "AB"};
				std::vector<std::string> header = {"AI_avg", "P0_avg", "infect_avg", "MG_avg", "G_avg", "MSP_avg", "cyt_avg", "AB_avg"};

				csv_writer.addHeader(header);

				std::vector<std::string> avg_vals;

				double* AI_avg = new double[total_steps];
				double* P0_avg = new double[total_steps];
				double* infect_avg = new double[total_steps];
				double* MG_avg = new double[total_steps];
				double* G_avg = new double[total_steps];
				double* MSP_avg = new double[total_steps];
				double* cyt_avg = new double[total_steps];
				double* AB_avg = new double[total_steps];


				for (int i = 0; i < total_steps; i++) {

					AI_avg[i] = (agent1.AI[i]+ agent2.AI[i]+ agent3.AI[i]+ agent4.AI[i]+ agent5.AI[i]+ agent6.AI[i]
								+ agent7.AI[i]+ agent8.AI[i]+ agent9.AI[i]+ agent10.AI[i])/num_agents;

					P0_avg[i] = (agent1.P0[i]+ agent2.P0[i]+ agent3.P0[i]+ agent4.P0[i]+ agent5.P0[i]+ agent6.P0[i]
								+ agent7.P0[i]+ agent8.P0[i]+ agent9.P0[i]+ agent10.P0[i])/num_agents;

					infect_avg[i] = (agent1.infect[i]+ agent2.infect[i]+ agent3.infect[i]+ agent4.infect[i]+ agent5.infect[i]+ agent6.infect[i]
								+ agent7.infect[i]+ agent8.infect[i]+ agent9.infect[i]+ agent10.infect[i])/num_agents;

					MG_avg[i] = (agent1.MG[i]+ agent2.MG[i]+ agent3.MG[i]+ agent4.MG[i]+ agent5.MG[i]+ agent6.MG[i]
								+ agent7.MG[i]+ agent8.MG[i]+ agent9.MG[i]+ agent10.MG[i])/num_agents;

					G_avg[i] = (agent1.G[i]+ agent2.G[i]+ agent3.G[i]+ agent4.G[i]+ agent5.G[i]+ agent6.G[i]
								+ agent7.G[i]+ agent8.G[i]+ agent9.G[i]+ agent10.G[i])/num_agents;

					MSP_avg[i] = (agent1.MSP[i]+ agent2.MSP[i]+ agent3.MSP[i]+ agent4.MSP[i]+ agent5.MSP[i]+ agent6.MSP[i]
								+ agent7.MSP[i]+ agent8.MSP[i]+ agent9.MSP[i]+ agent10.MSP[i])/num_agents;

					cyt_avg[i] = (agent1.cyt[i]+ agent2.cyt[i]+ agent3.cyt[i]+ agent4.cyt[i]+ agent5.cyt[i]+ agent6.cyt[i]
								+ agent7.cyt[i]+ agent8.cyt[i]+ agent9.cyt[i]+ agent10.cyt[i])/num_agents;

					AB_avg[i] = (agent1.AB[i]+ agent2.AB[i]+ agent3.AB[i]+ agent4.AB[i]+ agent5.AB[i]+ agent6.AB[i]
								+ agent7.AB[i]+ agent8.AB[i]+ agent9.AB[i]+ agent10.AB[i])/num_agents;


				//uncomment/comment the value you want to record. Do same for csv_writer above
				
					avg_vals.push_back(to_string(AI_avg[i]));
					avg_vals.push_back(to_string(P0_avg[i]));
					avg_vals.push_back(to_string(infect_avg[i]));
					avg_vals.push_back(to_string(MG_avg[i]));
					avg_vals.push_back(to_string(G_avg[i]));
					avg_vals.push_back(to_string(MSP_avg[i]));
					avg_vals.push_back(to_string(cyt_avg[i]));
					avg_vals.push_back(to_string(AB_avg[i]));

					csv_writer.addRow(avg_vals);

					avg_vals.clear();
				}

				std::cout << std::fixed;
				std::cout << std::setprecision(4);
			
				std::cout << "AI_avg= \n";
				for (int i = 0; i < total_steps; i++) 
					std::cout << AI_avg[i] << " ";
				std::cout << " \n ";

				std::cout << "P0_avg= \n";
				for (int i = 0; i < total_steps; i++) 
					std::cout << P0_avg[i] << " ";
				std::cout << " \n ";

				std::cout << "infect_avg= \n";
				for (int i = 0; i < total_steps; i++) 
					std::cout << infect_avg[i] << " ";
				std::cout << " \n ";

				std::cout << "MG_avg= \n";
				for (int i = 0; i < total_steps; i++) 
					std::cout << MG_avg[i] << " ";
				std::cout << " \n ";

				std::cout << "G_avg= \n";
				for (int i = 0; i < total_steps; i++) 
					std::cout << G_avg[i] << " ";
				std::cout << " \n ";

				std::cout << "MSP_avg= \n";
				for (int i = 0; i < total_steps; i++) 
					std::cout << MSP_avg[i] << " ";
				std::cout << " \n ";

				std::cout << "cyt_avg= \n";
				for (int i = 0; i < total_steps; i++) 
					std::cout << cyt_avg[i] << " ";
				std::cout << " \n ";

				std::cout << "AB_avg= \n";
				for (int i = 0; i < total_steps; i++) 
					std::cout << AB_avg[i] << " ";
				std::cout << " \n ";

			}

			void print_all_arrays() {

				std::cout << std::fixed;
				std::cout << std::setprecision(4);
			
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

				std::cout << "cyt= \n";
				for (int i = 0; i < total_steps; i++) 
					std::cout << cyt[i] << " ";
				std::cout << " \n ";

				std::cout << "MSP= \n";
				for (int i = 0; i < total_steps; i++) 
					std::cout << MSP[i] << " ";
				std::cout << " \n ";

				std::cout << "G= \n";
				for (int i = 0; i < total_steps; i++) 
					std::cout << G[i] << " ";
				std::cout << " \n ";

				std::cout << "AB= \n";
				for (int i = 0; i < total_steps; i++) 
					std::cout << AB[i] << " ";
				std::cout << " \n ";			
			}
		};
		
		//---------------------------------

		int total_steps = 360;  //simulation length

		Agent agent1 =  Agent(total_steps);
		Agent agent2 =  Agent(total_steps);
		Agent agent3 =  Agent(total_steps);
		Agent agent4 =  Agent(total_steps);
		Agent agent5 =  Agent(total_steps);
		Agent agent6 =  Agent(total_steps);
		Agent agent7 =  Agent(total_steps);
		Agent agent8 =  Agent(total_steps);
		Agent agent9 =  Agent(total_steps);
		Agent agent10 =  Agent(total_steps);

	    common::WH_Dynamics  whd = common::WH_Dynamics
	                        (
		                        total_steps, //timepoints
	                            0, //level
								2000000000, //II50, 
								11.0, //pmf, 
								0.75, // gamma,
                                0.225, // maxcyt, 
								1.175, //kfever, 
								0.1, // kfever_multiplier
								0.35, //memory, 
								150, //decay, 
								0.03, //switchrate,
                                2.0, //immunestack, 
								22.0, //switcheffect, 
								50000, // maxinit, 
								5.0, //etac, 
								3.0, //etah, 
                                50000000000, //x50ab, 
								0.5,    //x50ab_multiplier
								0.115, //kab, 
								0.02, // kab_multiplier
								1.0, //abprod, 
								20.0, // abdecay, 
								3.15, // killcyt, 
                                2.15, //killab, 
								0.05, //killab_multiplier
								0.0853, // kmsp,
								0.358, // Cmer, 
								9,    //kg_divisor
                                10,   //kmg_divisor
					            0.01, // kgam,
                                0.0025, //kgametocyte, 
								4000000000, //pt, 
								0.005, //C50g, 
								6000000, //gam50, 
								4000.0, //dosingpk,
                                960.0, //dosingpip, 
								3.25, //killratedha, 
								1.75, //killratepip, 
								300.0, //ce50dha, 
								450.0, //ce50pip,
                                3.0, //hdha, 
								3.0, //hpip, 
								9, //npk, 
								6, //npip, 
								4000, //p0_initial_val
                                true, //with_drug
								3.2500, //CL,  
								5.3750, //V2, 
								0.982, //THETA_3_divident
                                24,    //THETA_3_divisor
                                7.0, //NN, 
								0.75, //ACL2
                                54, // M_WE2
                                25, //AGE2
                                54, //WT2
                                1329.6, //THETA2_1
                                0.575, //EM502
                                5.51, //HILL2
								7440.0, //Q12, 
								2520.0, //Q22, 
								69840.0, //V22, 
								117840.0, //V32, 
								741600.0,//V42,
							    2.11, //THETA2_7_divident
                                24,  //THETA2_7_divisor
								2.0 //NN2
	                        );


		for (int t=0; t<total_steps; t++) 
		{
			
		    whd.step(t, agent1.tt, agent1.mda1, agent1.mda2, agent1.mda3, 
                     agent1.kfeverj, agent1.x50abj, agent1.kabj, agent1.killabj,
                     agent1.vpk, agent1.vvpk, agent1.vpip, agent1.vvpip,
                     agent1.P0, agent1.AI, agent1.infect, agent1.MG,
                     agent1.G, agent1.MSP, agent1.cyt, agent1.AB);

			/*whd.step(t, agent2.tt, agent2.mda1, agent2.mda2, agent2.mda3, 
                     agent2.kfeverj, agent2.x50abj, agent2.kabj, agent2.killabj,
                     agent2.vpk, agent2.vvpk, agent2.vpip, agent2.vvpip,
                     agent2.P0, agent2.AI, agent2.infect, agent2.MG,
                     agent2.G, agent2.MSP, agent2.cyt, agent2.AB); */
		
																						
		} 
		
		
		cout << "---Agent1:\n"; 
		agent1.print_all_arrays();

		cout << "---Agent2:\n"; 
		agent2.print_all_arrays();

		cout << "---Agent3:\n"; 
		agent3.print_all_arrays();

		cout << "---Agent4:\n"; 
		agent4.print_all_arrays();

		cout << "---Agent5:\n"; 
		agent5.print_all_arrays();

		cout << "---Agent6:\n"; 
		agent6.print_all_arrays();

		cout << "---Agent7:\n"; 
		agent7.print_all_arrays();

		cout << "---Agent8:\n"; 
		agent8.print_all_arrays();

		cout << "---Agent9:\n"; 
		agent9.print_all_arrays();

		cout << "---Agent10:\n"; 
		agent10.print_all_arrays();


		cout << " +++++++ AVG ++++++++++++++\n"; 
		agent1.calculate_averages(10, agent1, agent2, agent3, agent4, agent5, agent6, agent7, agent8, agent9, agent10);
		std::cout << "Writing avg values in csv file "<<"\n ";

    //------------------------------

	std::cout << " PASS\n";
}
