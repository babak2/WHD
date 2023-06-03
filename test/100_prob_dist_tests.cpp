#include <iostream>
 
#include <vector>

#include "third_party/catch2/catch.hpp"

#include <random>

#include <iomanip>

#include "util/CSVWriter.h"


using namespace std;

TEST_CASE( "100.1: Prob dist tests", "[prob_dist:tests]" ) {

    std::random_device rd;

    // Engines 
    std::mt19937 r_engin(rd());

    // Gaussian/ Normal Distribtuions
    std::normal_distribution<double> randn_1(0.0, 1.0);  //uniformly distributed on the interval (0, 1)
	std::normal_distribution<double> randn_2(0, 1); //normally distributed having zero mean and variance one.
	std::normal_distribution<double> randn_3{0.0,1.0}; //normally distributed having zero mean and variance one.
    std::normal_distribution<double> randn{0,1}; //normally distributed having zero mean and variance one.


    // Uniform Distribtuion
    //std::uniform_real_distribution<double> rand(0.0,1.0); //uniformly distributed on the interval (0, 1)
	//std::uniform_real_distribution<double> rand(0,1); //uniformly distributed on the interval (0, 1)
	//std::uniform_real_distribution<double> rand{0.0,1.0}; //uniformly distributed on the interval (0, 1)
	std::uniform_real_distribution<double> rand{0,1}; //uniformly distributed on the interval (0, 1)

	//std::uniform_int_distribution<double> unif(0,1);
	//std::uniform_real_distribution<double> unif(0.0,1.0);

	std::minstd_rand stdrand;

	int nb_of_draws = 1000;
	std::cout << "generate random numbers from prob dists. Nb of draws: " << nb_of_draws <<"\n";
	//int count[10] = {0};

    {
		util::CSVWriter  randn_csv_writer("randn_test_cpp.csv");
		std::vector<std::string> randn_vals;

		std::cout << "draw from randn (normal dist), write to: randn_test_cpp.csv \n";
		for (int i = 0; i < nb_of_draws; i++)  {

			randn_vals.push_back(to_string(randn(r_engin)));
			randn_csv_writer.addRow(randn_vals);
			randn_vals.clear();

			//std::cout << randn(r_engin) << " ";
			//++count[( randn(r_engin) % 100) / 10];
			//++count[( randn(r_engin) % 100) / 10];
			//++count[randn(r_engin) * 100)];
		}
		//std::cout << " \n ";
    }


    {
		util::CSVWriter  rand_csv_writer("rand_test_cpp.csv");
		std::vector<std::string> rand_vals;

		std::cout << "draw from rand (uniform dist), write to: rand_test_cpp.csv \n";
		for (int i = 0; i < nb_of_draws; i++)  {

			rand_vals.push_back(to_string(rand(r_engin)));
			rand_csv_writer.addRow(rand_vals);
			rand_vals.clear();
		}
		//std::cout << " \n ";
    }

	{
        //P0(j,1)=max(round(normrnd(P0(j,1),20000) ),5*10^2);
        double p0 = 4*pow(10,4);
		//int p0 = 4*pow(10,4);
		util::CSVWriter  normrnd_csv_writer("normrnd_test_cpp.csv");
		std::vector<std::string> normrnd_vals;
		std::normal_distribution<> normrnd(p0, 20000.0);  
            
		std::cout << "draw from normrnd (normal dist), write to: normrnd_test_cpp.csv \n";
		for (int i = 0; i < nb_of_draws; i++)  {

			normrnd_vals.push_back(to_string(normrnd(r_engin)));
			normrnd_csv_writer.addRow(normrnd_vals);
			normrnd_vals.clear();
		}
		//std::cout << " \n ";

	}


//==============================

	std::cout << " PASS\n";

}