/*********************************************************
 * File: malaria_model.cpp
 * Description: This is the main file for malaria_model (MATLAB conversion)
 * 
 * Author: Babak Mahdavi Ardestani
 * Email: babak.m.ardestani@gmail.com
 * Created:  2021-10-15
 * 
 *********************************************************/
#define DO_QUOTE(X) #X
#define QUOTE(X) DO_QUOTE(X)

#include <iomanip>
#include "common/wh_dynamics.hpp"
#include "human/human.h"
#include "util/util.h"

#include <experimental/filesystem> // -lstdc++fs in LDLIBS
namespace fs = std::experimental::filesystem;

#include "rapidjson/document.h"



bool validate_config(const std::string &configuration_file_name)
{

   bool check_pass = false;
   rapidjson::Document config_doc = util::get_json_from_file(configuration_file_name);

   // JSON Validation
   assert(config_doc["simulation"]["schema_file_name"].IsString());
   rapidjson::Document rj_doc_config_schema =
       util::get_json_from_file(config_doc["simulation"]["schema_file_name"].GetString());

   assert(util::validate_json_against_schema(&rj_doc_config_schema, &config_doc));

   assert(config_doc["simulation"]["total_steps"].IsInt());

   check_pass = true;
   return check_pass;
}

//--------  MAIN() --------------
int main(int argc, char **argv)
{

   for (int i = 0; i < argc; ++i)
   {
      std::cout << argv[i] << " ";
   }
   std::cout << "\n";

   util::InputParser input_parser(argc, argv);

   const std::string &config_file_name = input_parser.getCmdOption("-c");

   std::cout << "config_file_name: " << config_file_name << "\n";

   if (!config_file_name.empty())
   {

      if (!fs::exists(config_file_name.c_str()))
      {
         std::cout << "File (" << config_file_name << ") does not exist." << std::endl;
         return EXIT_FAILURE;
      }

      std::cout << "Reading configuration from ";
      std::cout << config_file_name << "" << std::endl;
   }
   else
   {

      std::cout << "Please give name of configuration file (JSON) using the \"-c [file]\" commandline option.";
      std::cout << " (e.g. -c data/config.json)" << std::endl;
      return EXIT_FAILURE;
   }

   const std::string &output_folder_name = input_parser.getCmdOption("-o");
   std::cout << "output_folder_name: " << output_folder_name << "\n";

   if (!output_folder_name.empty())
   {

      if (!fs::exists(output_folder_name.c_str()))
      {
         if (!fs::create_directory(output_folder_name.c_str()))
         {
            std::cout << "Could not create directory at (" << output_folder_name << ")." << std::endl;
            return EXIT_FAILURE;
         }
      }
   }
   else
   {
      // TODO: give a default value and proceed

      std::cout << "Please give name of output directory using the \"-o [path]\" commandline option.";
      std::cout << " (e.g. -o outputs)" << std::endl;

      return EXIT_FAILURE;
   }

   // -i [prefix], optional, if not given a timestamp will be used to prefix this run
   const std::string &output_prefix_name = input_parser.getCmdOption("-i");
   std::cout << "output_prefix_name: " << output_prefix_name << "\n";

   if (!validate_config(config_file_name))
   {
      std::cout << "Error in configuration file (" << config_file_name << ")\n";
      return EXIT_FAILURE;
   }
   else
   {

      std::cout << "Schema check passed\n";
      // return EXIT_SUCCESS;
   }

   // const int total_steps  =  atoi(argv[1]);  //simulation length
   // const int nb_of_agents =  atoi(argv[2]);

   const rapidjson::Document config = util::get_json_from_file(config_file_name);
   int total_steps = config["simulation"]["total_steps"].GetInt();
   int nb_of_agents = config["simulation"]["total_population"].GetInt();

   //------------

   human::Human_Manager human_mgr(nb_of_agents, total_steps,
                                   config["whd"]["whd_level"].GetInt(),
                                   config["whd"]["npk"].GetInt(),
                                   config["whd"]["npip"].GetInt());

   common::WHD_Manager whd_mgr(config_file_name,
                               output_folder_name,
                               output_prefix_name,
                               nb_of_agents, &human_mgr);

   for (int t = 0; t < total_steps; t++)
   {
      whd_mgr.step_within_host_dynamics(t);
   }

   std::cout << " -------------------\n";

   // human_mgr.print_WH_vals_for_all_humans();

   //whd_mgr.print_P0_patent_for_all_humans();

   //whd_mgr.print_P0_duration_for_all_humans();

   whd_mgr.export_to_csv_whd_param_vals_for_all_humans(output_folder_name + util::kPathSeparator + output_prefix_name);
   //whd_mgr.export_to_csv_whd_param_avg_vals_of_all_humans(output_folder_name + util::kPathSeparator + output_prefix_name);
   whd_mgr.export_to_csv_P0_vals_for_all_humans(output_folder_name + util::kPathSeparator + output_prefix_name);
   //whd_mgr.export_to_csv_P0_patent_for_all_humans(output_folder_name + util::kPathSeparator + output_prefix_name);
   //whd_mgr.export_to_csv_P0_peakpar_for_all_humans(output_folder_name + util::kPathSeparator + output_prefix_name);
   
   whd_mgr.export_to_csv_P0_duration_for_all_humans(output_folder_name + util::kPathSeparator + output_prefix_name);

   whd_mgr.export_to_csv_recieved_drug_timestamps_for_all_humans(output_folder_name + util::kPathSeparator + output_prefix_name);
   whd_mgr.export_to_csv_scheduled_drug_timestamps_for_all_humans(output_folder_name + util::kPathSeparator + output_prefix_name);
   
   whd_mgr.export_to_csv_agents_info_missed_drug(output_folder_name + util::kPathSeparator + output_prefix_name);

   
   //-----to be removed?:

   //whd_mgr.export_to_csv_total_MDA_vals(output_folder_name + util::kPathSeparator + output_prefix_name);
   //whd_mgr.export_to_cvs_P0_patent_for_no_drug_receivers(output_folder_name + util::kPathSeparator + output_prefix_name );

   std::cout << " ===================\n";

   return 0;
}
