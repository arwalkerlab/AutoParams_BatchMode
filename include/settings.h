#ifndef SETTINGS_H
#define SETTINGS_H

#include "utilities.h"

class Settings
{
    public:
        // Logfiles and Job Working Directory
        std::string outfile;
        std::string errfile;
        std::string job_dir;
        std::string mol2file;
        std::string frcmodfile;

        // Program Flags
        bool TC_EXISTS;
        bool PSI4_EXISTS;
        bool AIMNET_EXISTS;
        bool TLEAP_EXISTS;
        bool PARMED_EXISTS;
        bool CPPTRAJ_EXISTS;
        bool ANTECHAMBER_EXISTS;

        // Job Flags
        bool OPTIMIZE_FIRST;
        bool INCLUDE_TIMESTAMPS;
        
        // Molecule Information;
        int mol_charge;
        int mol_spin;
        std::string forcefield;
        std::string inputfile;
        std::string head_atom_name;
        std::string tail_atom_name;
        std::vector <std::string> dummy_atom_names;

        // Functions
        void CheckPrograms();
        void SetUpLogFiles();
        void parse_command_line(int argc,char **argv);
        void Output(std::string output);
        void Error(std::string error);

        // Constructor/Destructor
        Settings();
        ~Settings();
};

#endif