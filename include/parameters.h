#ifndef PARAMETERS_H
#define PARAMETERS_H


#include "classes.h"
#include "utilities.h"


// std::map <std::string, std::string> BASIC_ATOM_TYPES;
// std::map <std::string, std::string> mass_library;
// std::vector <std::string> mass_library_keys;
// std::map <std::string, std::string> bond_library;
// std::vector <std::string> bond_library_keys;
// std::map <std::string, std::string> angle_library;
// std::vector <std::string> angle_library_keys;
// std::map <std::string, std::string> dihedral_library;
// std::vector <std::string> dihedral_library_keys;
// std::map <std::string, std::string> vdw_library;
// std::vector <std::string> vdw_library_keys;


std::string FindMissingMass(std::string atom_name);
std::string FindMissingBond(std::string a1, std::string a2);
std::string FindMissingAngle(std::string a1, std::string a2, std::string a3);
std::string FindMissingDihedral(std::string a1, std::string a2, std::string a3, std::string a4);
std::string FindMissingTorsion(std::string a1, std::string a2, std::string a3, std::string a4);
std::string FindMissingNonbonded(std::string a1);

#endif