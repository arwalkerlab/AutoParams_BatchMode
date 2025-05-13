#ifndef AUTOPARAMS_H
#define AUTOPARAMS_H

#include "utilities.h"
#include "classes.h"
#include "settings.h"

void write_TC_resp_input(Settings settings, Molecule mol);
bool validate_TC_resp_output(Settings settings);
void parse_TC_resp_output(Molecule &mol,Settings settings);
void Generate_Mol2_File(Settings settings);
void Validate_Mol2_File(Settings settings);
void Add_Charges_To_Mol2(Settings settings, Molecule &mol);
void Check_For_Missing_Parameters(Settings settings, Molecule &mol);
void Build_Frcmod(Frcmod_File &frcmod, Molecule mol, Settings settings);

#endif