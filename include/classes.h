#ifndef CLASSES_H
#define CLASSES_H

#include "utilities.h"
#include "settings.h"

class Atom
{
    public:
        std::string atom_name;
        int atom_number;
        std::string residue_name;
        int residue_number;
        double xx;
        double yy;
        double zz;
        std::string element;
        int formal_charge;
        double resp_charge;
        double mass;
        int electron_counter;
        double vdw_radius;
        std::string atom_type;
        std::vector <int> bonded_to_indexes;
        std::vector <std::string> bonded_to_elements;
        
        Atom(std::string line);
        ~Atom();
        void set_resname(std::string resname);
        void set_atomname(std::string atomname);
        void set_resp_charge(double respcharge);
        std::string print_PDB_line();
        std::string print_mol2_line();
};


class Molecule
{
public:
    int n_atoms;
    int charge;
    int spin_mult;
    int total_electrons;
    std::string head_atom_name;
    std::string tail_atom_name;
    std::string res_name;
    std::vector <Atom> atoms;
    std::vector <std::vector<int>> bonds;
    std::vector <std::vector<int>> angles;
    std::vector <std::vector<int>> dihedrals;
    std::vector <std::vector<int>> torsions;
    std::vector <std::vector<int>> four_rings;
    std::vector <std::vector<int>> five_rings;
    std::vector <std::vector<int>> six_rings;
    
    // Missing Parameters
    std::vector<std::string> missing_masses;
    std::vector<std::string> missing_bonds;
    std::vector<std::string> missing_angles;
    std::vector<std::string> missing_torsions;
    std::vector<std::string> missing_dihedrals;
    std::vector<std::string> missing_vdw_params;
    
    Molecule(Settings settings);
    ~Molecule();
    void Write_PDB(std::string job_dir, std::string inputfile);
    void move_to_COM();
    bool SpinChargeValidate();
    void FindBonds();
    void FindAngles();
    void FindDihedrals();
    void FindTorsions();
    void FindRings();
};


class Frcmod_File
{
    private:
        std::vector <std::string> params_added;
        std::string masses_section;
        std::string bonds_section;
        std::string angles_section;
        std::string torsions_section;
        std::string dihedrals_section;
        std::string vdw_section;

    public:
        std::string frcmod_filename;
        std::ofstream frcmod_write;
        void Add_Mass(std::string line);
        void Add_Bond(std::string line);
        void Add_Angle(std::string line);
        void Add_Dihedral(std::string line);
        void Add_Torsion(std::string line);
        void Add_vdW(std::string line);
        Frcmod_File(std::string filename);
        ~Frcmod_File();
        void Write_Frcmod_File(Settings settings);

};

#endif