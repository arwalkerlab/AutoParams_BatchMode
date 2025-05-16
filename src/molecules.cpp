#include "classes.h"

Molecule::Molecule(Settings settings)
{
    std::string filename = settings.inputfile;
    int formal_charge = settings.mol_charge;
    int spin = settings.mol_spin;
    head_atom_name = settings.head_atom_name;
    tail_atom_name = settings.tail_atom_name;
    
    std::vector<std::string> atom_name_set;
    std::vector<std::string> residue_name_set;
    
    charge = formal_charge;
    spin_mult = spin;
    res_name = "UNK";
    std::ifstream pdbfile(filename);
    std::string line;
    
    // Parse PDB lines
    while (getline(pdbfile,line))
    {
        if ((line.find("ATOM") == std::string::npos) && (line.find("HETATM") == std::string::npos))
        {
            continue;
        }
        Atom curr_atom(line);
        // if (std::find(settings.dummy_atom_names.begin(),settings.dummy_atom_names.end(),curr_atom.atom_name) != settings.dummy_atom_names.end())
        // {
        //     continue;
        // }
        int name_counter=1;
        if (res_name == "UNK")
        {
            res_name = curr_atom.residue_name;
        }
        while (std::find(atom_name_set.begin(),atom_name_set.end(),curr_atom.atom_name) != atom_name_set.end())
        {
            std::stringstream buffer;
            std::string newname;
            buffer.str("");
            buffer << curr_atom.element << name_counter;
            name_counter++;
            newname = trim_whitespace(buffer.str());
            curr_atom.set_atomname(newname);
        }
        curr_atom.residue_number = 1;
        curr_atom.residue_name = res_name;
        atom_name_set.push_back(curr_atom.atom_name);
        residue_name_set.push_back(curr_atom.residue_name);
        atoms.push_back(curr_atom);
    }

    n_atoms = atoms.size();

    // check unique atom names
    std::set<std::string> s(residue_name_set.begin(), residue_name_set.end());
    if (s.size() != 1)
    {
        for (Atom atom : atoms)
        {
            atom.residue_name = "LIG";
        }
    }
}

Molecule::~Molecule()
{
}

void Molecule::Write_PDB(std::string job_dir,std::string inputfile)
{
    std::ofstream outfile;
    outfile.open(job_dir+"/"+inputfile,std::ios::out);
    for (Atom atom : atoms)
    {
        outfile << atom.print_PDB_line();
    }
    outfile << "END\n";
    outfile.close();
}

void Molecule::move_to_COM()
{
    double com_x = 0.0;
    double com_y = 0.0;
    double com_z = 0.0;
    for (Atom atom : atoms)
    {
        com_x += atom.xx * atom.mass;
        com_y += atom.yy * atom.mass;
        com_z += atom.zz * atom.mass;
    }
    com_x /= n_atoms;
    com_y /= n_atoms;
    com_z /= n_atoms;
    for (Atom atom : atoms)
    {
        atom.xx -= com_x;
        atom.yy -= com_y;
        atom.zz -= com_z;
    }
}

bool Molecule::SpinChargeValidate()
{
    total_electrons = 0;
    for (Atom atom : atoms)
    {
        total_electrons += atom.electron_counter;
    }
    total_electrons -= charge;
    if (total_electrons % 2)
    {
        if (spin_mult % 2)
        {
            return false;
        }
        return true;
    }
    if (spin_mult % 2)
    {
        return true;
    }
    return false;
}

void Molecule::FindBonds()
{
    //bonds (vector of vector of ints)
    for (unsigned int i = 0; i < atoms.size()-1; i++)
    {
        double x1 = atoms[i].xx;
        double y1 = atoms[i].yy;
        double z1 = atoms[i].zz;
        double vdw_1 = atoms[i].vdw_radius;
        for (unsigned int j = i+1; j < atoms.size(); j++)
        {
            double x2 = atoms[j].xx;
            double y2 = atoms[j].yy;
            double z2 = atoms[j].zz;
            double vdw_2 = atoms[i].vdw_radius;
            double dist = sqrt( pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2) );
            double vdw_dist = (vdw_2 + vdw_1)/2; //average of vdw radii?
            if (dist < vdw_dist)
            {
                // based on VDW radii distances, 
                bonds.push_back({(int)i,(int)j});
            }
        }
    }
    for (std::vector<int> bond : bonds)
    {
        int a = bond[0];
        int b = bond[1];
        atoms[a].bonded_to_indexes.push_back(b);
        atoms[a].bonded_to_elements.push_back(atoms[b].element);
        atoms[b].bonded_to_indexes.push_back(a);
        atoms[b].bonded_to_elements.push_back(atoms[a].element);
    }
}

void Molecule::FindAngles()
{
    std::vector<int> angle_tmp;
    for (unsigned int i = 0; i < bonds.size()-1; i++)
    {
        int a1 = bonds[i][0];
        int a2 = bonds[i][1];
        
        for (unsigned int j = i+1; j < bonds.size(); j++)
        {
            int a3 = bonds[j][0];
            int a4 = bonds[j][1];
            if ((a1 == a3)||(a2 == a3)||(a1==a4)||(a2==a4))
            {
                if (a1 == a3)
                {
                    angle_tmp = {a2,a1,a4};
                }
                else if (a2 == a3)
                {
                    angle_tmp = {a1,a2,a4};
                }
                else if (a1 == a4)
                {
                    angle_tmp = {a2,a1,a3};
                }
                else if (a2 == a4)
                {
                    angle_tmp = {a1,a2,a3};
                }
                angles.push_back(angle_tmp);
            }
        }
    }
}

void Molecule::FindTorsions()
{
    std::vector <int> torsion_tmp;
    //torsions (vector of vector of ints)

    //iterate through angles, match ends to bonds with non-middle atom names.
    for (unsigned int a=0; a < angles.size(); a++)
    {   
        int a1 = angles[a][0];
        int a2 = angles[a][1];
        int a3 = angles[a][2];
        for (unsigned int b=0;b<bonds.size(); b++)
        {
            int b1 = bonds[b][0];
            int b2 = bonds[b][1];
            if ( (a1 == b1) && (a2 != b2) ) //check the middle atom in the angle is NOT the same as the second atom in the bond.
            {
                torsion_tmp = {b2,a1,a2,a3};
            }
            else if ( (a1 == b2) && (a2 != b1) ) //check the middle atom in the angle is NOT the same as the second atom in the bond.
            {
                torsion_tmp = {b1,a1,a2,a3};
            }
            else if ( (a3 == b1) && (a2 != b2) )
            {
                torsion_tmp = {a1,a2,a3,b2};
            }
            else if ( (a3 == b2) && (a2 != b1) )
            {
                torsion_tmp = {a1,a2,a3,b1};
            }
            else
            {
                continue;
            }
            torsions.push_back(torsion_tmp);
        }
    }
}

void Molecule::FindDihedrals()
{
    // angles with matching middle to other bonds.
    std::vector <int> dihedral_tmp;
    //dihedrals (vector of vector of ints)
       //iterate through angles, match ends to bonds with non-middle atom names.
    for (unsigned int a=0; a < angles.size(); a++)
    {   
        int a1 = angles[a][0];
        int a2 = angles[a][1];
        int a3 = angles[a][2];
        for (unsigned int b=0;b<bonds.size(); b++)
        {
            int b1 = bonds[b][0];
            int b2 = bonds[b][1];
            if ( (a2 != b1) && (a2 != b2) ) // If neither atom in the bond is the middle atom of the angle, continue to next bond.
            {
                continue;
            }
            else if ( (a1==b2) || (a3==b2) || (a1==b1) || (a3==b1) ) // if any bond atoms are in any of the ends, continue to next bond.
            {
                continue;
            }
            else if (a2 == b1) // if b1 is the middle atom, add b2 to the list for dihedrals.
            {
                dihedral_tmp = {a1,a2,a3,b2};
            }
            else if (a2 == b2) // if b2 is the middle atom, add b1 to the list for dihedrals.
            {
                dihedral_tmp = {a1,a2,a3,b1};
            }
            else
            {
                continue;
            }
            dihedrals.push_back(dihedral_tmp);
        }
    }
}

void Molecule::FindRings()
{
    //iterate through all torsions (A-B-C-D) in the molecule
    for (unsigned int i=0; i < torsions.size(); i++)
    {   
        //assign atom numbers in torsion
        int at_A = torsions[i][0];
        int at_B = torsions[i][1];
        int at_C = torsions[i][2];
        int at_D = torsions[i][3];
        
        for (std::vector <int> bond : bonds) // iterate through bonds to identify all 4-member rings.
        {
            int at_E = bond[0];
            int at_F = bond[1];
            if (((at_A == at_E) && (at_D == at_F)) || ((at_A == at_F) && (at_D == at_E)) )
            {
                // If torsion is 4-membered ring, skip subsequent comparisons for this bond.
                four_rings.push_back({at_A,at_B,at_C,at_D});
                break;
            }
        }

        for (std::vector <int> angle : angles) // iterate through angles against torsions to identify all 5 member rings.
        {
            int at_E = angle[0];
            int at_F = angle[1];
            int at_G = angle[2];
            if (((at_A == at_E) && (at_D == at_G)) || ((at_A == at_G) && (at_D == at_E)))
            {
                if ((at_F != at_B) && (at_F != at_C))
                {
                        five_rings.push_back({at_A,at_B,at_C,at_D,at_F});
                }
            }
        }

        for (unsigned int j=i+1; j < torsions.size(); j++) // iterate through torsions to identify all 6-member rings.
        {
            /*
               B--C
              /    \           B--C
             A      D         /    \
                vs.    -->   A      D --> A,B,C,D,G,F
             E      H         \    /
              \    /           F--G
               F--G
            */
            int at_E = torsions[j][0];
            int at_F = torsions[j][1];
            int at_G = torsions[j][2];
            int at_H = torsions[j][3];
            if ((at_A == at_E) && (at_D == at_H))
            {
                if ((at_F != at_B) && (at_F != at_C) && (at_G != at_B) && (at_G != at_C))
                {
                    six_rings.push_back({at_A,at_B,at_C,at_D,at_G,at_F});

                }
            }
            else if ((at_A == at_H) && (at_D == at_E))
            {
                if ((at_F != at_B) && (at_F != at_C) && (at_G != at_B) && (at_G != at_C))
                {
                    six_rings.push_back({at_A,at_B,at_C,at_D,at_F,at_G});
                }
            }
        }
    }

    // gotta check if any rings are duplicates...
}