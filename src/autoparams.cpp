#include "autoparams.h"
#include "parameters.h"

void write_TC_resp_input(Settings settings, Molecule mol)
{
    std::stringstream buffer;
    if (settings.dummy_atom_names.size() > 0)
    {
        std::string line;
        buffer.str("");
        buffer << "cp " << settings.inputfile << " capped.pdb";
        silent_shell(buffer.str().c_str());
        
        std::ofstream outfile(settings.job_dir + settings.inputfile);
        for (Atom atom : mol.atoms)
        {   if (find(settings.dummy_atom_names.begin(),settings.dummy_atom_names.end(),atom.atom_name)  != settings.dummy_atom_names.end())
            {
                continue;
            }
            outfile << atom.print_PDB_line();
        }
        outfile << "END\n";
        outfile.close();
    }
    std::vector<std::string> tmp ={"coordinates", "charge", "spinmult", "basis", "method", "convthre", "threall", "precision", "maxit", "scf", "gpus", "gpumem", "scrdir","run", "resp"};
    std::map<std::string,std::string> tc_keywords={  
        {"coordinates",settings.inputfile},
        {"charge",std::to_string(settings.mol_charge)},
        {"spinmult",std::to_string(settings.mol_spin)},
        {"basis","6-31gss"},
        {"method","b3lyp"},
        {"convthre","1e-7"},
        {"threall","1e-14"},
        {"precision","mixed"},
        {"maxit","200"},
        {"scf","diis+a"},
        {"gpus","1"},
        {"gpumem","256"},
        {"scrdir","scr/"},
        {"run","energy"},
        {"resp","yes"}};
    //Write the TeraChem input to calculate resp charges.
    
    
    std::ofstream outfile;
    outfile.open(settings.job_dir + "/resp.in",std::ios::out);
    for (std::string key : tmp)
    {
        buffer.str("");
        buffer << std::setw(20) << std::left << key;
        buffer << std::setw(20) << std::right << tc_keywords[key] << std::endl;
        outfile << buffer.str();
    }
    // for ( int i =0; i < tmp.size(); i++ )
    // {
    //     buffer.str("");
    //     buffer << std::setw(20) << std::left << tmp[i];
    //     buffer << std::setw(20) << std::right << tc_keywords[tmp[i]] << std::endl;
    //     outfile << buffer.str();
    // } 
    outfile.close();
}

bool validate_TC_resp_output(Settings settings)
{
    std::ifstream respfile(settings.job_dir + "/resp.out",std::ios::in);
    std::string line;
    std::stringstream buffer;
    if ( !respfile.is_open())
    {
        settings.Error("Unable to open output file from RESP calculation.");
        return false;
    }
    while (getline(respfile,line))
    {
        if (line.find("Job terminated") != std::string::npos)
        {
            settings.Error("RESP Calculation did not complete properly.  Be sure to check charges in the mol2 file.");
            return false;
        }
        else if (line.find("Job finished") != std::string::npos)
        {
            settings.Output("RESP Calculation completed. Parsing RESP charges now.");
            DeleteIfEmpty(settings.job_dir + "/resp.err");
            return true;
        }
    }
    return false;
}

void parse_TC_resp_output(Molecule &mol,Settings settings)
{
    std::vector<double> resp_charges;
    std::ifstream respfile(settings.job_dir + "/resp.out",std::ios::in);
    std::string line;
    std::stringstream buffer;
    std::string atom,x,y,z,charge;

    if ( !respfile.is_open())
    {
        settings.Error("Unable to open output file from RESP calculation.");
        return;
    }
    // Parse the output file for RESP charges.
    while (getline(respfile,line))
    {
        if (line.find("ESP restraint charges") != std::string::npos)
        {
            getline(respfile,line);
            getline(respfile,line);
            break;
        }
    }
    int atom_count=0;
    while (getline(respfile,line))
    {
        if (line.find("------") != std::string::npos)
        {
            break;
        }
        buffer.str(line);
        buffer >> atom >> x >> y >> z >> charge;
        buffer.flush();
        buffer.str("");
        mol.atoms[atom_count].set_resp_charge(stof(charge));
        atom_count++;
    }
    return;
}

void Generate_Mol2_File(Settings settings)
{
    std::stringstream buffer;
    std::string curr_path = fs::current_path();
    fs::current_path(settings.job_dir);
    
    // run antechamber
    buffer.str("");
    buffer << "antechamber -i " << settings.inputfile << " -fi pdb -at amber -o tmp.mol2 -fo mol2 -pf y 1> antechamber.log 2> antechamber.err";
    DeleteIfEmpty(settings.job_dir + "/antechamber.err");
    silent_shell(buffer.str().c_str());
    Validate_Mol2_File(settings);
    fs::current_path(curr_path);
    return;
}

void Validate_Mol2_File(Settings settings)
{
    settings.Output("Validating mol2 file...");
    std::ifstream tmp_mol2;
    std::ofstream new_mol2;
    std::string line;
    std::stringstream buffer;
    int atom_count=0;
    std::string curr_path = fs::current_path();
    fs::current_path(settings.job_dir);
    tmp_mol2.open("tmp.mol2",std::ios::in);
    new_mol2.open("tmp2.mol2",std::ios::out);
    if (!tmp_mol2.is_open())
    {
        settings.Error("Unable to open original mol2.");
        return;
    }
    if (!new_mol2.is_open())
    {
        settings.Error("Unable to write new mol2.");
        return;
    }
    // Starting line:  @<TRIPOS>ATOM
    while (getline(tmp_mol2,line))
    {
        new_mol2 << line << std::endl;
        if (line.find("@<TRIPOS>ATOM") != std::string::npos)
        {
            break;
        }
    }
    while (getline(tmp_mol2,line))
    {
        if (line.find("@<TRIPOS>BOND") != std::string::npos)
        {
            // Ending line:  @<TRIPOS>BOND
            new_mol2 << line << std::endl;
            break;
        }
        if (line.find(" DU ") != std::string::npos)
        {
            settings.Output("Found dummy atom (DU).  Attempting to replace with valid atom type.");
            settings.Output(line);
            // Fix the DU atom type
            std::string dummy;
            std::stringstream new_atom_type;
            new_atom_type.str("");
            buffer.str(line);
            buffer >> dummy >> dummy; // first item in the line is the atom index, second is the atom name.
            if (isdigit(dummy[1]))
            {
                new_atom_type << dummy.substr(0,1) << " ";
            }
            else
            {
                new_atom_type << dummy.substr(0,2);
            }
            new_mol2 << substring_replace(line, "DU", new_atom_type.str()) << std::endl;
            settings.Output("Replaced DU with "+new_atom_type.str());
        }
        else
        {
            new_mol2 << line << std::endl;
        }
        
    }
  
    while (getline(tmp_mol2,line))
    {
        new_mol2 << line << std::endl;
    }
    tmp_mol2.close();
    new_mol2.close();
    silent_shell("mv tmp2.mol2 tmp.mol2");
    fs::current_path(curr_path);
    return;
}

void Add_Charges_To_Mol2(Settings settings, Molecule &mol)
{   
    std::ifstream tmp_mol2;
    std::ofstream new_mol2;
    std::string line;
    std::stringstream buffer;
    int atom_count=0;
    std::string curr_path = fs::current_path();
    fs::current_path(settings.job_dir);
    tmp_mol2.open("tmp.mol2",std::ios::in);
    new_mol2.open(settings.mol2file,std::ios::out);
    if (!tmp_mol2.is_open())
    {
        settings.Error("Unable to open original mol2.");
        return;
    }
    if (!new_mol2.is_open())
    {
        settings.Error("Unable to write new mol2.");
        return;
    }

    while (getline(tmp_mol2,line))
    {   
        int location = line.find("0.000000");
        if (location != std::string::npos)
        {
            buffer.str("");
            buffer << line.substr(0,location) << std::fixed << std::right << mol.atoms[atom_count].resp_charge;
            mol.atoms[atom_count].atom_type = line.substr(50,2);
            new_mol2 << buffer.str() << std::endl;
            atom_count++;
        }
        else
        {
            new_mol2 << line << std::endl;
        }
        
    }
    if ((mol.head_atom_name != "0") || (mol.tail_atom_name != "0"))
    {
        new_mol2 << "@<TRIPOS>HEADTAIL" << std::endl;
        if (mol.head_atom_name == "0")
        {
            new_mol2 << "0 0" << std::endl;
        }
        else
        {
            new_mol2 << mol.head_atom_name << " 1" << std::endl;
        }
        if (mol.tail_atom_name == "0")
        {
            new_mol2 << "0 0" << std::endl;
        }
        else
        {
            new_mol2 << mol.tail_atom_name << " 1"  << std::endl;
        }
        new_mol2 << "@<TRIPOS>RESIDUECONNECT" << std::endl;
        new_mol2 << "1 " << mol.head_atom_name << " " << mol.tail_atom_name << " 0 0 0 0" << std::endl << std::endl;
    }
    tmp_mol2.close();
    new_mol2.close();
    silent_shell("rm tmp.mol2");

    // Strip dummy atoms from mol2 using CPPTRAJ
    //     parm o37-with-charges-reordered.mol2 name sugar
    // trajin o37-with-charges-reordered.mol2 parm sugar
    // strip @O43,C20,1H20,2H20,3H20
    // strip @C17,1H17,2H17,3H17 charge -1.0
    // trajout o37-fragment.mol2 mol2
    // run
    fs::current_path(curr_path);
    return;
}

void Check_For_Missing_Parameters(Settings settings, Molecule &mol)
{
    std::string curr_path = fs::current_path();
    fs::current_path(settings.job_dir);

    std::map<std::string,std::string> leaprc_dict={{"DNA","source leaprc.DNA.OL15"},{"RNA","source leaprc.RNA.OL3"},{"PROTEIN","source leaprc.protein.ff14SB"},{"CARBOHYDRATE","source leaprc.GLYCAM_06j-1"}};
    std::stringstream buffer;
    std::ofstream tleap_in;
    tleap_in.open("tleap.in",std::ios::out);
    if (!tleap_in.is_open())
    {
        settings.Error("Unable to write tleap.in file.");
        return;
    }

    // Write tleap.in file.
    tleap_in << leaprc_dict[settings.forcefield] << std::endl;
    tleap_in << "source leaprc.water.tip3p" << std::endl;
    tleap_in << mol.res_name << " = loadmol2 " << settings.mol2file << std::endl;
    if (CheckFileExists(settings.frcmodfile))
    {
        tleap_in << "loadamberparams " << settings.frcmodfile << std::endl;
    }    
    
    // Load system depends on connecting atoms
    if ((settings.head_atom_name == "0") && (settings.tail_atom_name == "0"))// if no connecting atoms
    {
        settings.Output("Validating non-connecting molecule.");
        tleap_in << "test = loadpdb " << settings.inputfile << std::endl;
    }
    else if (settings.tail_atom_name == "0")// if only head connecting atom
    {
        settings.Output("Validating head-connecting molecule.");
        if (settings.forcefield == "PROTEIN")
        {
            tleap_in << "test = sequence { ALA " << mol.res_name << " } " << std::endl;
        }
        else if (settings.forcefield == "RNA")
        {
            tleap_in << "test = sequence { A5 " << mol.res_name << " } " << std::endl;
        }
        else if (settings.forcefield == "DNA")
        {
            tleap_in << "test = sequence { DA5 " << mol.res_name << " } " << std::endl;
        }
        else
        {
            tleap_in << "test = sequence { " << mol.res_name << " } " << std::endl;
        }
    }
    else if (settings.head_atom_name == "0")// if only tail connecting atom
    {
        settings.Output("Validating tail-connecting molecule.");
        if (settings.forcefield == "PROTEIN")
        {
            tleap_in << "test = sequence { " << mol.res_name << " ALA } " << std::endl;
        }
        else if (settings.forcefield == "RNA")
        {
            tleap_in << "test = sequence { " << mol.res_name << " A3 } " << std::endl;
        }
        else if (settings.forcefield == "DNA")
        {
            tleap_in << "test = sequence { " << mol.res_name << " DA3 } " << std::endl;
        }
        else
        {
            tleap_in << "test = sequence { " << mol.res_name << " } " << std::endl;
        }
    }
    else  // if both head and tail connecting atoms
    {
        settings.Output("Validating multi-connecting molecule.");
        if (settings.forcefield == "PROTEIN")
        {
            tleap_in << "test = sequence { ALA " << mol.res_name << " ALA } " << std::endl;
        }
        else if (settings.forcefield == "RNA")
        {
            tleap_in << "test = sequence { A5 " << mol.res_name << " A3 } " << std::endl;
        }
        else if (settings.forcefield == "DNA")
        {
            tleap_in << "test = sequence { DA5 " << mol.res_name << " DA3 } " << std::endl;
        }
        else
        {
            tleap_in << "test = sequence { " << mol.res_name << " } " << std::endl;
        }
    }
   
    tleap_in << "addions test Na+ 0" << std::endl;
    tleap_in << "addions test Cl- 0" << std::endl;
    tleap_in << "solvatebox test TIP3PBOX 12.0" << std::endl;
    tleap_in << "saveamberparm test tmp.prmtop tmp.rst7" << std::endl;
    tleap_in << "quit" << std::endl;
    tleap_in.close();

    // Run tleap with tleap.in
    buffer.str("");
    buffer << "tleap -f tleap.in 1> tleap.out 2> tleap.err";
    silent_shell(buffer.str().c_str());
    DeleteIfEmpty(settings.job_dir + "/tleap.err");
    buffer.str("");

    std::ifstream leaplog;
    leaplog.open("tleap.out",std::ios::in);
    if (!leaplog.is_open())
    {
        settings.Error("Unable to open tleap output file.");
        return;
    }
    std::string line;
    while(getline(leaplog,line))
    {
        if (line.find("Could not find bond parameter for atom types:") != std::string::npos)
        {
            std::stringstream bondline;
            bondline.str(line.substr(46,line.size()-46));
            // std::cout << "found bond string: " << bondline.str() << std::endl;  //DEBUGGING OUTPUT
            std::string tmp = trim_whitespace(bondline.str());
            std::vector<std::string> atom_segments = string_split(tmp, '-');
            atom_segments[0] = trim_whitespace(atom_segments[0]);
            atom_segments[1] = trim_whitespace(atom_segments[1]);
            // std::cout << atom_segments[0] << ", " << atom_segments[1] << std::endl;  //DEBUGGING OUTPUT

            std::stringstream bond_pair;
            bond_pair.str("");

            if (atom_segments[0] < atom_segments[1])
            {
                bond_pair << std::setw(2) << std::left  << atom_segments[0] << "-" << std::setw(2) << std::left  << atom_segments[1];
            }
            else
            {
                bond_pair << std::setw(2) << std::left  << atom_segments[1] << "-" << std::setw(2) << std::left  << atom_segments[0];
            }
            // std::cout << bond_pair.str() << std::endl;  //DEBUGGING OUTPUT
            if (count(mol.missing_bonds.begin(), mol.missing_bonds.end(), bond_pair.str()) == 0)
            {
                mol.missing_bonds.push_back(bond_pair.str());
                buffer.str("");
                buffer << "missing bond between " << atom_segments[0] << " and " << atom_segments[1] << std::endl;
                settings.Output(buffer.str());
                buffer.str("");
            }
        }
        else if (line.find("Could not find bond parameter for:") != std::string::npos)
        { // bonds
            std::stringstream bondline;
            bondline.str(line.substr(35,line.size()-35));
            // std::cout << "found bond string: " << bondline.str() << std::endl;  //DEBUGGING OUTPUT
            std::string tmp = trim_whitespace(bondline.str());
            std::vector<std::string> atom_segments = string_split(tmp, '-');
            atom_segments[0] = trim_whitespace(atom_segments[0]);
            atom_segments[1] = trim_whitespace(atom_segments[1]);
            // std::cout << atom_segments[0] << ", " << atom_segments[1] << std::endl;  //DEBUGGING OUTPUT

            std::stringstream bond_pair;
            bond_pair.str("");

            if (atom_segments[0] < atom_segments[1])
            {
                bond_pair << std::setw(2) << std::left  << atom_segments[0] << "-" << std::setw(2) << std::left  << atom_segments[1];
            }
            else
            {
                bond_pair << std::setw(2) << std::left  << atom_segments[1] << "-" << std::setw(2) << std::left  << atom_segments[0];
            }
            // std::cout << bond_pair.str() << std::endl;  //DEBUGGING OUTPUT
            if (count(mol.missing_bonds.begin(), mol.missing_bonds.end(), bond_pair.str()) == 0)
            {
                mol.missing_bonds.push_back(bond_pair.str());
                buffer.str("");
                buffer << "missing bond between " << atom_segments[0] << " and " << atom_segments[1] << std::endl;
                settings.Output(buffer.str());
                buffer.str("");
            }
        }
        else if (line.find("Could not find angle parameter for atom types:") != std::string::npos)
        {
            std::stringstream bondline;
            bondline.str(line.substr(47,line.size()-47));
            // std::cout << "found angle string: " << bondline.str() << std::endl;
            std::vector<std::string> atom_segments = string_split(bondline.str(), '-');
            std::stringstream bond_pair;
            bond_pair.str("");

            
            if (atom_segments[0] < atom_segments[2])
            {
                bond_pair << std::setw(2) << std::left  << atom_segments[0] << "-" << std::setw(2) << std::left  << atom_segments[1] << "-" << std::setw(2) << std::left  << atom_segments[2];
            }
            else
            {
                bond_pair << std::setw(2) << std::left  << atom_segments[2] << "-" << std::setw(2) << std::left  << atom_segments[1] << "-" << std::setw(2) << std::left  << atom_segments[0];
            }
            if (count(mol.missing_angles.begin(), mol.missing_angles.end(), bond_pair.str()) == 0)
            {
                mol.missing_angles.push_back(bond_pair.str());
                buffer.str("");
                buffer << "missing angle between " << atom_segments[0] << ", " << atom_segments[1] << ", and " << atom_segments[2] <<  std::endl;
                settings.Output(buffer.str());
                buffer.str("");
            }
        }
        else if (line.find("Could not find angle parameter:") != std::string::npos)
        { // angles
            std::stringstream bondline;
            bondline.str(line.substr(31,line.size()-31));
            // std::cout << "found angle string: " << bondline.str() << std::endl;
            std::vector<std::string> atom_segments = string_split(bondline.str(), '-');
            std::stringstream bond_pair;
            bond_pair.str("");

            
            if (atom_segments[0] < atom_segments[2])
            {
                bond_pair << std::setw(2) << std::left  << atom_segments[0] << "-" << std::setw(2) << std::left  << atom_segments[1] << "-" << std::setw(2) << std::left  << atom_segments[2];
            }
            else
            {
                bond_pair << std::setw(2) << std::left  << atom_segments[2] << "-" << std::setw(2) << std::left  << atom_segments[1] << "-" << std::setw(2) << std::left  << atom_segments[0];
            }
            if (count(mol.missing_angles.begin(), mol.missing_angles.end(), bond_pair.str()) == 0)
            {
                mol.missing_angles.push_back(bond_pair.str());
                buffer.str("");
                buffer << "missing angle between " << atom_segments[0] << ", " << atom_segments[1] << ", and " << atom_segments[2] <<  std::endl;
                settings.Output(buffer.str());
                buffer.str("");
            }
        }
        
        else if (line.find(" ** No torsion terms for atom types") != std::string::npos)
        {  // dihedrals
            std::stringstream bondline;
            bondline.str(line.substr(36,line.size()-36));
            // std::cout << "found dihedral string: " << bondline.str() << std::endl;
            std::vector<std::string> atom_segments = string_split(bondline.str(), '-');

            std::stringstream bond_pair;
            bond_pair.str("");
            if (atom_segments[1] < atom_segments[2])
            {
                bond_pair << std::setw(2) << std::left << atom_segments[0] << "-" << std::setw(2) << std::left  << atom_segments[1] << "-" << std::setw(2) << std::left  << atom_segments[2] << "-" << std::setw(2) << std::left  << atom_segments[3];
            }
            else
            {
                bond_pair << std::setw(2) << std::left << atom_segments[3] << "-" << std::setw(2) << std::left  << atom_segments[2] << "-" << std::setw(2) << std::left  << atom_segments[1] << "-" << std::setw(2) << std::left  << atom_segments[0];
            }
            // std::cout << "Processed dihedral stringstream is: " << bond_pair.str() << std::endl;
            if (count(mol.missing_dihedrals.begin(), mol.missing_dihedrals.end(), bond_pair.str()) == 0)
            {
                mol.missing_dihedrals.push_back(bond_pair.str());
                buffer.str("");
                buffer << "missing dihedral between " << atom_segments[0] << ", " << atom_segments[1] <<  ", " << atom_segments[2] << ", and " << atom_segments[3] <<  std::endl;
                settings.Output(buffer.str());
                buffer.str("");
            }
        }
        else if (line.find(" ** No torsion terms for  ") != std::string::npos)
        {  // dihedrals
            std::stringstream bondline;
            bondline.str(line.substr(26,line.size()-26));
            // std::cout << "found dihedral string: " << bondline.str() << std::endl;
            std::vector<std::string> atom_segments = string_split(bondline.str(), '-');
            std::stringstream bond_pair;
            bond_pair.str("");
            if (atom_segments[1] < atom_segments[2])
            {
                bond_pair << std::setw(2) << std::left << atom_segments[0] << "-" << std::setw(2) << std::left  << atom_segments[1] << "-" << std::setw(2) << std::left  << atom_segments[2] << "-" << std::setw(2) << std::left  << atom_segments[3];
            }
            else
            {
                bond_pair << std::setw(2) << std::left << atom_segments[3] << "-" << std::setw(2) << std::left  << atom_segments[2] << "-" << std::setw(2) << std::left  << atom_segments[1] << "-" << std::setw(2) << std::left  << atom_segments[0];
            }
            // std::cout << "Processed dihedral stringstream is: " << bond_pair.str() << std::endl;
            if (count(mol.missing_dihedrals.begin(), mol.missing_dihedrals.end(), bond_pair.str()) == 0)
            {
                mol.missing_dihedrals.push_back(bond_pair.str());
                buffer.str("");
                buffer << "missing dihedral between " << atom_segments[0] << ", " << atom_segments[1] <<  ", " << atom_segments[2] << ", and " << atom_segments[3] <<  std::endl;
                settings.Output(buffer.str());
                buffer.str("");
            }
        }

        else if (line.find("No sp2 improper torsion term for") != std::string::npos)
        { // improper torsions

        }

        else if (line.find("could not find vdW (or other) parameters for type (") != std::string::npos)
        { // vdw parameters
            std::stringstream vdwline;
            vdwline.str(line.substr(line.find_last_of('(')+1,2));
            std::cout << "VDW ATOM FOUND: " << vdwline.str() << std::endl;
            mol.missing_vdw_params.push_back(vdwline.str());
        }
        
    }
    fs::current_path(curr_path);
    return;
}

void Build_Frcmod(Frcmod_File &frcmod, Molecule mol, Settings settings)
{
    // process molecule's missing masses
    for (int i=0; i < mol.missing_masses.size(); i++)
    {
        
        // use mass library
        // identify most likely mass line to use OR add a new mass and alert user.
    }

    // process molecule's missing bonds
    for (int i=0; i < mol.missing_bonds.size(); i++)
    {
        // std::cout << mol.missing_bonds[i] << std::endl;
        std::string a1 = mol.missing_bonds[i].substr(0,2);
        std::string a2 = mol.missing_bonds[i].substr(3,2);
        // std::cout << "Atoms: " << a1 << " and " << a2 << "... searching" << std::endl;
        std::string newbond = FindMissingBond(a1,a2);
        // std::cout << "got newbond: " << newbond << std::endl;
        frcmod.Add_Bond(newbond);
        
        //use bond library first
        //if not able to find exact match, use BOND_LENGTH_DICT and BOND_FORCE_CONSTANT_DICT against BASIC_ATOM_TYPES string pair.
    }

    // process molecule's missing angles
    for (int i=0; i < mol.missing_angles.size(); i++)
    {
        // std::cout << mol.missing_angles[i] << std::endl;
        std::vector<std::string> atom_set = string_split(mol.missing_angles[i],'-');

        std::string a1 = trim_whitespace(atom_set[0]);
        std::string a2 = trim_whitespace(atom_set[1]);
        std::string a3 = trim_whitespace(atom_set[2]);
        if (a1.size() == 1)
        {
            a1 = a1+" ";
        }
        if (a2.size() == 1)
        {
            a2 = a2+" ";
        }
        if (a3.size() == 1)
        {
            a3 = a3+" ";
        }
        // std::cout << "Atoms: " << a1 << ", " << a2 << ", and " << a3 << "... searching" << std::endl;
        std::string newangle = FindMissingAngle(a1,a2,a3);
        // std::cout << "got newangle: " << newangle << std::endl;
        frcmod.Add_Angle(newangle);
        //use angle library first
        //if not able to find exact match, use ANGLE_MEASURE_DICT and ANGLE_FORCE_CONSTANT_DICT against BASIC_ATOM_TYPES string pair.
    }

    // process molecule's missing dihedrals
    for (int i=0; i < mol.missing_dihedrals.size(); i++)
    {
        //use dihedrals library first
        // std::cout << mol.missing_dihedrals[i] << std::endl;
        std::vector<std::string> atom_set = string_split(mol.missing_dihedrals[i],'-');

        std::string a1 = trim_whitespace(atom_set[0]);
        std::string a2 = trim_whitespace(atom_set[1]);
        std::string a3 = trim_whitespace(atom_set[2]);
        std::string a4 = trim_whitespace(atom_set[3]);
        if (a1.size() == 1)
        {
            a1 = a1+" ";
        }
        if (a2.size() == 1)
        {
            a2 = a2+" ";
        }
        if (a3.size() == 1)
        {
            a3 = a3+" ";
        }
        if (a4.size() == 1)
        {
            a4 = a4+" ";
        }
        // std::cout << "Atoms: " << a1 << ", " << a2 << ", " << a3 << ", and " << a4 << "... searching" << std::endl;
        std::string newdihedral = FindMissingDihedral(a1,a2,a3,a4);
        // std::cout << "got newdihedral: " << newdihedral << std::endl;
        frcmod.Add_Dihedral(newdihedral);
        //if not able to find exact match, use DIHEDRAL_SET_DICT against BASIC_ATOM_TYPES string pair.
    }

    // process molecule's missing torsions
    for (int i=0; i < mol.missing_torsions.size(); i++)
    {
        //use dihedrals library first
        //if not able to find exact match, use DIHEDRAL_SET_DICT against BASIC_ATOM_TYPES string pair.
    }

    // process molecule's missing vdW params
    for (int i=0; i < mol.missing_vdw_params.size(); i++)
    {
        //use nonbonded library first
        frcmod.Add_vdW(FindMissingNonbonded(mol.missing_vdw_params[i]));
        frcmod.Add_Mass(FindMissingMass(mol.missing_vdw_params[i]));
        //if not able to find exact match, create new nonbonded term.
    }


    // write frcmod file
    frcmod.Write_Frcmod_File(settings);
}