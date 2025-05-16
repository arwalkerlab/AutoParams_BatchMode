#include "autoparams.h"
#include "utilities.h"
#include "classes.h"
#include "settings.h"
#include "parameters.h"
#include "optimize.h"
#include "mol2.h"

int main(int argc, char **argv)
{
    // Variable Declarations
    std::stringstream buffer;

    // Initialize program, figure out what is available for use.
    Settings settings;
    settings.parse_command_line(argc, argv);

    // Initialize secondary variables.
    Frcmod_File frcmod(settings.frcmodfile);
    
    // read input file (PDB)
    settings.Output("Reading PDB...");
    Molecule mol(settings);
    settings.Output("Checked PDB for unique atom names.");
    settings.Output("Moving molecule to center of mass...");

    // move molecule to center of mass.
    mol.move_to_COM();
    mol.Write_PDB(settings.job_dir, settings.inputfile);

    // check total charge/spin combination
    if (mol.SpinChargeValidate())
    {
        settings.Output("Spin/charge combination validated.\n\n");
        buffer.str("");
        buffer << "Spin:   " << mol.spin_mult << "\n";
        buffer << "Charge: " << mol.charge << "\n";
        buffer << "Total Electron Count:" << mol.total_electrons << " \n";
        settings.Output(buffer.str());
    }
    else
    {
        buffer.str("");
        buffer << "Given spin/charge combination is not valid. \n";
        buffer << "Spin:   " << mol.spin_mult << "\n";
        buffer << "Charge: " << mol.charge << "\n";
        buffer << "Total Electron Count:" << mol.total_electrons << " \n";
        settings.Error(buffer.str());
    }

    settings.Output("Analyzing molecular structure...\n\n");
    // identify all bonds
    // by distance?
    mol.FindBonds();
    buffer.str("");
    buffer << "Found " << mol.bonds.size() << " bonds in molecule." << std::endl;
    settings.Output(buffer.str());
    
    // identify all angles
    // by bondpairs with one matching atom
    mol.FindAngles();
    buffer.str("");
    buffer << "Found " << mol.angles.size() << " angles in molecule." << std::endl;
    settings.Output(buffer.str());

    // identify all dihedrals
    // angles with matching middle to other bonds.
    mol.FindTorsions();
    buffer.str("");
    buffer << "Found " << mol.torsions.size() << " torsions in molecule." << std::endl;
    settings.Output(buffer.str());

    // identify all torsions
    // angles with matching ends
    mol.FindDihedrals();
    buffer.str("");
    buffer << "Found " << mol.dihedrals.size() << " dihedrals in molecule." << std::endl;
    settings.Output(buffer.str());
   
    // identify all rings
    mol.FindRings();
    buffer.str("");
    int total_rings = mol.four_rings.size() + mol.five_rings.size() + mol.six_rings.size();
    buffer << "Found " << total_rings << " rings in molecule.  Ring sizes are: " << std::endl;
    buffer << "Four-membered rings: " << mol.four_rings.size() << std::endl;
    buffer << "Five-membered rings: " << mol.five_rings.size() << std::endl;
    buffer << "Six-membered rings:  " << mol.six_rings.size() << std::endl;
    settings.Output(buffer.str());

    // OPTIONAL optimize structure
    if (settings.OPTIMIZE_FIRST) //TODO:  Add TeraChem optimization function.  Quick and dirty, nothing fancy, just get it done.
    {
        // Run TeraChem optimization (flag for GS/ES opt?)
        // confirm successful job completion, if failed, bail out?
        // copy optimized structure to original filename in job_dir.
        settings.Output("Running QM optimization before parametrizing.  This may take some time.");
        TeraChemOpt(settings);
    }
        
    // calculate RESP charges
    // calculate RESP charges, apply results to atoms in molecule.
    write_TC_resp_input(settings,mol);
    buffer.str("");
    buffer << "Running TeraChem for RESP charges.  This may take a while..." << std::endl;
    settings.Output(buffer.str());
    std::string curr_path = fs::current_path();
    fs::current_path(settings.job_dir);
    silent_shell("terachem -i resp.in 1> resp.out 2> resp.err");
    fs::current_path(curr_path);

    // validate successful resp job before trying to parse output.
    // If validation fails, there could be an issue with convergence, but the charges may still be useful.  
    // We will alert the end user and continue on.  
    if (!validate_TC_resp_output(settings))
    {   
        settings.Output("RESP calculation did not complete successfully.  Attempting to continue with current charge values.");
    }
    
    // parse RESP charges
    parse_TC_resp_output(mol, settings);

    // ************************* Here's where I can put in my molecule atom-typer code...
    AtomTyping(mol);
    
    // Generate mol2 file
    Generate_Mol2_File(settings);
    Add_Charges_To_Mol2(settings, mol);

    // check for missing parameters with tleap, parse tleap error file.
    settings.Output("\nIdentifying missing parameters...\n#########################################\n");
    Check_For_Missing_Parameters(settings, mol);

    // OPTIONAL calculate parameters with TeraChem optimizations
        // Each bond should be tested with 0.1A increments to get bond strength
        // Each angle should be tested with 10 degree increments (0-180)
        // Each dihedral should be tested with 10 degree increments (0-360)


    // Generate frcmod file
    settings.Output("Generating frcmod file.");
    Build_Frcmod(frcmod, mol, settings);

    // Test mol2/frcmod file combination in tleap
    Check_For_Missing_Parameters(settings, mol);
    std::string prmtop = settings.job_dir + "/tmp.prmtop";
    std::string restart = settings.job_dir + "/tmp.rst7";

    if (!CheckFileExists(restart) || !CheckFileExists(prmtop))
    {
        settings.Output("Attempting to regenerate frcmod file for new missing parameters.");
        Build_Frcmod(frcmod, mol, settings);
        // Test mol2/frcmod file combination in tleap
        Check_For_Missing_Parameters(settings, mol);
    }

    // if successful, return results.
    if (CheckFileExists(restart) && CheckFileExists(prmtop))
    {
        settings.Output("Successfully generated prmtop/rst7 with new parameters.");
        buffer.str("");
        buffer << "Mol2 File:   " << settings.mol2file << std::endl;
        buffer << "Frcmod File: " << settings.frcmodfile << std::endl;
        settings.Output(buffer.str());
    }
    else
    {
        settings.Error("Unable to generate prmtop/rst7 files, check leap.log output for possible issues.");
    }
    return 0;
}