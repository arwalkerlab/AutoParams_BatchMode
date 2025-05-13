#include "settings.h"

void Settings::Output(std::string output)
{
    std::ofstream outfile_stream;
    if (!CheckFileExists(outfile))
    {
        outfile_stream.open(outfile,std::ios::out);
        outfile_stream << FormattedTimeStamp();
        outfile_stream << "Initializing AutoParams Output..." << std::endl;
    }
    else
    {
        outfile_stream.open(outfile,std::ios::app);
    }    
    outfile_stream << output;
    if (output.back() != '\n')
    {
        outfile_stream << "\n";
    }
    outfile_stream.close();
}

void Settings::Error(std::string error)
{
    std::ofstream errfile_stream;    
    if (!CheckFileExists(outfile))
    {
        errfile_stream.open(errfile,std::ios::out);
        errfile_stream << FormattedTimeStamp();
        errfile_stream << "Initialize AutoParams Error Logging." << std::endl;
    }
    else
    {
        errfile_stream.open(errfile,std::ios::app);
    }    
    errfile_stream << error;
    if (error.back() != '\n')
    {
        errfile_stream << "\n";
    }
    errfile_stream.close();
}

void Settings::CheckPrograms()
{
    Output("Locating External Programs...\n");
    TC_EXISTS = CheckProgramExists("terachem");
    PSI4_EXISTS = CheckProgramExists("psi4");
    AIMNET_EXISTS = CheckProgramExists("aimnet2");
    TLEAP_EXISTS = CheckProgramExists("tleap");
    PARMED_EXISTS = CheckProgramExists("parmed");
    CPPTRAJ_EXISTS = CheckProgramExists("cpptraj");
    ANTECHAMBER_EXISTS = CheckProgramExists("antechamber");

    if (TC_EXISTS) { Output("TeraChem located."); }
    if (PSI4_EXISTS) { Output("Psi4 located."); }
    if (AIMNET_EXISTS) { Output("AimNET2 located."); }
    if (TLEAP_EXISTS) { Output("Tleap located."); }
    if (PARMED_EXISTS) { Output("Parmed located."); }
    if (CPPTRAJ_EXISTS) { Output("CPPTRAJ located."); }
    if (ANTECHAMBER_EXISTS) { Output("Antechamber located."); }
    Output("\n");
}

void Settings::SetUpLogFiles()
{
    // Set up output and error files, ensuring numerical naming and matched pairs with job folder.
    // If either the output or error files exist at a given number, the number is advanced for both to maintain matching pairs.
    int i=1;
    outfile = fs::absolute("autoparams.0000.out");
    errfile = fs::absolute("autoparams.0000.err");
    job_dir = fs::absolute("autoparams.0000/");
    
    if (CheckFileExists(outfile) | CheckFileExists(errfile))
    {
        do
        {
            std::stringstream buffer;
            buffer.str("");
            buffer << "autoparams." << std::setw(4) << std::setfill('0') << i << ".out";
            outfile = buffer.str();
            buffer.str("");
            buffer << "autoparams." << std::setw(4) << std::setfill('0') << i << ".err";
            errfile = buffer.str();
            buffer.str("");
            buffer << "autoparams." << std::setw(4) << std::setfill('0') << i << "/";
            job_dir = buffer.str();
            i++;
        } while (CheckFileExists(outfile) | CheckFileExists(errfile));
    }
    // Create job folder (scratch, etc.)
    if (!fs::exists(job_dir)) 
    {
        fs::create_directory(job_dir);
    }
}

void Settings::parse_command_line(int argc,char **argv)
{
    Output("Parsing command line arguments...\n");
    mol_charge = 0;
    mol_spin = 1;
    forcefield = "PROTEIN";
    inputfile = "";
    OPTIMIZE_FIRST = false;
    head_atom_name = "0";
    tail_atom_name = "0";
    dummy_atom_names = {};


    /*
    valid command line flags:
        -c,--charge = formal charge of molecular fragment
        -s,--spin   = multiplicity of molecular fragment
        -i,--input  = coordinates/structure file
        -f,--forcefield = forcefield choice
        -o,--optimize = run TeraChem gas-phase optimization first
        -h,--head = head atom name
        -t,--tail = tail atom name
    */
   
    int i = 0;
    while (i < argc)
    {
        if (((std::string)argv[i] == "-c") || ((std::string)argv[i] == "--charge" ))
        {
            mol_charge = atoi(argv[i+1]);
            i+=2;
            continue;
        }
        if (((std::string)argv[i] == "-s") || ((std::string)argv[i] == "--spin") )
        {
            mol_spin = atoi(argv[i+1]);
            i+=2;
            continue;
        }
        if (((std::string)argv[i] == "-f") || ((std::string)argv[i] == "--forcefield") )
        {
            forcefield = argv[i+1];
            transform(forcefield.begin(), forcefield.end(), forcefield.begin(),::toupper);
            i+=2;
            continue;
        }
        if (((std::string)argv[i] == "-i") || ((std::string)argv[i] == "--input") )
        {
            inputfile = argv[i+1];
            i+=2;
            continue;
        }
        if (((std::string)argv[i] == "-o") || ((std::string)argv[i] == "--optimize") )
        {
            if (TC_EXISTS)
            {
                Output("Gas phase QM optimization will be performed before RESP fitting.\n\tThis may take a while.");
                OPTIMIZE_FIRST = true;
            }
            else
            {
                Error("Optimization requested but TeraChem is not available.  Skipping optimization.");
                OPTIMIZE_FIRST = false;
            }
            i++;
            continue;
        }
        if (((std::string)argv[i] == "-h") || ((std::string)argv[i] == "--head" ))
        {
            head_atom_name = (std::string)argv[i+1];
            i+=2;
            continue;
        }
        if (((std::string)argv[i] == "-t") || ((std::string)argv[i] == "--tail" ))
        {
            tail_atom_name = (std::string)argv[i+1];
            i+=2;
            continue;
        }
        if (((std::string)argv[i] == "-d") || ((std::string)argv[i] == "--dummy" ))
        {
            dummy_atom_names.push_back((std::string)argv[i+1]);
            i+=2;
            continue;
        }
                
        i++;
    }
    if (inputfile == "")
    {
        std::stringstream buffer;
        buffer.str("");
        buffer << "Input coordinates not provided.  Please specify input coordinate file." << std::endl;
        Error(buffer.str());
        return;
    }
    if (!CheckFileExists(inputfile))
    {
        std::stringstream buffer;
        buffer.str("");
        buffer << "Input coordinates not found.  Please check the filename: " << inputfile << std::endl;
        Error(buffer.str());
        return;
    }
    mol2file = inputfile.substr(0,inputfile.find_last_of('.')) + ".mol2";
    mol2file = job_dir + mol2file;
    frcmodfile = job_dir + inputfile.substr(0,inputfile.find_last_of('.')) + ".frcmod";

    std::stringstream buffer;
    buffer.str("");
    buffer << "Input coordinates: " << inputfile << std::endl;
    buffer << "Molecular charge:  " << mol_charge << std::endl;
    buffer << "Molecular spin:    " << mol_spin << std::endl;
    buffer << "Base Forcefield:   " << forcefield << std::endl;
    buffer << "Mol2 filename:     " << string_split(mol2file,'/').back() << std::endl;
    buffer << "Frcmod filename:   " << string_split(frcmodfile,'/').back() << std::endl;
    buffer << "\n";
    Output(buffer.str());
    buffer.str("");
}

Settings::Settings()
{
    SetUpLogFiles();
    CheckPrograms();
    return;
}

Settings::~Settings()
{

}