#include "classes.h"
#include "settings.h"

Frcmod_File::Frcmod_File(std::string filename)
{
    std::vector <std::string> params_added = {};
    frcmod_filename = filename;
    masses_section = "MASS\n";
    bonds_section = "BOND\n";
    angles_section = "ANGLE\n";
    dihedrals_section = "DIHE\n";
    torsions_section = "IMPROPER\n";
    vdw_section = "NONBON\n";

}

Frcmod_File::~Frcmod_File()
{
    if (frcmod_write.is_open())
    {
        frcmod_write.close();
    }
}

void Frcmod_File::Add_Mass(std::string line)
{
    std::string _tmpline = line;
    if (line.back() != '\n')
    {
        _tmpline += '\n';
    }
    if (count(params_added.begin(), params_added.end(), _tmpline) == 0)
    {
        params_added.push_back(_tmpline);
        masses_section += _tmpline;
    }
    
}

void Frcmod_File::Add_Bond(std::string line)
{
    std::string _tmpline = line;
    if (line.back() != '\n')
    {
        _tmpline += '\n';
    }
    if (count(params_added.begin(), params_added.end(), _tmpline) == 0)
    {
        params_added.push_back(_tmpline);
        bonds_section += _tmpline;
    }
}

void Frcmod_File::Add_Angle(std::string line)
{
    std::string _tmpline = line;
    if (line.back() != '\n')
    {
        _tmpline += '\n';
    }
    if (count(params_added.begin(), params_added.end(), _tmpline) == 0)
    {
        params_added.push_back(_tmpline);
        angles_section += _tmpline;
    }
    
}

void Frcmod_File::Add_Dihedral(std::string line)
{
    std::string _tmpline = line;
    if (line.back() != '\n')
    {
        _tmpline += '\n';
    }
    if (count(params_added.begin(), params_added.end(), _tmpline) == 0)
    {
        params_added.push_back(_tmpline);
        dihedrals_section += _tmpline;
    }
}

void Frcmod_File::Add_Torsion(std::string line)
{
    std::string _tmpline = line;
    if (line.back() != '\n')
    {
        _tmpline += '\n';
    }
    if (count(params_added.begin(), params_added.end(), _tmpline) == 0)
    {
        params_added.push_back(_tmpline);
        torsions_section += _tmpline;
    }
}

void Frcmod_File::Add_vdW(std::string line)
{
    std::string _tmpline = line;
    if (line.back() != '\n')
    {
        _tmpline += '\n';
    }
    if (count(params_added.begin(), params_added.end(), _tmpline) == 0)
    {
        params_added.push_back(_tmpline);
        vdw_section += _tmpline;
    }
}


void Frcmod_File::Write_Frcmod_File(Settings settings)
{
    // Prepare buffer for data entry.
    std::stringstream buffer;
    buffer.str("");

    frcmod_write.open(frcmod_filename,std::ios::out);
    if (!frcmod_write.is_open())
    {
        settings.Error("Unable to open FRCMOD file to write.");
        return;
    }
    std::map<std::string,std::string> leaprc_dict={{"DNA","leaprc.DNA.OL15"},{"RNA","leaprc.RNA.OL3"},{"PROTEIN","leaprc.protein.ff14SB"},{"CARBOHYDRATE","leaprc.GLYCAM_06j-1"}};
    // Write title line
    buffer << settings.frcmodfile.substr(settings.frcmodfile.find_last_of('/')+1,settings.frcmodfile.find_last_of('.') - settings.frcmodfile.find_last_of('/')-1) << ", use with " << leaprc_dict[settings.forcefield] << "." << std::endl;
    frcmod_write << buffer.str();
    buffer.str("");

    // write MASS
    frcmod_write << masses_section << std::endl;

    // write BOND
    frcmod_write << bonds_section << std::endl;

    // write ANGLE
    frcmod_write << angles_section << std::endl;

    // write DIHE
    frcmod_write << dihedrals_section << std::endl;
    
    // write IMPROPER
    frcmod_write << torsions_section << std::endl;
    
    // write NONBON
    frcmod_write << vdw_section << std::endl << std::endl;


    frcmod_write.close();
}
