#include "classes.h"

std::map<std::string, double> AtomicMasses = { 
    {"H", 1.0079 },   {"HE", 4.0026 },  {"LI", 6.941 },   {"BE", 9.0122 },  {"B", 10.811 },
    {"C", 12.0107 },  {"N", 14.0067 },  {"O", 15.9994 },  {"F", 18.9984 },  {"NE", 20.1797 },
    {"NA", 22.9897 }, {"MG", 24.305 },  {"AL", 26.9815 }, {"SI", 28.0855 }, {"P", 30.9738 },
    {"S", 32.065 },   {"CL", 35.453 },  {"AR", 39.948 },  {"K", 39.0983 },  {"CA", 40.078 },
    {"SC", 44.9559 }, {"TI", 47.867 },  {"V", 50.9415 },  {"CR", 51.9961 }, {"MN", 54.938 },
    {"FE", 55.845 },  {"CO", 58.9332 }, {"NI", 58.6934 }, {"CU", 63.546 },  {"ZN", 65.39 },
    {"GA", 69.723 },  {"GE", 72.64 },   {"AS", 74.9216 }, {"SE", 78.96 },   {"BR", 79.904 },
    {"KR", 83.8 },    {"RB", 85.4678 }, {"SR", 87.62 },   {"Y", 88.9059 },  {"ZR", 91.224 },
    {"NB", 92.9064 }, {"MO", 95.94 },   {"TC", 98 },      {"RU", 101.07 },  {"RH", 102.9055 },
    {"PD", 106.42 },  {"AG", 107.8682 },{"CD", 112.411 }, {"IN", 114.818 }, {"SN", 118.71 },
    {"SB", 121.76 },  {"TE", 127.6 },   {"I", 126.9045 }, {"XE", 131.293 }, {"CS", 132.9055 },
    {"BA", 137.327 }, {"LA", 138.9055 },{"CE", 140.116 }, {"PR", 140.9077 },{"ND", 144.24 },
    {"PM", 145 },     {"SM", 150.36 },  {"EU", 151.964 }, {"GD", 157.25 },  {"TB", 158.9253 },
    {"DY", 162.5 },   {"HO", 164.9303 },{"ER", 167.259 }, {"TM", 168.9342 },{"YB", 173.04 },
    {"LU", 174.967 }, {"HF", 178.49 },  {"TA", 180.9479 },{"W", 183.84 },   {"RE", 186.207 },
    {"OS", 190.23 },  {"IR", 192.217 }, {"PT", 195.078 }, {"AU", 196.9665 },{"HG", 200.59 },
    {"TL", 204.3833 },{"PB", 207.2 },   {"BI", 208.9804 },{"PO", 209 },     {"AT", 210 },
    {"RN", 222 },     {"FR", 223 },     {"RA", 226 },     {"AC", 227 },     {"TH", 232.0381 },
    {"PA", 231.0359 },{"U", 238.0289 }, {"NP", 237 },     {"PU", 244 },     {"AM", 243 },
    {"CM", 247 },     {"BK", 247 },     {"CF", 251 },     {"ES", 252 },     {"FM", 257 },
    {"MD", 258 },     {"NO", 259 },     {"LR", 262 },     {"RF", 261 },     {"DB", 262 },
    {"SG", 263 },     {"BH", 262 },     {"HS", 265 },     {"MT", 266 },     {"DS", 281 },
    {"RG", 272 },     {"CN", 285 },     {"UT", 284 },     {"FL", 289 },     {"UP", 288 },
    {"LV", 292 },     {"US", 0 },       {"UO", 294 } };

std::map<std::string, int> AtomicNumbers = {
    {"H", 1 },{"HE", 2 },{"LI", 3 },{"BE", 4 },{"B", 5 },
    {"C", 6 },{"N", 7 },{"O", 8 },{"F", 9 },{"NE", 10 },
    {"NA", 11 },{"MG", 12 },{"AL", 13 },{"SI", 14 },{"P", 15 },
    {"S", 16 },{"CL", 17 },{"AR", 18 },{"K", 19 },{"CA", 20 },
    {"SC", 21 },{"TI", 22 },{"V", 23 },{"CR", 24 },{"MN", 25 },
    {"FE", 26 },{"CO", 27 },{"NI", 28 },{"CU", 29 },{"ZN", 30 },
    {"GA", 31 },{"GE", 32 },{"AS", 33 },{"SE", 34 },{"BR", 35 },
    {"KR", 36 },{"RB", 37 },{"SR", 38 },{"Y", 39 },{"ZR", 40 },
    {"NB", 41 },{"MO", 42 },{"TC", 43 },{"RU", 44 },{"RH", 45 },
    {"PD", 46 },{"AG", 47 },{"CD", 48 },{"IN", 49 },{"SN", 50 },
    {"SB", 51 },{"TE", 52 },{"I", 53 },{"XE", 54 },{"CS", 55 },
    {"BA", 56 },{"LA", 57 },{"CE", 58 },{"PR", 59 },{"ND", 60 },
    {"PM", 61 },{"SM", 62 },{"EU", 63 },{"GD", 64 },{"TB", 65 },
    {"DY", 66 },{"HO", 67 },{"ER", 68 },{"TM", 69 },{"YB", 70 },
    {"LU", 71 },{"HF", 72 },{"TA", 73 },{"W", 74 },{"RE", 75 },
    {"OS", 76 },{"IR", 77 },{"PT", 78 },{"AU", 79 },{"HG", 80 },
    {"TL", 81 },{"PB", 82 },{"BI", 83 },{"PO", 84 },{"AT", 85 },
    {"RN", 86 },{"FR", 87 },{"RA", 88 },{"AC", 89 },{"TH", 90 },
    {"PA", 91 },{"U", 92 },{"NP", 93 },{"PU", 94 },{"AM", 95 },
    {"CM", 96 },{"BK", 97 },{"CF", 98 },{"ES", 99 },{"FM", 100 },
    {"MD", 101 },{"NO", 102 },{"LR", 103 },{"RF", 104 },{"DB", 105 },
    {"SG", 106 },{"BH", 107 },{"HS", 108 },{"MT", 109 },{"DS", 110 },
    {"RG", 111 },{"CN", 112 },{"UT", 113 },{"FL", 114 },{"UP", 115 },
    {"LV", 116 },{"US", 117 },{"UO", 118 } };


std::map<std::string, int> vdWRadii = { {"H", 120 },{"HE", 140 },{"LI", 182 },{"BE", 153 },{"B", 192 },
{"C", 170 },{"N", 155 },{"O", 152 },{"F", 147 },{"NE", 154 },
{"NA", 227 },{"MG", 173 },{"AL", 184 },{"SI", 210 },{"P", 180 },
{"S", 180 },{"CL", 175 },{"AR", 188 },{"K", 275 },{"CA", 231 },
{"SC", 211 },{"TI", 0 },{"V", 0 },{"CR", 0 },{"MN", 0 },
{"FE", 0 },{"CO", 0 },{"NI", 163 },{"CU", 140 },{"ZN", 139 },
{"GA", 187 },{"GE", 211 },{"AS", 185 },{"SE", 190 },{"BR", 185 },
{"KR", 202 },{"RB", 303 },{"SR", 249 },{"Y", 0 },{"ZR", 0 },
{"NB", 0 },{"MO", 0 },{"TC", 0 },{"RU", 0 },{"RH", 0 },
{"PD", 163 },{"AG", 172 },{"CD", 158 },{"IN", 193 },{"SN", 217 },
{"SB", 206 },{"TE", 206 },{"I", 198 },{"XE", 216 },{"CS", 343 },
{"BA", 268 },{"LA", 0 },{"CE", 0 },{"PR", 0 },{"ND", 0 },
{"PM", 0 },{"SM", 0 },{"EU", 0 },{"GD", 0 },{"TB", 0 },
{"DY", 0 },{"HO", 0 },{"ER", 0 },{"TM", 0 },{"YB", 0 },
{"LU", 0 },{"HF", 0 },{"TA", 0 },{"W", 0 },{"RE", 0 },
{"OS", 0 },{"IR", 0 },{"PT", 175 },{"AU", 166 },{"HG", 155 },
{"TL", 196 },{"PB", 202 },{"BI", 207 },{"PO", 197 },{"AT", 202 },
{"RN", 220 },{"FR", 348 },{"RA", 283 },{"AC", 186 },{"TH", 0 },
{"PA", 0 },{"U", 186 },{"NP", 0 },{"PU", 0 },{"AM", 0 },
{"CM", 0 },{"BK", 0 },{"CF", 0 },{"ES", 0 },{"FM", 0 },
{"MD", 0 },{"NO", 0 },{"LR", 0 },{"RF", 0 },{"DB", 0 },
{"SG", 0 },{"BH", 0 },{"HS", 0 },{"MT", 0 },{"DS", 0 },
{"RG", 0 },{"CN", 0 },{"NH", 0 },{"FL", 0 },{"MC", 0 },
{"LV", 0 },{"TS", 0 },{"OG", 0 }, };

Atom::Atom(std::string line)
{
    atom_name = trim_whitespace(line.substr(12,4));
    atom_number = atoi(line.substr(6,5).c_str());
    residue_name = line.substr(17,3);
    residue_number = atoi(line.substr(22,4).c_str());
    xx = atof(line.substr(30,8).c_str());
    yy = atof(line.substr(38,8).c_str());
    zz = atof(line.substr(46,8).c_str());
    element = trim_whitespace(line.substr(76,2));
    formal_charge = atoi(line.substr(78,2).c_str());
    mass = AtomicMasses[element];
    electron_counter = AtomicNumbers[element];
    vdw_radius = vdWRadii[element]/100.0;
    bonded_to_elements = {};
    bonded_to_indexes = {};
}

Atom::~Atom()
{

}

void Atom::set_resname(std::string resname)
{
    residue_name = resname;
}

void Atom::set_atomname(std::string atomname)
{
    atom_name = atomname;
}

void Atom::set_resp_charge(double respcharge)
{
    resp_charge = respcharge;
}

std::string Atom::print_PDB_line()
{
    std::stringstream line;
    line.str("");
    line << "ATOM  ";  //ATOM  
    line << std::setw(5) << std::right << atom_number;  // "    1" //atom index (start at 1)
    
    if (atom_name.size() == 1)
    {
        line << "  " << atom_name << "  "; // " P1 "; // atom name (space + 4 letters)
    }
    if (atom_name.size() == 2)
    {
        line << "  " << atom_name << " "; // " P1 "; // atom name (space + 4 letters)        
    }
    if (atom_name.size() >= 3)
    {
        line << " " << std::setw(4) << std::right << atom_name; // " P1 "; // atom name (space + 4 letters)        
    }
    
    line << " " << residue_name; // "UNK"; // residue name (space + 3 letters)
    line << " " << std::right << "    1"; // residue number (space + 5 numbers)
    line << "    "; // code for insertion (leave blank, I don't care.)
    line << std::fixed << std::setprecision(3) << std::setw(8) << xx; // "   8.908"; // x coordinate (8 digits, 3 decimals, includes - sign)
    line << std::fixed << std::setprecision(3) << std::setw(8) << yy; // "  -1.813"; // y coordinate (8 digits, 3 decimals, includes - sign)
    line << std::fixed << std::setprecision(3) << std::setw(8) << zz; // "  -0.503"; // z coordinate (8 digits, 3 decimals, includes - sign)
    line << "  1.00  0.00         "; // whatever, don't care.
    line << "  "; // formal atomic charge if any
    line << std::setw(2) << std::left << element; // "P "; // element, left-justified
    // line << std::setw(2) << std::left << formal_charge; // formal charge.
    // line << "  ";
    line << std::endl;
    return line.str();
}

std::string Atom::print_mol2_line()
{
    return "mol2 line\n";
}

