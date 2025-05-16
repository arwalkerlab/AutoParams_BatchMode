#include "mol2.h"


Mol2File::Mol2File(Molecule mol)    //Constructor
{  
    for (Atom atom : mol.atoms)
    {
        std::cout << std::endl;
    }

}

Mol2File::~Mol2File()   //Destructor
{

}

// Atom Typing:
// Identify element (H, N, C, O, S, Cl, Br, F, I, P, etc.)
// Identify # of bonds vs. max total bonds (count number of bonds that include the atom index, then identify if those are single, double, triple, etc.)
// Identify characteristics (in ring, bridging, etc.)
// Select atom type from base forcefield requested. (if statements, etc.)

std::vector<std::string> get_bonded_elements(Atom atom, Molecule mol)
{
    std::vector <std::string> bonded_atoms = {};
    for (std::vector<int> bond : mol.bonds)
    {
        if (bond[0] == atom.atom_number - 1)
        {
            bonded_atoms.push_back(mol.atoms[bond[1]].element);
        }
        else if (bond[1] == atom.atom_number - 1)
        {
            bonded_atoms.push_back(mol.atoms[bond[0]].element);
        }
        else
        {
            continue;
        }
    }
    return bonded_atoms;
}

std::string process_S(Atom atom,Molecule mol)
{
    std::vector <std::string> bonded_atoms = get_bonded_elements(atom, mol);
    if (std::find(bonded_atoms.begin(), bonded_atoms.end(), "H") != bonded_atoms.end())
    {
        // If there's a hydrogen bonded to the sulfur, it's an "SH" atom type.
        return "SH";
    }
    // Otherwise, it's an "S " atom type.
    return "S ";
}

std::string process_C(Atom atom,Molecule mol)
{
    std::vector <std::string> bonded_atoms = get_bonded_elements(atom, mol);
    if (bonded_atoms.size() == 4)
    {
        // If it's an SP3 carbon, it's "CT"
        return "CT";
    }
    if (bonded_atoms.size() == 1)
    {
    // If it's an SP carbon, it's "CY" (bonded to nitrogen) or "CZ" (not bonded to nitrogen).
        if (std::find(bonded_atoms.begin(), bonded_atoms.end(),"N") != bonded_atoms.end())
        {
            return "CY";
        }
        return "CZ";
    }

    // // Identify ring structures

    // If 2 or more rings, it's easy.

    // If 5-membered ring with 2 nitrogen bonds, it's 
/*

def process_C(atom,ring_groups,bridging_atoms):
    n_bonds = len(atom.bonds)
    #### SP3 CARBONS
    if n_bonds == 4:
        bond_atoms = get_bonded_elements(atom)
        if all(["H" in bond_atoms,"N" in bond_atoms,"C" in bond_atoms]):
            # CX 12.01         0.360               protein C-alpha (new to ff10)
            return "CX"
        # CT 12.01         0.878               sp3 aliphatic C
        return "CT"
    
    #### SP CARBONS
    if n_bonds == 2:
        bond_atoms = get_bonded_elements(atom)
        # CY 12.01         0.360               nitrile C (Howard et al.JCC,16,243,1995)
        if "N" in bond_atoms:
            return "CY"
        # CZ 12.01         0.360               sp C (Howard et al.JCC,16,243,1995)
        else:
            return "CZ"

    #### RING CARBONS
    ring_sizes=[]
    for group in ring_groups:
        if atom.idx in group:
            ring_sizes.append(len(group))
            
    ## Bridging carbons
    if len(ring_sizes) > 1:
        # CB 12.01         0.360               sp2 aromatic C, 5&6 membered ring junction
        return "CB"
    ## All Other Rings
    elif len(ring_sizes) == 1:
        ## 5-member rings
        if ring_sizes[0] == 5:
            bond_atoms = get_bonded_elements(atom)
            n_nitrogens = sum([bool(x == "N" for x in bond_atoms)])
            if n_nitrogens < 2:
                # CC 12.01         0.360               sp2 aromatic C, 5 memb. ring HIS
                return "CC"
            else:    
                # CQ 12.01         0.360               sp2 C in 5 mem.ring of purines between 2 N
                return "CQ"

        ## 6-membered rings
        elif ring_sizes[0] == 6:
            # CA 12.01         0.360               sp2 C pure aromatic (benzene)
            return "CA"
    #### SP2 CARBONS
    elif len(ring_sizes) == 0:
        bond_atoms = get_bonded_elements(atom)
        if "O" not in bond_atoms:
            # CD 12.01         0.360               sp2 C atom in the middle of: C=CD-CD=C
            return "CD"
        else:
            # C  12.01         0.616  !            sp2 C carbonyl group
            return "C "
    # If all else fails, return generic SP3 carbon...
    return "CT"

*/
}

/*
def process_H(atom,ring_groups,bridging_atoms):
    bond_atoms = get_bonded_elements(atom)
    if bond_atoms[0] == "N":
        # H  1.008         0.161               H bonded to nitrogen atoms
        return "H "
    if bond_atoms[0] == "O":
        # HO 1.008         0.135               hydroxyl group
        return "HO"
    if bond_atoms[0] == "S":
        # HS 1.008         0.135               hydrogen bonded to sulphur (pol?)
        return "HO"
    if bond_atoms[0] == "C":
        if atom.bonds[0].atom1.element == "C":
            carb_at_type = process_C(atom.bonds[0].atom1,ring_groups,bridging_atoms)
        else:
            carb_at_type = process_C(atom.bonds[0].atom2,ring_groups,bridging_atoms)
        if carb_at_type in ["CY","CZ"]:
            # HZ 1.008         0.161               H bond sp C (Howard et al.JCC,16,243,1995)
            return "HZ"
        elif carb_at_type == "CQ":
            # H4 1.008         0.167               H arom. bond. to C with 1 electrwd. group
            return "H4"
        elif carb_at_type == "CA":
            # HA 1.008         0.167               H arom. bond. to C without elctrwd. groups    
            return "HA"
        elif carb_at_type == "CC":
            # H5 1.008         0.167               H arom.at C with 2 elctrwd. gr,+HCOO group
            return "H5"
            
        elif carb_at_type == "CD":
            # HC 1.008         0.135               H aliph. bond. to C without electrwd.group
            return "HC"
        return "HC"
    return "HC"

def process_N(atom,ring_groups,bridging_atoms):
    bond_atoms = get_bonded_elements(atom)
    n_bonds = len(bond_atoms)
    if n_bonds == 1:
        # NY 14.01         0.530               nitrile N (Howard et al.JCC,16,243,1995)
        return "NY"
    elif n_bonds == 3:
        # NT 14.01         0.530               sp3 N for amino groups amino groups
        return "N "
    elif n_bonds == 4:
        # N3 14.01         0.530               sp3 N for charged amino groups (Lys, etc)
        return "N3"
    else:
        ring_sizes=[]
        for group in ring_groups:
            if atom.idx in group:
                ring_sizes.append(len(group))
        if len(ring_sizes) == 1:
            if ring_sizes[0] == 5:
                if "H" in bond_atoms:
                    # NA 14.01         0.530               sp2 N in 5 memb.ring w/H atom (HIS)
                    return "NA"
                # NB 14.01         0.530               sp2 N in 5 memb.ring w/LP (HIS,ADE,GUA)
                return "NB"
            elif ring_sizes[0] == 6:
                # NC 14.01         0.530               sp2 N in 6 memb.ring w/LP (ADE,GUA)
                return "NC"
        elif len(ring_sizes) == 0:
            if "C" in bond_atoms:
                for bond in atom.bonds:
                    if bond.atom1.element_name == "C":
                        c_bonded = get_bonded_elements(bond.atom1)
                        if "O" in c_bonded:
                            # N  14.01         0.530               sp2 nitrogen in amide groups
                            return "N "
        # N* 14.01         0.530               sp2 N
        return "N*"
            
def process_O(atom):
    n_bonds = len(atom.bonds)
    bond_atoms = get_bonded_elements(atom)
    if "P" in bond_atoms:
        if n_bonds == 1:
            # OP 16.00         0.465               2- phosphate oxygen
            return "OP"
        # O2 16.00         0.434               carboxyl and phosphate group oxygen
        return "O2"
    if "H" in bond_atoms:
        # OH 16.00         0.465               oxygen in hydroxyl group
        return "OH"
    if "C" in bond_atoms:
        if n_bonds == 1:
            # O  16.00         0.434               carbonyl group oxygen
            return "O "
        # OS 16.00         0.465               ether and ester oxygen
        return "OS"
*/