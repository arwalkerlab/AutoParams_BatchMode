#include "mol2.h"
// Atom Typing:
// Identify element (H, N, C, O, S, Cl, Br, F, I, P, etc.)
// Identify # of bonds vs. max total bonds (count number of bonds that include the atom index, then identify if those are single, double, triple, etc.)
// Identify characteristics (in ring, bridging, etc.)
// Select atom type from base forcefield requested. (if statements, etc.)


std::string process_S(Atom atom, Molecule mol)
{
    if (count_element_in_array(atom.bonded_to_elements,"H") > 0)
    {
        return "SH"; // If there's a hydrogen bonded to the sulfur, it's an "SH" atom type.
    }
    if (count_element_in_array(atom.bonded_to_elements,"S") > 0)
    {
        return "SS"; // If part of a disulfide bond, it's "SS" atom type.
    }
    // Otherwise, it's an "S " atom type.
    return "S ";
}

std::string process_C(Atom atom, Molecule mol)
{
    if (atom.bonded_to_elements.size() == 4)
    {
        // If it's an SP3 carbon, it's "CT"
        return "CT";
    }
    if (atom.bonded_to_elements.size() == 1)
    {
    // If it's an SP carbon, it's "CY" (bonded to nitrogen) or "CZ" (not bonded to nitrogen).
        if (std::find(atom.bonded_to_elements.begin(), atom.bonded_to_elements.end(),"N") != atom.bonded_to_elements.end())
        {
            return "CY";
        }
        return "CZ";
    }

    // Identify ring structures
    int n_rings = 0;
    bool five_member = false;
    bool six_member = false;

    for (std::vector<int> ring : mol.five_rings)
    {
        if (std::find(ring.begin(), ring.end(), atom.atom_number-1) != ring.end())
        {
            n_rings++;
            five_member = true;
        }
    }
    for (std::vector<int> ring : mol.six_rings)
    {
        if (std::find(ring.begin(), ring.end(), atom.atom_number-1) != ring.end())
        {
            n_rings++;
            six_member = true;
        }
    }

    // If 2 or more rings, it's easy.
    if (n_rings > 1)
    {
        return "CB";
    }

    // If 5-membered ring with 2 nitrogen bonds, it's 
    if (five_member)
    {
        if (count_element_in_array(atom.bonded_to_elements,"N") < 2)
        {
            // if there are fewer than two nitrogens bonded to this carbon, return CC
            return "CC";
        }
        
        // else if there are two or more nitrogens bonded to this carbon, return CQ
        return "CQ";
    }
    if (six_member)
    {
        return "CA";
    }
    // Now we're out of ring structures entirely and should be purely sp2 carbons.
    if (count_element_in_array(atom.bonded_to_elements,"O") == 0)
    {
        return "CD"; // atom in the middle of C=CD-CD=C (conjugated double-bonds)
    }
    if (count_element_in_array(atom.bonded_to_elements,"O") > 0)
    {
        return "C "; // carbonyl group
    }
    
    return "CT"; // in failure, return an sp3 carbon...
}

std::string process_H(Atom atom, Molecule mol)
{
    if (count_element_in_array(atom.bonded_to_elements,"N") == 1)
        return "H ";
    if (count_element_in_array(atom.bonded_to_elements,"O") == 1)
        return "HO";
    if (count_element_in_array(atom.bonded_to_elements,"S") == 1)
        return "HS";
    if (count_element_in_array(atom.bonded_to_elements,"C") == 1)
    {
        // identify carbon type.
        int bonded_atom = atom.bonded_to_indexes[0];
        std::string bonded_atom_type = process_C(mol.atoms[bonded_atom],mol);
        std::map <std::string,std::string> hydrogen_map = {{"CY","HZ"},{"CZ","HZ"},{"CQ","H4"},{"CA","HA"},{"CC","H5"},{"CD","HC"}};
        if (hydrogen_map.count(bonded_atom_type))
        {
            return hydrogen_map[bonded_atom_type];
        }
        return "HC";
    }

    return "HC";    
}

std::string process_N(Atom atom, Molecule mol)
{
    int n_bonds = atom.bonded_to_indexes.size();
    if (n_bonds == 1)
    {
        return "NY";
    }
    if (n_bonds == 3)
    {
        return "N ";
    }
    if (n_bonds == 4)
    {
        return "N3";
    }

    // Identify ring structures
    int n_rings = 0;
    bool five_member = false;
    bool six_member = false;

    for (std::vector<int> ring : mol.five_rings)
    {
        if (std::find(ring.begin(), ring.end(), atom.atom_number-1) != ring.end())
        {
            n_rings++;
            five_member = true;
        }
    }
    for (std::vector<int> ring : mol.six_rings)
    {
        if (std::find(ring.begin(), ring.end(), atom.atom_number-1) != ring.end())
        {
            n_rings++;
            six_member = true;
        }
    }
    if (n_rings > 1)
    {
        return "N*";
    }
    if (n_rings == 1)
    {
        if (five_member)
        {
            if (count_element_in_array(atom.bonded_to_elements,"H") > 0)
            {
                return "NA";
            }
            return "NB";
        }
        if (six_member)
        {
            return "NC";
        }
    }
    if (count_element_in_array(atom.bonded_to_elements,"C") > 0)
    {
        for (int idx : atom.bonded_to_indexes)
        {
            if (mol.atoms[idx].element == "C")
            {
                if (count_element_in_array(mol.atoms[idx].bonded_to_elements,"O") > 0)
                {
                    return "N ";
                }
            }
        }
    }
    return "N*";
}

std::string process_O(Atom atom, Molecule mol)
{
    int n_bonds = atom.bonded_to_elements.size();
    if (count_element_in_array(atom.bonded_to_elements,"P") > 0)
    {
        if (n_bonds == 1)
        {
            return "OP";
        }
        return "O2";
    }
    if (count_element_in_array(atom.bonded_to_elements,"H") > 0)
    {
        return "OH";
    }
    if (count_element_in_array(atom.bonded_to_elements,"C") > 0)
    {
        if (n_bonds == 1)
        {
            return "O ";
        }
        return "OS";
    }
    return "OH";
}


void AtomTyping(Molecule &mol)
{
    for (Atom atom : mol.atoms)
    {
        if (atom.element == "S")
        {
            atom.atom_type = process_S(atom,mol);
            continue;
        }
        if (atom.element == "C")
        {
            atom.atom_type = process_C(atom,mol);
            continue;
        }
        if (atom.element == "H")
        {
            atom.atom_type = process_H(atom,mol);
            continue;
        }
        if (atom.element == "N")
        {
            atom.atom_type = process_N(atom,mol);
            continue;
        }
        if (atom.element == "H")
        {
            atom.atom_type = process_H(atom,mol);
            continue;
        }
        std::stringstream buf;
        buf.str("");
        buf << atom.element;
        if (atom.element.size() < 2)
            buf << " ";
        atom.atom_type = buf.str();
    }
}
