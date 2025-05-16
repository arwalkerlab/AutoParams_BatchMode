#ifndef MOL2_H
#define MOL2_H

#include "utilities.h"
#include "classes.h"

std::string process_C(Atom atom, Molecule mol);
std::string process_N(Atom atom, Molecule mol);
std::string process_S(Atom atom, Molecule mol);
std::string process_O(Atom atom, Molecule mol);
std::string process_H(Atom atom, Molecule mol);
void AtomTyping(Molecule &mol);

#endif