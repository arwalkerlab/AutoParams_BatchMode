#ifndef MOL2_H
#define MOL2_H

#include "utilities.h"
#include "classes.h"

class Mol2File
{
    private:
        std::vector <Atom> atoms;
    public:
        void WriteMol2();
        Mol2File(Molecule mol);
        ~Mol2File();
};

#endif