#include "globals.h"
#include <iostream>

namespace LoopCGAL
{
    bool verbose = false; // Definition of the verbose flag

    void set_verbose(bool value)
    {
        verbose = value;
        std::cout << "Verbose flag set to: " << verbose << std::endl;
    }
}