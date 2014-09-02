#ifndef gmx2pqr_read_parm_hpp
#define gmx2pqr_read_parm_hpp
#include "my_structs.hpp"
#include <fstream>
#include <sstream>

void read_DAT( const char *filename, std::vector<gmx2amb> &dat, const bool &bVerbose);

#endif