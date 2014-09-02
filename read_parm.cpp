#include "read_parm.hpp"

void read_DAT( const char *filename, std::vector<gmx2amb> &dat, const bool &bVerbose)
{
    fprintf(stderr,"\nWill read reference charges and radii from %s\n",filename);
    
    std::string line;
    std::ifstream file(filename);
    if (file.is_open()){
        while (file.good()){
            getline(file, line);
            if (not line.empty() && line.substr(0,1) != ";" && line.substr(0,1) != "#") {
                gmx2amb eachline;
                bool add_to_library = true;
                std::stringstream linestream(line);
                linestream >> eachline.resname >> eachline.gmxname >> eachline.charge >> eachline.radius >> eachline.ambername;
                /* GLY H0 in gromacs is GLY h1 in AMBER.DAT */
                if (strncmp(eachline.resname,"GLY",strlen(eachline.resname)+1) == 0 &&
                    strncmp(eachline.ambername,"H1",strlen(eachline.ambername)+1) == 0) {
                    memcpy(eachline.ambername,"H0",strlen("HO")+1);
                }
                /* Check to see if the entry already exists in dat */
                for (int i=0; (int) i<dat.size(); i++){
                    if (strncmp(eachline.ambername, dat[i].ambername, strlen(eachline.ambername)+1) == 0){
                        add_to_library = false;
                        break;
                    }
                }
                if (add_to_library){
                    dat.push_back(eachline);
                }
            }
        }
    }
    else if (!file)
    {
        fprintf(stderr,"\nError reading %s\n", filename);
        std::exit(1);
    }
    file.close();
   
    if (bVerbose) {
        for (int i=0; i<(int) dat.size(); i++) {
            fprintf(stderr,"\n%i out of %i: ambername = %s, radius = %.4f",i,(int)dat.size(), dat[i].ambername, dat[i].radius);
        }
        fprintf(stderr,"\n");
    }
    return;
}