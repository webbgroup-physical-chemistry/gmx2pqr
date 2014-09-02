//
//  main.cpp
//  apbs_numeric_rxn_field
//
//  Created by Andrew Ritchie on 5/19/14.
//  Copyright (c) 2014 Andrew Ritchie. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include <cmath>

void print_selection( const std::vector<int> &selection, std::string name );

double parse_selection( const std::vector<float> &rxnpot, const std::vector<std::vector<float> > &coord, const std::vector<int> &selection ); 

int main(int argc, const char * argv[])
{
    const char *pqr, *solv, *ref;
    if (argc == 4)
    {
        pqr = argv[1];
        solv = argv[2];
        ref = argv[3];
    }
    else
    {
        std::cerr << "\nUsage: <pqr> <apbs atom potentials in solvent> <apbs atom potentials in vacuum>\n" << std::endl;
        std::exit(1);
    }
    
    /* Read dummy atom coords from PQR file */
    std::string pqrline,throwAway,atom,resname;
    std::vector<float> coords(3,0);

    int n = 0;
    std::vector<std::vector<float> > coordArray;
    std::vector<int> DumLineNumbers;

    std::vector<int> midzN;
    std::vector<int> midyN;
    std::vector<int> midxN;

    std::vector<int> a1zN;
    std::vector<int> a1yN;
    std::vector<int> a1xN;

    std::vector<int> a2zN;
    std::vector<int> a2yN;
    std::vector<int> a2xN;

    std::vector<int> pN;

    std::ifstream file (pqr);
    if (file.is_open())
    {
        while (file.good())
        {
            getline(file, pqrline);
            if (not pqrline.empty())
            {
                std::stringstream linestream(pqrline);
                linestream >> throwAway >> throwAway >> atom >> resname >> throwAway >> coords[0] >> coords[1] >> coords[2];
                if (std::strncmp(resname.c_str(),"DUM",4) == 0)
                {
                    coordArray.push_back(coords);
                    DumLineNumbers.push_back(n);
                    if (std::strncmp(atom.c_str(),"MIDz",5) == 0)
                    {
                        midzN.push_back(n);
                    }
                    else
                    if (std::strncmp(atom.c_str(),"MIDy",5) == 0)
                    {
                        midyN.push_back(n);
                    }
                    else
                    if (std::strncmp(atom.c_str(),"MIDx",5) == 0)
                    {
                        midxN.push_back(n);
                    }
                    else
                    if (std::strncmp(atom.c_str(),"CDz",4) == 0)
                    {
                        a1zN.push_back(n);
                    }
                    else
                    if (std::strncmp(atom.c_str(),"CDy",4) == 0)
                    {
                        a1yN.push_back(n);
                    }
                    else
                    if (std::strncmp(atom.c_str(),"CDx",4) == 0)
                    {
                        a1xN.push_back(n);
                    }
                    else
                    if (std::strncmp(atom.c_str(),"NEz",4) == 0)
                    {
                        a2zN.push_back(n);
                    }
                    else
                    if (std::strncmp(atom.c_str(),"NEy",4) == 0)
                    {
                        a2yN.push_back(n);
                    }
                    else
                    if (std::strncmp(atom.c_str(),"NEx",4) == 0)
                    {
                        a2xN.push_back(n);
                    }
                    else
                    if (std::strncmp(atom.c_str(),"p",1) == 0)
                    {
                        pN.push_back(n);
                    }
                }
            }
            n++;
        }
    }
    int nDum = (int)coordArray.size();
    
    /* Read potentials from APBS atom potentials in solvent */
    std::string solvline;
    std::vector<float> solvpot(nDum,0);
    std::ifstream solvfile (solv);
    int i = 0, dumN = 0;
    if (solvfile.is_open())
    {
        while (solvfile.good())
        {
            getline(solvfile,solvline);
            if (not solvline.empty())
            {
                std::stringstream solvstream(solvline);
                solvstream >> throwAway;
                if (std::strncmp(throwAway.c_str(),"#",1) != 0)
                {
                    if (i == DumLineNumbers[dumN])
                    {
                        std::stringstream pot(solvline);
                        pot >> solvpot[dumN];
                        dumN++;
                    }
                    i++;
                }
            }
        }
    }
    
    /* Read potentials from APBS atom potentials in vacuum */
    std::string refline;
    std::vector<float> refpot(nDum,0);
    std::ifstream reffile (ref);
    i = 0;
    dumN = 0;
    if (reffile.is_open())
    {
        while (reffile.good())
        {
            getline(reffile,refline);
            if (not refline.empty())
            {
                std::stringstream refstream(refline);
                refstream >> throwAway;
                if (std::strncmp(throwAway.c_str(),"#",1) != 0)
                {
                    if (i == DumLineNumbers[dumN])
                    {
                        std::stringstream pot(refline);
                        pot >> refpot[dumN];
                        dumN++;
                    }
                    i++;
                }
            }
        }
    }
    
    /* Rxn potential is (potential in solvent) - (potential in vacuum) */
    std::vector<float> rxnpot(nDum,0);
    for (int i=0; i<nDum; i++)
    {
        rxnpot[i] = solvpot[i] - refpot[i];
    }

    /* print out some info to stderr */
    print_selection(midzN, "MIDz");
    print_selection(midyN, "MIDy");
    print_selection(midxN, "MIDx");
    
    print_selection(a1zN, "CDz");
    print_selection(a1yN, "CDy");
    print_selection(a1xN, "CDx");
    
    print_selection(a2zN, "NEz");
    print_selection(a2yN, "NEy");
    print_selection(a2xN, "NEx");

    int psize = pN.size();
    if ( psize > 0 ) {
        for (int i=0 ; i<psize; i++) {
            std::cerr << " p" << i;
        }
    }

    std::cerr << std::endl;

    parse_selection(rxnpot, coordArray, midzN);
    parse_selection(rxnpot, coordArray, midyN);
    parse_selection(rxnpot, coordArray, midxN);
    parse_selection(rxnpot, coordArray, a1zN);
    parse_selection(rxnpot, coordArray, a1yN);
    parse_selection(rxnpot, coordArray, a1xN);
    parse_selection(rxnpot, coordArray, a2zN);
    parse_selection(rxnpot, coordArray, a2yN);
    parse_selection(rxnpot, coordArray, a2xN);

    if ( psize > 0 ) {
        for (int i=0 ; i<psize; i++) {
            std::cerr << " " << rxnpot[pN[i]];
        }
    }

    std::cout << std::endl;

    return 0;
}

void print_selection( const std::vector<int> &selection, std::string name )
{
    int n = selection.size();
    if (n>0) 
    {
        std::cerr << " " << name;
    }
}

double parse_selection( const std::vector<float> &rxnpot, const std::vector<std::vector<float> > &coord, const std::vector<int> &selection ) 
{
    int n = selection.size();
    if (n>0) 
    {
        /* 
         *  Assuming equally spaced dummy atoms, the middle 2 values are 
         * the closest to the midpoint, and will therefore have the 
         * smallest error in the first derivative. 
         */
        int n0 = selection[0], n1 = selection[n-1];
        while (n0+1 < n1 && n0 < n1-1 && n0+1 < n1-1 && n1+1 > n0-1)
        {
            n0++;
            n1--;
        }
        /* Field is -dV/dR */
        double dx2 = (coord[n0][0] - coord[n1][0]) * 
                     (coord[n0][0] - coord[n1][0]);
        double dy2 = (coord[n0][1] - coord[n1][1]) *
                     (coord[n0][1] - coord[n1][1]);
        double dz2 = (coord[n0][2] - coord[n1][2]) *
                     (coord[n0][2] - coord[n1][2]);
        double dR = sqrt(dx2+dy2+dz2);
        double dV = rxnpot[n1] - rxnpot[n0];

        std::cout << " " << -dV/dR;

        return -dV/dR;
    }
    return 0;
}
