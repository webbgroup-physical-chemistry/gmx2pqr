/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009,2010,2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "my_structs.hpp"
#include "read_parm.hpp"
#include "linalg.hpp"
#include "electrostatics.hpp"
#ifdef __cplusplus
extern "C"
#endif
#include <gromacs/copyrite.h>
#include <gromacs/filenm.h>
#include <gromacs/macros.h>
#include <gromacs/pbc.h>
#include <gromacs/statutil.h>
#include <gromacs/xvgr.h>
#include <gromacs/gmx_fatal.h>
#include <gromacs/trajana.h>


/*! \brief
 * Template analysis data structure.
 */
typedef struct
{
    gmx_ana_selection_t     *refsel;
    FILE                    *fp;
    const char              *pqr;
    int                     framen;
    int                     a1, a2;
    real                    delta, ring_dist;
    double                  rpdie;
    std::vector<int>        site_ndx, exclude_ndx;
    bool                    bVerbose;
    std::vector<gmx2amb>    amber;
} t_analysisdata;

typedef struct
{
    gmx_bool                bDump;
    gmx_bool                bFracNorm;
    const char             *routt;
    int                    *size;
    FILE                   *sfp;
    FILE                   *cfp;
    FILE                   *ifp;
    t_blocka               *block;
    char                  **gnames;
    FILE                   *mfp;
    gmx_ana_indexmap_t     *mmap;
} t_dsdata;

/*! \brief
 * Function that does the analysis for a single frame.
 *
 * It is called once for each frame.
 */
static int analyze_frame(t_topology *top, t_trxframe *fr, t_pbc *pbc,
              int nr, gmx_ana_selection_t *sel[], void *data)
{
    t_analysisdata      *d = (t_analysisdata *)data;
    t_electrostatics    dat;
    
    /* Print the current frame number to the xvg file */
    if (d->fp) {
        fprintf(d->fp, "%10i ", d->framen);
    }
    
    /* Index the pqr file to the frame number */
    std::string pqrname(d->pqr);
    pqrname = pqrname.substr(0,pqrname.size()-4);
    std::stringstream pqr;
    pqr << pqrname << d->framen << ".pqr";
    FILE *pqrout = ffopen(pqr.str().c_str(), "w");
    if ( pqrout == NULL ) {
        gmx_fatal(FARGS, "\nFailed to open output file, '%s'\n",pqrout);
    }
    
    /* Get the bond vector and bond midpoint */
    rvec bondvector, midpoint;
    rvec_sub(fr->x[d->a2],fr->x[d->a1],bondvector);
    double bondlength = pow( iprod(bondvector,bondvector), 0.5);
    for (int i=0; i<3; i++){
        midpoint[i] = fr->x[d->a1][i] + 0.5 * bondvector[i];
    }
    
    /* 
     * Defined axes based on the nitrile bond vector
     * The z-axis is the bond axis
     * The y-axis is from solving y0 = -(z1*y1+z2*y2)/z0, where (y1,y2)=(z2,-z0), then normalizing
     * The x-axis is the cross product of the y-axis with the z-axis
     */
    double zvec[3] = { bondvector[0] / bondlength, bondvector[1] / bondlength, bondvector[2] / bondlength };
    double yvec[3];
    if (zvec[1] != zvec[2] ) {
        yvec[0] = 1 ; yvec[1] = zvec[2]; yvec[2] = -zvec[1];
    }
    else if (zvec[0] != zvec[1]) {
        yvec[0] = zvec[1] ; yvec[1] = -zvec[0]; yvec[2] = 1;
    }
    else if ( zvec[0] != zvec[2]) {
        yvec[0] = zvec[2] ; yvec[1] = 1; yvec[2] = -zvec[0];
    }
    else {
        std::cerr << "\nERROR!  Cannot make vector perpendicular to bond axis!\n";
        std::cout << "\nBVEC: " << bondvector[0] << " " << bondvector[1] << " " << bondvector[2] << std::endl;
        std::cout << "ZVEC: " << zvec[0] << " " << zvec[1] << " " << zvec[2] << std::endl;
        std::exit(1);
    }
    if (dot(zvec,yvec) > 1e-5 ) {
        yvec[0] = -1*(zvec[1]*yvec[1]+zvec[2]*yvec[2])/zvec[0];
    }
    double ylen = pow(dot(yvec,yvec),0.5);
    for (int i=0; i<3;i++) {
        yvec[i] = yvec[i]/ylen;
    }
    double xvec[3];
    cross(yvec,zvec,xvec);
    double xlen = pow(dot(xvec,xvec),0.5);
    for (int i=0; i<3;i++) {
        xvec[i] = xvec[i]/xlen;
    }
    
    /* Find the points that are Cho group-sites for calculating the potential */
    std::vector<std::vector<double> > each_point;
    int nsites = 0;
    if (d->site_ndx.size() > 0) {
        /* 
         * There are 8 sites from -a1, 8 sites from -a2, 1 site along the bond vector
         * and 1 site for each atom in site_ndx
         */
        nsites = d->site_ndx.size() + 17;
        int n = 0;
        each_point = std::vector<std::vector<double> > (nsites, std::vector<double> (3,0));
        //std::vector<std::vector<double> > each_point(nsites,std::vector<double> (3,0));
        /* Assign coordinates for the atoms in site_ndx */
        for (int i=0; i<(int)d->site_ndx.size(); i++) {
            for (int j=0; j<3; j++) {
                each_point[n][j] = fr->x[d->site_ndx[i]][j];
            }
            n++;
        }
        /* Assign coordinate for the site along the bond vector */
        for (int i=0; i<3; i++) {
            each_point[n][i] = fr->x[d->a2][i] + bondvector[i] / bondlength * d->ring_dist;
        }
        n++;
        /* 8 points around -a1 */
        ring_points( fr->x[d->a1], &xvec[0], &yvec[0], &zvec[0], d->ring_dist, each_point,n);
        /* 8 points around -a2 */
        ring_points( fr->x[d->a2], &xvec[0], &yvec[0], &zvec[0], d->ring_dist, each_point,n);
        std::vector<double> minv(3,10);
        std::vector<double> maxv(3,-10);
    }
    
    /* Find the dummy atom points */
    /* midpoint xyz vectors */
    double pzcomp[3], nzcomp[3], pycomp[3], nycomp[3], pxcomp[3], nxcomp[3];
    /* -a1 xyz vectors */
    double pzcomp1[3], nzcomp1[3], pycomp1[3], nycomp1[3], pxcomp1[3], nxcomp1[3];
    /* -a2 xyz vectors */
    double pzcomp2[3], nzcomp2[3], pycomp2[3], nycomp2[3], pxcomp2[3], nxcomp2[3];

    for (int i=0; i<3; i++) {
        /* midpoint xyz vectors */
        pzcomp[i] = midpoint[i] + zvec[i] * d->delta;
        nzcomp[i] = midpoint[i] - zvec[i] * d->delta;
        pycomp[i] = midpoint[i] + yvec[i] * d->delta;
        nycomp[i] = midpoint[i] - yvec[i] * d->delta;
        pxcomp[i] = midpoint[i] + xvec[i] * d->delta;
        nxcomp[i] = midpoint[i] - xvec[i] * d->delta;
        /* -a1 xyz vectors */
        pzcomp1[i] = fr->x[d->a1][i] + zvec[i] * d->delta;
        nzcomp1[i] = fr->x[d->a1][i] - zvec[i] * d->delta;
        pycomp1[i] = fr->x[d->a1][i] + yvec[i] * d->delta;
        nycomp1[i] = fr->x[d->a1][i] - yvec[i] * d->delta;
        pxcomp1[i] = fr->x[d->a1][i] + xvec[i] * d->delta;
        nxcomp1[i] = fr->x[d->a1][i] - xvec[i] * d->delta;
        /* -a2 xyz vectors */
        pzcomp2[i] = fr->x[d->a2][i] + zvec[i] * d->delta;
        nzcomp2[i] = fr->x[d->a2][i] - zvec[i] * d->delta;
        pycomp2[i] = fr->x[d->a2][i] + yvec[i] * d->delta;
        nycomp2[i] = fr->x[d->a2][i] - yvec[i] * d->delta;
        pxcomp2[i] = fr->x[d->a2][i] + xvec[i] * d->delta;
        nxcomp2[i] = fr->x[d->a2][i] - xvec[i] * d->delta;
    }

    /* Write all the dummy atoms to the pqr file */
    write_dummy_atom( pqrout, nzcomp, "MIDz" );
    write_dummy_atom( pqrout, pzcomp, "MIDz" );
    write_dummy_atom( pqrout, nycomp, "MIDy" );
    write_dummy_atom( pqrout, pycomp, "MIDy" );
    write_dummy_atom( pqrout, nxcomp, "MIDx" );
    write_dummy_atom( pqrout, pxcomp, "MIDx" );
    
    write_dummy_atom( pqrout, nzcomp1, (std::string)*top->atoms.atomname[d->a1]+"z" );
    write_dummy_atom( pqrout, pzcomp1, (std::string)*top->atoms.atomname[d->a1]+"z" );
    write_dummy_atom( pqrout, nycomp1, (std::string)*top->atoms.atomname[d->a1]+"y" );
    write_dummy_atom( pqrout, pycomp1, (std::string)*top->atoms.atomname[d->a1]+"y" );
    write_dummy_atom( pqrout, nxcomp1, (std::string)*top->atoms.atomname[d->a1]+"x" );
    write_dummy_atom( pqrout, pxcomp1, (std::string)*top->atoms.atomname[d->a1]+"x" );
    
    write_dummy_atom( pqrout, nzcomp2, (std::string)*top->atoms.atomname[d->a2]+"z" );
    write_dummy_atom( pqrout, pzcomp2, (std::string)*top->atoms.atomname[d->a2]+"z" );
    write_dummy_atom( pqrout, nycomp2, (std::string)*top->atoms.atomname[d->a2]+"y" );
    write_dummy_atom( pqrout, pycomp2, (std::string)*top->atoms.atomname[d->a2]+"y" );
    write_dummy_atom( pqrout, nxcomp2, (std::string)*top->atoms.atomname[d->a2]+"x" );
    write_dummy_atom( pqrout, pxcomp2, (std::string)*top->atoms.atomname[d->a2]+"x" );
    
    if (d->site_ndx.size() > 0) {
        for (int i=0; i<nsites; i++) {
            std::stringstream key;
            key << "p" << i;
            write_dummy_atom( pqrout, &each_point[i][0], key.str());
        }
    }
    
    /* Initialize all the electrostatics, since they are the result of summing parts */
    int nexclude = d->exclude_ndx.size();
    dat.field_proj_total = 0.0f;
    dat.field_proj_each_exclude = std::vector<double> (nexclude,0.0f);
    dat.field_each_exclude = std::vector<std::vector<double> > (nexclude, std::vector<double> (3,0));
    dat.field_proj_exclude = 0.0f;
    for (int i=0; i<3; i++){
        dat.field_exclude[i] = 0.0f;
        dat.field_mid[i] = 0.0f;
        dat.field_a1[i] = 0.0f;
        dat.field_a2[i] = 0.0f;
    }
    dat.protein_site = std::vector<double> (nsites, 0.0f);
    dat.water_site = std::vector<double> (nsites, 0.0f);
    
    int nwater=0;
    /* Loop through everything */
    for (int g = 0; g < nr; ++g) { // for each group
        if (g > 1) {
            fprintf(stderr,"\nWarning: More than 1 group was specified.  Cowardly ignoring all additional groups\n");
            break;
        }
        for (int i = 0; i < sel[g]->p.nr; ++i) { // for each atom in the group
            int atom_ndx = sel[g]->g->index[i]; // how we get to the atom index
            /* Get .dat parameters and write the atom to the .pqr file */
            int resid = top->atoms.atom[atom_ndx].resind;
            double charge = top->atoms.atom[atom_ndx].q;
            char *atomname = *top->atoms.atomname[atom_ndx];
            char *resname = *top->atoms.resinfo[resid].name;
            double radius = -1;
            int namber = d->amber.size();
            for (int j=0; j<namber; j++) {
                if (d->bVerbose) {
                    fprintf(stderr, "\nTrying to match: %s %s to %s (%d/%d)", resname, *top->atoms.atomtype[atom_ndx], d->amber[j].ambername, j, namber);
                }
                if ( strncmp(*top->atoms.atomtype[atom_ndx], d->amber[j].ambername, strlen(*top->atoms.atomtype[atom_ndx])+1) == 0 ) {
                    radius = d->amber[j].radius;
                    if (d->bVerbose) {
                        fprintf(stderr,"\nMatched %s on %s to %s! (%d/%d)",resname, *top->atoms.atomtype[atom_ndx], d->amber[j].ambername, j, namber);
                    }
                    break;
                }
            }
            if (radius < 0) {
                fprintf(stderr,"\nError: Frame %d, could not match %s %s %s %i\n",d->framen, atomname, resname, *top->atoms.atomtype[atom_ndx],i+1);
                std::exit(1);
            }
            fprintf(pqrout,"ATOM %6i %4s %4s %5i %11.3f %7.3f %7.3f %7.4f %7.3f\n",atom_ndx+1,atomname,resname,resid+1,fr->x[atom_ndx][XX]*10,fr->x[atom_ndx][YY]*10,fr->x[atom_ndx][ZZ]*10,charge,radius);
            
            /* Do the field calculations now */
            bool keepAtom = true;
            calculate_field(midpoint, atom_ndx, *top, *fr, dat.field_mid, 1);
            /* Get the field due to each excluded atom */
            for (int j=0; j<nexclude ; j++) {
                if (atom_ndx == d->exclude_ndx[j]) {
                    calculate_field(midpoint, d->exclude_ndx[j], *top, *fr, &dat.field_each_exclude[j][0], 1);
                    keepAtom = false;
                    break;
                }
            }
            if (keepAtom) {
                calculate_field(midpoint, atom_ndx, *top, *fr, dat.field_exclude, 1);
                calculate_field(fr->x[d->a1], atom_ndx, *top, *fr, dat.field_a1, 1);
                calculate_field(fr->x[d->a2], atom_ndx, *top, *fr, dat.field_a2, 1);
            }
            /* Do the potential calculations now */
            /* We want to exclude atoms at the sites */
            /* We are also going to exclude atoms in the exclude group */
            if ((int)d->site_ndx.size() > 0) {
                if (keepAtom) {
                    for (int j=0; j<(int)d->site_ndx.size(); j++) {
                        if (atom_ndx == d->site_ndx[j] ) {
                            keepAtom = false;
                            break;
                        }
                    }
                }
                if (keepAtom) {
                    if (strncmp("SOL", resname, sizeof(resname)+1) == 0 ||
                        strncmp("HOH", resname, sizeof(resname)+1) == 0 ) {
                        for (int j=0; j<nsites; j++) {
                            double r = calculate_r( &each_point[j][0], atom_ndx, *top, *fr);
                            dat.water_site[j] += calculate_potential( r, atom_ndx, *top, *fr);
                        }
                        nwater++;
                    }
                    else {
                        for (int j=0; j<nsites; j++) {
                            double r = calculate_r( &each_point[j][0], atom_ndx, *top, *fr);
                            dat.protein_site[j] += calculate_potential( r, atom_ndx, *top, *fr);
                        }
                    }
                }
            }
        }
    }
    
    /* Write all output to xvg file */
    if (d->fp) {
        fprintf(d->fp, "%.6e ", d->rpdie * project_field(bondvector, dat.field_mid));
        for (int i=0; i<nexclude; i++) {
            fprintf(d->fp, "%.6e ", d->rpdie * project_field(bondvector, &dat.field_each_exclude[i][0]));
        }
        fprintf(d->fp, "%.6e ", d->rpdie * project_field(bondvector, dat.field_exclude));

        fprintf(d->fp, "%.6e ", d->rpdie * project_field(xvec,dat.field_exclude));
        fprintf(d->fp, "%.6e ", d->rpdie * project_field(yvec,dat.field_exclude));
        fprintf(d->fp, "%.6e ", d->rpdie * project_field(zvec,dat.field_exclude));

        fprintf(d->fp, "%.6e ", d->rpdie * project_field(xvec,dat.field_a1));
        fprintf(d->fp, "%.6e ", d->rpdie * project_field(yvec,dat.field_a1));
        fprintf(d->fp, "%.6e ", d->rpdie * project_field(zvec,dat.field_a1));
        
        fprintf(d->fp, "%.6e ", d->rpdie * project_field(xvec,dat.field_a2));
        fprintf(d->fp, "%.6e ", d->rpdie * project_field(yvec,dat.field_a2));
        fprintf(d->fp, "%.6e ", d->rpdie * project_field(zvec,dat.field_a2));
        
        if ((int)d->site_ndx.size() > 0) {
            for (int i=0; i<nsites; i++) {
                fprintf(d->fp, "%.6e %.6e ", d->rpdie * dat.protein_site[i], d->rpdie * dat.water_site[i]);
            }
        }
        fprintf(d->fp, "\n");
    }

    /* close the pqr file */
    fclose(pqrout);
    /* increment the frame number */
    d->framen++;
    /* We need to return 0 to tell that everything went OK */
    return 0;
}

/*! \brief
 * Function that implements the analysis tool.
 *
 * Following the style of Gromacs analysis tools, this function is called
 * \p gmx_something.
 */

int gmx_gmx2pqr(int argc, char *argv[])
{
    const char *desc[] = {
        "\tThis convert from gromacs to APBS pqr files.",
        "Use the -a1 flag for the start of the bond vector",
        "and the -a2 flag for the end of the bond vector.",
        "\n\tUse the -exclude flag to specify atoms which you",
        "would like to know the speficic contribution to the",
        "electrostatic field at the bond midpoint, or would",
        "like to exclude from the field calculation.",
        "Use the -site flag to add dummy atoms like those used",
        "by the Cho group {Choi, J. H.; Oh, K. I.; Lee, H.; Lee,",
        "C.; Cho, M.; Nitrile and thiocyanate IR",
        "probes: Quantum chemistry calculation studies and",
        "multivariate least-square fitting analysis; J. Phys.",
        "Chem. 2008}.  By default, both -site and -exclude flag",
        "will include the atoms indicated by -a1 and -a2.",
        "\n\tUse the -select flag to specify which atoms you",
        "want to include in the .pqr file and all calculations.",
        "For example, '(not resname SOL and not name Na) or",
        "resname SOL and same residue as within 0.5 of resname",
        "CNC and (name NE or name CD or name CB or name SG)",
        "would give protein and water near the heavy backbone",
        "atoms of any cyanocysteine residues. By default, when",
        "the -exclude flag is used, potentials will be reported",
        "in two group--water (resname HOH or resname SOL) and",
        "non-water contributions."
        "\n"
    };
    /* Command-line arguments */
    int             a1 = 0; // initialized for error checking
    int             a2 = 0;
    real            pdie = 1;
    real            delta = 0.010; // nm
    real            ring_dist = 0.07; // nm
    static char     *exclude = NULL;
    static char     *site = NULL;
    bool            bVerbose = false; // 0=false, 1=true
    
    t_pargs         pa[] = {
        { "-a1", TRUE, etINT,
            {&a1}, "Starting atom for bond vector--ie: CD in CNC"},
        { "-a2", TRUE, etINT,
            {&a2}, "Ending atom for bond vector--ie: NE in CNC"},
        { "-site", TRUE, etSTR,
            {&site}, "Atoms to be included as sites for the electrostatic model used by the Cho group.  Atoms pass by -a1 and -a2 are added to this group automatically.  Note that the whole selection string will need to be quoted so that your shell will pass it in as a string."},
        { "-exclude", TRUE, etSTR,
            {&exclude}, "Specific atoms to exclude from the field calculations.  Atoms pass by -a1 and -a2 are added to this group automatically.  Note that the whole selection string will need to be quoted so that your shell will pass it in as a string."},
        { "-pdie", TRUE, etREAL,
            {&pdie}, "Protein dielectric constant"},
        { "-delta", TRUE, etREAL,
            {&delta}, "(nm) Spacing between dummy atoms at the bond midpoint"},
        { "-v", FALSE, etBOOL,
            {&bVerbose}, "Be slightly more verbose"}
    };
    
    t_filenm        fnm[] = {
         { efDAT, "-d", "/Users/ritchie/Utilities/apbs/AMBER.DAT", ffREAD }
        ,{ efPQR, "-o", NULL, ffWRITE }
        ,{ efXVG, "-of", "coulomb_field", ffWRITE }
        ,{ efDAT, "-bw", "bw.dat", ffOPTWR}
        //,{ efNDX, NULL, NULL, ffREAD }
    };
    
#define NFILE asize(fnm)
    
    gmx_ana_traj_t      *trj;
    t_topology          *top;
    output_env_t        oenv;
    //t_dsdata            d;
    //t_analysisdata d;
    t_analysisdata      d;
    int                 ngrps;
    gmx_ana_selection_t **sel;
    int                 g;
    int                 rc;
    
    CopyRight(stderr, argv[0]);
    
    gmx_ana_traj_create(&trj, ANA_REQUIRE_TOP);
    gmx_ana_set_nanagrps(trj, -1);
    parse_trjana_args(trj, &argc, argv, 0,
                      NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL,
                      &oenv);
    gmx_ana_get_topology(trj, FALSE, &top, NULL);

    gmx_ana_get_nanagrps(trj, &ngrps);
    gmx_ana_get_anagrps(trj, &sel);
    gmx_ana_init_coverfrac(trj, CFRAC_SOLIDANGLE);
    
    /* Get output file names */
    const char          *datname, *framepqr;
    datname = opt2fn("-d", NFILE, fnm);
    framepqr = opt2fn("-o", NFILE, fnm);
    
    /* open xvg file */
    d.fp = xvgropen(opt2fn("-of", NFILE, fnm), "Coulomb field for each frame", "Frame Number", "Field[kbt/eA]", oenv);
    std::vector<std::string> legend;
    legend.push_back("frame");
    legend.push_back("total_field");
    /* I'm also going to make a helpful string for my bw_wham script */
    std::stringstream bw_wham_string;
    bw_wham_string << "-d coulomb_all_atoms -u kbt/eA ";
    
    /* Make sure -a1 and -a2 are included and increment by -1 to match internal numbering */
    if ( a1<1 || a2<1 ) {
        gmx_fatal(FARGS, "\nAtom numbers -a1 and -a2 defining the bond vector must be specified\n");
    }
    a1--; a2--;
    
    /* If any, get the atom numbers for the excluded atoms */
    std::vector<int> exclude_ndx;
    if (exclude) {
        int n;
        std::stringstream excluded_line(exclude);
        while ( excluded_line >> n ) {
            exclude_ndx.push_back(n-1);
        }
        exclude_ndx.push_back(a1);
        exclude_ndx.push_back(a2);
        fprintf(stderr,"\nWill calculate and remove field contributions due to atoms: ");
        for (int i=0; i<(int)exclude_ndx.size(); i++) {
            fprintf(stderr,"%i%s ",exclude_ndx[i]+1,*top->atoms.atomname[exclude_ndx[i]]);
            /* Add the excluded atoms to the xvg legend */
            std::stringstream key;
            key << exclude_ndx[i]+1 << *top->atoms.atomname[exclude_ndx[i]] << "_field";
            legend.push_back(key.str());
            bw_wham_string << "-d coulomb_" << *top->atoms.atomname[exclude_ndx[i]] << " -u kbt/eA ";
        }
    }
    /* Add entries to the xvg legend */
    legend.push_back("F-except");
    bw_wham_string << "-d coulomb_minus_ch2scn -u kbt/eA ";

    char xyz[] = "XYZ";
    /* midpoint field vector components */
    for (int i=0;i<3;i++){
        std::stringstream key;
        key << "MID" << xyz[i];
        legend.push_back(key.str());
        bw_wham_string << "-d coulomb_midpoint_" << xyz[i] << " -u kbt/eA ";
    }
    /* -a1 field vector components */
    for (int i=0;i<3;i++){
        std::stringstream key;
        key << *top->atoms.atomname[a1] << xyz[i];
        legend.push_back(key.str());
        bw_wham_string << "-d coulomb_" << *top->atoms.atomname[a1] << "_" << xyz[i] << " -u kbt/eA ";
    }
    /* -a2 field vector components */
    for (int i=0;i<3;i++){
        std::stringstream key;
        key << *top->atoms.atomname[a2] << xyz[i];
        legend.push_back(key.str());
        bw_wham_string << "-d coulomb_" << *top->atoms.atomname[a2] << "_" << xyz[i] << " -u kbt/eA ";
    }
    
    /* If any, get the atom numbers for the Cho site atoms */
    std::vector<int> site_ndx;
    if (site) {
        int n;
        std::stringstream site_line(site);
        while (site_line >> n ) {
            site_ndx.push_back(n-1);
        }
        site_ndx.push_back(a1);
        site_ndx.push_back(a2);
        fprintf(stderr,"\nWill calculate the electrostatic potential at sites on the following atoms: ");
        for (int i=0; i<(int)site_ndx.size(); i++) {
            fprintf(stderr,"%i%s ",site_ndx[i]+1,*top->atoms.atomname[site_ndx[i]]);
            /* Add the sites to the xvg legend */
            std::stringstream pkey, wkey;
            pkey << "protein_p" << i;
            legend.push_back(pkey.str());
            bw_wham_string << "-d coulomb_" << pkey.str() << " -u kbt/e ";
            wkey << "water_p" << i;
            legend.push_back(wkey.str());
            bw_wham_string << "-d coulomb_" << wkey.str() << " -u kbt/e ";
        }
        /* Add entries for each of the 8 sites on each -a1 and -a2 and the site straight out
           from the bond */
        for (int i=0; i<17; i++) {
            std::stringstream pkey, wkey;
            pkey << "protein_p" << i+(int)site_ndx.size();
            legend.push_back(pkey.str());
            bw_wham_string << "-d coulomb_" << pkey.str() << " -u kbt/e ";
            wkey << "water_p" << i+(int)site_ndx.size();
            legend.push_back(wkey.str());
            bw_wham_string << "-d coulomb_" << wkey.str() << " -u kbt/e ";
        }
    }
    
    /* convert the legend to a char array and include it in the xvg file */
    int nlegend_items = legend.size();
    char **flegend = new char* [nlegend_items];
    for (int i=0; i<nlegend_items; i++) {
        flegend[i] = new char [sizeof(legend[i].c_str())];
        strcpy(flegend[i],legend[i].c_str());
    }
    xvgr_legend(d.fp, nlegend_items, (const char**)flegend, oenv);
    
    /* Make sure that -a1, -a2, -exclude, and -site atoms exist in the trajectory */
    if (a1<0 || a1>top->atoms.nr){
        gmx_fatal(FARGS, "\nError: Atom -a1 is %d, which is outside of the range 0 to %d\n",a1+1,top->atoms.nr+1);
    }
    if (a2<0 || a2>top->atoms.nr){
        gmx_fatal(FARGS, "\nError: Atom -a2 is %d, which is outside of the range 0 to %d\n",a2+1,top->atoms.nr+1);
    }
    for (int i=0; i<(int)exclude_ndx.size(); i++){
        if (exclude_ndx[i]<0 || exclude_ndx[i]>top->atoms.nr){
            gmx_fatal(FARGS, "\nError: Excluded atom %d, %d, which is outside of the range 0 to %d\n",i,exclude_ndx[i]+1,top->atoms.nr+1);
        }
    }
    for (int i=0; i<(int)site_ndx.size(); i++){
        if (site_ndx[i]<0 || site_ndx[i]>top->atoms.nr){
            gmx_fatal(FARGS, "\nError: Cho site atom %d, %d, which is outside of the range 0 to %d\n",i,site_ndx[i]+1,top->atoms.nr+1);
        }
    }
    
    /* If asked for, print out the bw_wham_string to a file */
    if (opt2bSet("-bw", NFILE, fnm)) {
        fprintf(stderr,"\nWriting flags for bw_wham to %s",opt2fn("-bw", NFILE, fnm));
        std::ofstream bwfile;
        bwfile.open(opt2fn("-bw", NFILE, fnm));
        bwfile << bw_wham_string.str() << std::endl;
        bwfile.close();
    }
    
    /* Parse through the DAT file to make the gmx->dat mapping */
    read_DAT(datname, d.amber, bVerbose);

    /* Pass variables to structure for analysis */
    d.a1 = a1;
    d.a2 = a2;
    d.rpdie = 1./pdie;
    d.delta = delta * 0.5;
    d.ring_dist = ring_dist;
    d.site_ndx = site_ndx;
    d.exclude_ndx = exclude_ndx;
    d.framen = 0;
    d.pqr = framepqr;
    d.bVerbose = bVerbose;
    /* Parse through the frames */
    gmx_ana_do(trj, 0, &analyze_frame, &d);
    ffclose(d.fp);
}



/*! \brief
 * The main function.
 *
 * In Gromacs, most analysis programs are implemented such that the \p main
 * function is only a wrapper for a \p gmx_something function that does all
 * the work, and that convention is also followed here. */
int main(int argc, char *argv[])
{
    gmx_gmx2pqr(argc, argv);
    return 0;
}
