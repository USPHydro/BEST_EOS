
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <omp.h>

using std::cout;
using std::endl;
using std::ifstream;
using std::vector;
using std::isnan;
using std::isinf;
using std::istringstream;
using std::string;


// LET THE CODE FOR PRAIA FORMAT

#define hbarc 0.1973269631 // GeV*fm


int main(int argc, char** argv){

    string path = argv[1];
    //string string_tmp[] = {"1", "2", "3", "4", "5", "6", "7"};
    const int ntables = 7;

    //double  **p, **t, **u; //pressure, temperature and baryon potential




    for(int itable = 0; itable <= 6; ++itable){

        cout<<"Table "<<itable+1<<endl;

//        # e0= 3.60000000e-03 de= 1.93220339e-04 Ne= 60 rhob0= 0.00000000e+00 drhob0= 5.00000000e-05 Nrhob= 301 (only pressure archives)

        std::ifstream eos_p(path+"/BEST_eos_p_" + std::to_string(itable)+".dat");
        std::ifstream eos_t(path+"/BEST_eos_T_" + std::to_string(itable)+".dat");
        std::ifstream eos_u(path+"/BEST_eos_muB_" + std::to_string(itable)+".dat");
        std::ifstream eos_cs2(path+"/BEST_eos_cs2_" + std::to_string(itable)+".dat");
        
        string name = path+"/"+"visuEOS_"+std::to_string(itable+1)+".dat";
        FILE *output_file1 = fopen(name.c_str(), "w");

        name = path+"/"+"EOS_"+std::to_string(itable+1)+".dat";
        FILE *output_file2 = fopen(name.c_str(), "w");
    

        // Get the minimum values of rhob and e
        //getline(eos_p, line);
        //istringstream iss(line);
        //iss >> emin >> de >> ne >> nmin >> dn >> nn;

        double p, t, u, e, n, cs2;
        double nmin, emin, dn, de;
        int nn, ne;
        string line;

        double dummy;

        eos_p >> emin;
        eos_p >> de;
        eos_p >> ne;
        eos_p >> nmin;
        eos_p >> dn;
        eos_p >> nn;
        
        //nn++;
        //ne++;
std::cout<<"emin= "<<emin
                 <<" emax= "<<emin + double(ne)*de
                 <<" ne= "<<ne
                 <<" nmin= "<<nmin
                 <<" nmax= "<<nmin + double(nn)*dn
                 <<" nn= "<<nn<<std::endl;
        //std::cout<<"emin= "<<emin<<" de= "<<de<<" ne= "<<ne<<" nmin= "<<nmin<<" dn= "<<dn<<" nn= "<<nn<<std::endl;

        //cout<< emin << " \t " << de << " \t " << ne << " \t " << nmin << " \t " << dn << " \t " <<  nn << endl;
        fprintf(output_file2, "%d %g %g %d %g %g\n", ne, emin, emin + double(ne)*de, nn, nmin, nmin + double(nn)*dn);

        //fprintf(output_file1, "%g %g\n", nmin + dn*double(nn -1), emin + de*double(ne -1));
        for (int j = 0; j < ne; j++) {
            for (int i = 0; i < nn; i++) {
            
            
               
                n = dn*double(i);
                e = emin + de*double(j);

                eos_p >> p;
                eos_t >> t;
                eos_u >> u;
                eos_cs2 >> cs2;

                double s =(e+p-u*n)/t;

                if(t == 0.) fprintf(output_file2, "%g %g %g %g %g\n", t, u, 0., 0., 0.);
                else fprintf(output_file2, "%g %g %g %g %g\n", t, u, s, p, cs2);

                if(t == 0.) fprintf(output_file1, "%g %g %g %g %g %g %g\n", t, u, 0., 0., 0., 0., 0.);
                else fprintf(output_file1, "%g %g %g %g %g %g %g\n", t, u, e, n, s, p, cs2);
                
                //if(!isnan(s) && !isinf(s)) fprintf(output_file1, "%g %g\n", n, s);
            }
        }
        //cout<<"B: "<<l<<endl;

        eos_p.close();
        eos_t.close();
        eos_u.close();
        eos_cs2.close();
        
        fclose(output_file1);
        fclose(output_file2);

    }

    return 0;
}





/*void  initialize_eos() {

	bool flag_muS = false;
    bool flag_muC = false;
    string eos_file_string_array[7];
    if (eos_id == 12) {
        //music_message.info("reading EOS neos_b ...");
        spath << "/EOS/neos_b/";
        string string_tmp[] = {"1", "2", "3", "4", "5", "6", "7"};
        std::copy(std::begin(string_tmp), std::end(string_tmp),
                  std::begin(eos_file_string_array));
    } else if (eos_id == 13) {
        music_message.info("reading EOS neos_bs ...");
        spath << "/EOS/neos_bs/";
        string string_tmp[] = {"1s", "2s", "3s", "4s", "5s", "6s", "7s"};
        std::copy(std::begin(string_tmp), std::end(string_tmp),
                  std::begin(eos_file_string_array));
        flag_muS = true;
    } else if (eos_id == 14) {
        music_message.info("reading EOS neos_bqs ...");
        spath << "/EOS/neos_bqs/";
        string string_tmp[] = {"1qs", "2qs", "3qs", "4qs", "5qs", "6qs", "7qs"};
        std::copy(std::begin(string_tmp), std::end(string_tmp),
                  std::begin(eos_file_string_array));
        flag_muS = true;
        flag_muC = true;
    }
    
    const int ntables = 7;

    pressure_tb    = new double** [ntables];
    temperature_tb = new double** [ntables];
    mu_B_tb        = new double** [ntables];
    if (flag_muS) {
        mu_S_tb = new double** [ntables];
    }
    if (flag_muC) {
        mu_C_tb = new double** [ntables];
    }

    for (int itable = 0; itable < ntables; itable++) {
        std::ifstream eos_p("neos" + eos_file_string_array[itable]
                            + "_p.dat");
        std::ifstream eos_T("neos" + eos_file_string_array[itable]
                            + "_t.dat");
        std::ifstream eos_mub("neos" + eos_file_string_array[itable]
                              + "_mub.dat");
        std::ifstream eos_muS;
        std::ifstream eos_muC;
        if (flag_muS) {
            eos_muS.open("neos" + eos_file_string_array[itable]
                         + "_mus.dat");
        }
        if (flag_muC) {
            eos_muC.open("neos" + eos_file_string_array[itable]
                         + "_muq.dat");
        }
        // read the first two lines with general info:
        // first value of rhob, first value of epsilon
        // deltaRhob, deltaE, number of rhob points, number of epsilon points
        // the table size is
        // (number of rhob points + 1, number of epsilon points + 1)
        int N_e, N_rhob;
        eos_p >> nb_bounds[itable] >> e_bounds[itable];
        eos_p >> nb_spacing[itable] >> e_spacing[itable]
              >> N_rhob >> N_e;
        nb_length[itable] = N_rhob + 1;
        e_length[itable]  = N_e + 1;

        e_bounds[itable]  /= Util::hbarc;   // 1/fm^4
        e_spacing[itable] /= Util::hbarc;   // 1/fm^4

        // skip the header in T and mu_B files
        string dummy;
        std::getline(eos_T, dummy);
        std::getline(eos_T, dummy);
        std::getline(eos_mub, dummy);
        std::getline(eos_mub, dummy);
        if (flag_muS) {
            std::getline(eos_muS, dummy);
            std::getline(eos_muS, dummy);
        }
        if (flag_muC) {
            std::getline(eos_muC, dummy);
            std::getline(eos_muC, dummy);
        }

        // allocate memory for EOS arrays
        pressure_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                               e_length[itable]);
        temperature_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                                  e_length[itable]);
        mu_B_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                           e_length[itable]);
        if (flag_muS) {
            mu_S_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                               e_length[itable]);
        }
        if (flag_muC) {
            mu_C_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                               e_length[itable]);
        }

        // read pressure, temperature and chemical potential values
        for (int j = 0; j < e_length[itable]; j++) {
            for (int i = 0; i < nb_length[itable]; i++) {
                eos_p >> pressure_tb[itable][i][j];
                eos_T >> temperature_tb[itable][i][j];
                eos_mub >> mu_B_tb[itable][i][j];

                if (flag_muS) {
                    eos_muS >> mu_S_tb[itable][i][j];
                    mu_S_tb[itable][i][j] /= Util::hbarc;    // 1/fm
                }
                if (flag_muC) {
                    eos_muC >> mu_C_tb[itable][i][j];
                    mu_C_tb[itable][i][j] /= Util::hbarc;    // 1/fm
                }

                pressure_tb[itable][i][j]    /= Util::hbarc;    // 1/fm^4
                temperature_tb[itable][i][j] /= Util::hbarc;    // 1/fm
                mu_B_tb[itable][i][j]        /= Util::hbarc;    // 1/fm
            }
        }
    }
}*/



