// #define ARMA_ALLOW_FAKE_GCC
#include <mpi.h>
#include <tclap/CmdLine.h>
#include "xatu.hpp"

int main(int argc, char* argv[]){

    int procMPI_rank, procMPI_size, procMPI_name_length;
    char procMPI_name[10];

    MPI_Init (&argc,&argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &procMPI_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &procMPI_size);
    MPI_Get_processor_name (procMPI_name, &procMPI_name_length);
    if(procMPI_rank == 0){
        xatu::printHeaderGTFprovisional();
        std::cout << std::endl;

        std::cout << "+---------------------------------------------------------------------------+" << std::endl;
        std::cout << "|                              Parallelization                              |" << std::endl;
        std::cout << "+---------------------------------------------------------------------------+" << std::endl;
        std::cout << "Number of MPI processes: " << procMPI_size << std::endl;
        std::cout << "Number of OpenMP threads: " << omp_get_max_threads() << std::endl;

        std::cout << std::endl;
        std::cout << "+---------------------------------------------------------------------------+" << std::endl;
        std::cout << "|                              BSE k=k' TERMS                               |" << std::endl;
        std::cout << "+---------------------------------------------------------------------------+" << std::endl;
    }

    // INPUT ////////////////////////////////////////////////////////////////////////////////////
    int ncells = 149; 
    std::string outp_file  = "/mnt2/fscratch/users/uam20_res/resh000119/Cpp/xatu/xatu_MPI_v1/TEMPORARY/MoS2_base_PBE_QTrev2_MoP.outp"; 
    std::string exc_file = "/mnt2/fscratch/users/uam20_res/resh000119/Cpp/xatu/xatu_MPI_v1/examples/excitonconfig/MoS2_GAUSSIAN.txt"; 
    arma::icolvec nPointsLine = {150,150,150};  
    arma::mat kVertCoords = {{0,0,0},{0.5,0,0},{0.333333,0.333333,0},{0,0,0}};
    // END of the INPUT /////////////////////////////////////////////////////////////////////////

    xatu::ConfigurationCRYSTAL CRYSTALconfig(outp_file, ncells, true);
    xatu::ConfigurationExciton Excitonconfig(exc_file);
    xatu::ExcitonGTF GTFexc(Excitonconfig, CRYSTALconfig, "", "/mnt2/fscratch/users/uam20_res/resh000119/Cpp/xatu/xatu_MPI_v1/TEMPORARY/Integrals_Files/");

    int nLines = nPointsLine.size();
    arma::mat Klist      = arma::zeros<arma::mat>(arma::sum(nPointsLine)-(nLines-1), 3); 
    arma::mat Klist_frac = arma::zeros<arma::mat>(arma::sum(nPointsLine)-(nLines-1), GTFexc.ndim); 
    kVertCoords = kVertCoords.cols(0,GTFexc.ndim - 1);
    uint32_t countr = 0;
    for(int line = 0; line < nLines; line++){
        arma::rowvec k0 {kVertCoords.row(line)};
        arma::rowvec kend {kVertCoords.row(line + 1)};
        for(uint32_t nk = (0 + (line!=0)); nk < nPointsLine(line); nk++){
            arma::rowvec k_frac {k0 + (kend-k0)*(nk/(double)(nPointsLine(line)-1))};
            Klist.row(countr)      = k_frac*(GTFexc.Gbasis.t());
            Klist_frac.row(countr) = k_frac;
            countr++;
        }
    }
    Klist      = (Klist.cols(0,GTFexc.ndim -1)).t();
    Klist_frac = Klist_frac.t();
    uint32_t Nk_total = countr;  //total number of k-points in the whole trajectory
   
    uint32_t Nk_per_proc = Nk_total / procMPI_size;
    uint Nk_remainder    = Nk_total % procMPI_size;
    if(Nk_per_proc == 0){
        throw std::invalid_argument("Some MPI processes are idle");
    }

    uint p = 0;   // number of previous processes containing one more remainding k-point 
    while(p < static_cast<uint>(procMPI_rank) && p < Nk_remainder){
        p++;
    }
    uint B = static_cast<uint>(procMPI_rank) < Nk_remainder;
    arma::mat kpoints_thisproc = Klist_frac.cols( Nk_per_proc*procMPI_rank + p, Nk_per_proc*(procMPI_rank+1) - 1 + B + p );
    uint32_t nk_thisproc = kpoints_thisproc.n_cols; 

    arma::mat energ_thisproc(GTFexc.norbitals, nk_thisproc);
    MPI_Barrier (MPI_COMM_WORLD);
    for(uint32_t kind = 0; kind < nk_thisproc; kind++){ // iterate over k-points 
        arma::colvec k = kpoints_thisproc.col(kind);

        arma::colvec eigval1;
        arma::cx_mat eigvec1;
        GTFexc.solveBands(k,  eigval1, eigvec1);
        energ_thisproc.col(kind) = eigval1;
    }
    MPI_Barrier (MPI_COMM_WORLD);

    for(int r = 0; r < procMPI_size; r++){
        if(procMPI_rank == r){
           if(procMPI_rank == 0){
               std::ofstream output_file0("energMPI", std::ios::trunc);
               output_file0.close();
               std::ofstream output_file("energMPI", std::ios::app);
               output_file << GTFexc.norbitals << ", " << Nk_total << ", " << procMPI_size << std::endl;
            }
            std::ofstream output_file("energMPI", std::ios::app);
            output_file << nk_thisproc << std::endl;
            output_file.precision(12);
            output_file << std::scientific;
            for(uint32_t i = 0; i < GTFexc.norbitals; i++){
                for(uint32_t j = 0; j < nk_thisproc; j++){
                    output_file << energ_thisproc(i,j) << " ";
                }
                output_file << std::endl;
            }
        }
        MPI_Barrier (MPI_COMM_WORLD);
    }    

    MPI_Finalize();
    return 0;

}