#include<fstream>
#include <cassert>
#include<iostream>
#include <stdlib.h>

#include <mpi.h>

#include"Heat.h"  

#include "Timer.h"


using namespace std;


int main(int argc, char *argv[]){ 
    
  (void) argc; //to suppress Warnings about unused argc
  int size(0); 
  int rank(0);
	
	MPI_Init( &argc, &argv );
	
	MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	
  assert(argc>8);
  //mpirun -np N ./heat nx ny c eps timesteps τ κ α vtk_spacing,
  int nx   = atoi(argv[1]);
  int ny   = atoi(argv[2]);
  int c    = atoi(argv[3]);	
  double eps  = stod(argv[4]);
  int timesteps = atoi(argv[5]);
  double tau  = stod(argv[6]);
  double k  = stod(argv[7]);
  double alfa  = stod(argv[8]);  
  int vtk_spacing = atoi(argv[9]);


    
   
  siwir::Timer timer;    
	
	//using template expretions
  run(nx, ny, c, eps, timesteps, tau, k, alfa, vtk_spacing, size);

    
  double time = timer.elapsed();


    
   if(rank == 0)
    cout<<"Total time:"<<time<<"\n";
    
 MPI_Finalize();
}
