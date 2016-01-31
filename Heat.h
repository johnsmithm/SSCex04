#include"Vector.h" //Expr_CG function
#include <fstream>
#include <string>

using namespace std;

string dir = "../vtk/";

string upStep = "<?xml version=\"1.0\"?>\n<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\">\n <PUnstructuredGrid>\n\n<PPoints>\n<PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"/>\n</PPoints>\n\n <PPointData>\n <PDataArray type=\"Float64\" Name=\"Temperature\" NumberOfComponents=\"1\" format=\"ascii\"/>\n </PPointData>\n";
string downStep = " </PUnstructuredGrid>\n</VTKFile>";



void write_files(int rank, int timestep, Vector* u){
	ofstream out(dir+"time_step_"+to_string(timestep)+"_"+to_string(rank)+".vtu");
    out&(*u);
    out.close();
}

void run(int nx, int ny, int c, double eps, int timesteps, double tau, double k, double alfa, int vtk_spacing, int nrpr){
	  // Create Cart
	  int py = getHeight(nrpr);
      int px = nrpr/py;      
	  int dims[2] = {px,py};
	  int periods[2] = {0,0}; // not periodic!
      const int reorder = 1;  // allow reordering of process ranks      
      int cartrank(0);
    // Coords[0] - y, coords[1] - x.
      int coords[2] = {0,0};
    // up,down,left, right.
      int nbrs[4] = {0,0,0,0};
      MPI_Cart_create( MPI_COMM_WORLD, 2, dims, periods, reorder, &cartcomm );
      MPI_Comm_rank( cartcomm, &cartrank );
      MPI_Cart_coords( cartcomm, cartrank, 2, coords );
      MPI_Cart_shift( cartcomm, 0, 1, &nbrs[LEFT], &nbrs[RIGHT] );
      MPI_Cart_shift( cartcomm, 1, 1, &nbrs[UP], &nbrs[DOWN] );
	  std::swap(coords[0],coords[1]);

	  // Write files
	  ofstream out(dir+"solution.pvd");
	  if(cartrank == 0){
        std::cout<<"Height="<<py<<", Width="<<px<<";\n";
        
        out<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\">\n <Collection>\n";
        out<<"<DataSet timestep=\"0\" file=\"time_step_0.pvtu\"/>\n";

        ofstream outT(dir+"time_step_0.pvtu");
		outT<<upStep;
		for(int j=0;j<nrpr; ++j){
			 outT<<" <Piece Source=\"time_step_0_"+to_string(j)+".vtu\"/>\n";
			}
		outT<<downStep;
		outT.close();
      }

    // Initiate vectors  
    double pi = 3.141592653589793;
    // Old solution
    Vector * u = get_vector(cartrank, nx, ny, px, py);
    double offsetx = u->offsetx__();
    double offsety = u->offsety__();
    //sin(πx) sin(πy)
    u->set_distace(1.0/nx,1.0/ny);
    u->set_value([pi,nx,ny,offsetx,offsety](size_t x,size_t y)->double{
    	return sin(pi*(x+offsetx)*(1.0/(nx+1)))*sin((1.0/(ny+1))*pi*(y+offsety));});

    //write state 0
    write_files(cartrank, 0, u);

    // explicit method condition
	if(tau > ((1.0/(nx*ny))/(2*k)) && cartrank == 0 && alfa == 0){
    		cout<<"maybe unstable!!\n";
    }

    //New solution
    Vector * uNew = get_vector(cartrank, nx, ny, px, py);
    uNew->set_distace(1.0/nx,1.0/ny);

    //Stencil(int nx1,int ny1, int * nbrs_,double tau_,double k_,double alfa_, int sign=1)
    Stencil A(nx,ny, nbrs, tau,k, alfa, -1);   	

    // For each timeSpep    
    	for(int i=1; i<=timesteps; ++i){
    		//return;
    		 if(cartrank == 0)
    		 	cout<<"Step:"<<i<<"\n";

    		if(alfa == 1){//implicit
				Expr_CG(nx, ny, c, eps, py, px, nbrs, coords, cartrank, tau, k, alfa, *u);//kind of ok
			}else if(alfa == 0){//explicit
				*uNew = A*(*u);
			 	*u = *uNew;//should decrease not increase??
			 	if(cartrank == 0)
			 		cout<<"number iteration:0;\n";
			}else{//Grank-Nicolson
				*uNew = A*(*u);
			 	*u = *uNew;//stay the same??
			 	Expr_CG(nx, ny, c, eps, py, px, nbrs, coords, cartrank, tau, k, alfa, *u);
			}

			 if((i != 0 && ((i) % (vtk_spacing+1)) == 0)){
               //write files			 	
			 	if(cartrank == 0){
			 		out<<"<DataSet timestep=\""+to_string(i)+"\" file=\"time_step_"+to_string(i)+".pvtu\"/>\n";
			 		ofstream outT(dir+"time_step_"+to_string(i)+".pvtu");
			 		outT<<upStep;
			 		for(int j=0;j<nrpr; ++j){
			 			outT<<" <Piece Source=\"time_step_"+to_string(i)+"_"+to_string(j)+".vtu\"/>\n";
			 		}
			 		outT<<downStep;
			 		outT.close();			 		
			 	}
			 	u->set_distace(1.0/nx,1.0/ny);
			 	write_files(cartrank, i, u);
			 }
		}
    
    if(cartrank == 0){
    	out<<"</Collection>\n</VTKFile>\n";
    	out.close();
	}
	
}



