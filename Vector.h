
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <math.h>
#include <functional>

const int UP    = 0;
const int DOWN  = 1;
const int LEFT  = 2;
const int RIGHT = 3;
MPI_Comm cartcomm( MPI_COMM_NULL );


template<class A>
	struct  Expr {
    operator const A& () const{
	return *static_cast<const A*>(this);
	}
};

//template <class A, class B>	class Add : public Expr<Add<A,B>>;


class Vector : public Expr<Vector > {

public :
    /**
    * set the derived datatype.
    */
	void commit_datatypes(){
       //             #blocks  bl-len   el.bet.bl   old-dt      new-dt 
	   MPI_Type_vector(   nyi,       1,       nxs,  MPI_DOUBLE, &columntype );
	   MPI_Type_commit( &columntype );
	
	   MPI_Type_vector(    ny,      nx,     nxs,   MPI_DOUBLE,  &printtype );
	   MPI_Type_commit( &printtype );	 
	}
	
    /**
    * Initiate with a value all grid points.
    *
    * n - total number of grid points(+ goast layers).
    * nx - Number inner grid points on x axe.
    * ny - Number inner grid points on y axe.
    * nxi - Number inner grid points that are updated on x axe.
    * nyi - Number inner grid points that are updated on y axe.
    * nxs - Total Number  grid points on x axe.
    * nys - Total Number  grid points on y axe.
    * x - Coordonate of fist inner point(even at order)
    * y - Coordonate of fist inner point(even at order)
    */
	Vector(int n_,int nx_,int ny_,int nxi_,int nyi_,int nxs_,int nys_,int x_,int y_, double w = 0)//d,r,z
		:n(n_),x(x_),y(y_),nx(nx_),ny(ny_),nxs(nxs_),nys(nys_),nxi(nxi_),nyi(nyi_){  		
       data = new double[n];
	   for(size_t i=0;i<n;++i)data[i]=w;
	   commit_datatypes();		
	}
	
    /**
    * Initiate with a function all grid points based on the position on the original grid.
    *
    * n - total number of grid points(+ goast layers).
    * nx - Number inner grid points on x axe.
    * ny - Number inner grid points on y axe.
    * nxi - Number inner grid points that are updated on x axe.
    * nyi - Number inner grid points that are updated on y axe.
    * nxs - Total Number  grid points on x axe.
    * nys - Total Number  grid points on y axe.
    * x - Coordonate of fist inner point(even at order)
    * y - Coordonate of fist inner point(even at order)
    * lambda - 4π^2 sin(2πx) sinh(2πy).
    */
	 Vector(int n_,int nx_,int ny_,int nxi_,int nyi_,int nxs_,int nys_,int x_,int y_,
			std::function<double(size_t,size_t)> f) //f
		 :n(n_),x(x_),y(y_),nx(nx_),ny(ny_),nxs(nxs_),nys(nys_),nxi(nxi_),nyi(nyi_)
	{	
	  data = new double[n];	 
	  for(size_t i=0;i<n;++i)data[i]=0;		 
	  for(size_t i=0;i<ny;++i)
		      for(size_t j=0;j<nx;++j)
			    data[j+x+nxs*(i+y)] = f(j,i);
	  commit_datatypes();	  
	}
	
    /**
    * Initiate with a function last row of grid points based on the position on the original grid.
    *
    * n - total number of grid points(+ goast layers).
    * nx - Number inner grid points on x axe.
    * ny - Number inner grid points on y axe.
    * nxi - Number inner grid points that are updated on x axe.
    * nyi - Number inner grid points that are updated on y axe.
    * nxs - Total Number  grid points on x axe.
    * nys - Total Number  grid points on y axe.
    * x - Coordonate of fist inner point(even at order)
    * y - Coordonate of fist inner point(even at order)
    * lambda - sin(2πx) sinh(2πy).
    */
	    Vector(int n_,int nx_,int ny_,int nxi_,int nyi_,int nxs_,int nys_,int x_,int y_, std::function<double(size_t)> f)//u
			:n(n_),x(x_),y(y_),nx(nx_),ny(ny_),nxs(nxs_),nys(nys_),nxi(nxi_),nyi(nyi_)
	{   	
	      data = new double[n];	
	      for(size_t i = 0;i < n; ++i){
              data[i] = 0;
          }
				
	      for(size_t i = 0;i < nx; ++i){
			  data[nxs * (nys - 1) + x + i] = f(i);
          }
		  commit_datatypes();
	}
    //end Vector Constructors
    
    
    //Vector(Vector&& o) noexcept : data(std::move(o.data)) {std::cout<<'|';}
	~Vector(){
		delete [] data;
		MPI_Type_free( &columntype );
		MPI_Type_free( &printtype );
	}
	
	//vector operators
	double operator[] (int i)const{
	return data[i];
	}
	
	double& operator[] (int i){
	return data[i];
	}
	
    /**
    * Assign operator(used with expresions).
    * Assign just inner point that are updated.
    */
	template<class A>
	void operator = (const Expr<A>& a_){
	const A& a(a_);
		
		for(size_t i = 0;i < nyi; ++i) {
			for(size_t j = 0;j < nxi; ++j) {
                data[j+1+(1+i)*nxs] = a[j+1+(1+i)*nxs];
            }
        }
	}

	double& operator()(size_t j, size_t i){
		return data[(i+y)*nxs+x+j];
	}

	double operator()(size_t j, size_t i)const{
		return data[(i+y)*nxs+x+j];
	}
	
    /*
    * Assign operator(with vector).
    */
	void operator = (const Vector& a)
	{
		n = a.n;
        x = a.x;
        y = a.y;
        nx = a.nx;
        ny = a.ny;
        nxs = a.nxs;
        nys = a.nys;
        for(size_t i = 0; i < n; ++i){
			data[i] = a[i];
        }
	}
	
    /*
    * Move operator.
    */
	void operator = ( Vector&& a)
	{
		   // std::cout<<"-";
		    n=a.n;x=a.x;y=a.y;nx=a.nx;ny=a.ny;nxs=a.nxs;nys=a.nys;
		    nxi=a.nxi;nyi=a.nyi;
		    offsetx=a.offsetx;offsety=a.offsety;
			data = a.data;
			a.data = NULL;
	}
	
    /**
    * Scalar product of two vectors.
    */
	double operator^(const Vector & a)const{
		double l = 0;		
		for(size_t i=0;i<nyi;++i)
			for(size_t j=0;j<nxi;++j)			
			      l+=data[j+1+(1+i)*nxs]*a.data[j+1+(1+i)*nxs];
		double l1=0;		
        //           sent rec.  #el   type      operation
	    MPI_Allreduce( &l, &l1 ,1 , MPI_DOUBLE ,MPI_SUM, cartcomm );
		return l1;
	}
	//end vector operation
    
	//begin vector members
    
    /**
    * Calculate Residual per all processes(not sqrt).
    */
	double LNorm(){
		double l = 0;		
		for(size_t i=0;i<nyi;++i)
			for(size_t j=0;j<nxi;++j)
			    l+=data[j+1+(1+i)*nxs]*data[j+1+nxs*(i+1)];
		double l1;	
        // Reduce to all residual from all processes.
        //           sent rec.  #el   type      operation
	    MPI_Allreduce( &l, &l1 ,1 , MPI_DOUBLE ,MPI_SUM, cartcomm );	
		return l1;
	}

	/**
	 *
	 */
	void set_value( std::function<double(size_t,size_t)> f){
		for(size_t i = 0;i < nyi; ++i) {
			for(size_t j = 0;j < nxi; ++j) {
                data[j+x+(y+i)*nxs] = f( j+offsetx, i+offsety);
            }
        }
	}
	
    /**
    * Get partial solution from the process with rank - rank.
    */
	void get_info(int rank){
	               MPI_Status status;
				   MPI_Recv( &data[nxs+x], 1, printtype, rank, 10, cartcomm, &status );
	}
    
    /**
    * Sent partial solution to 0 process.
    */
	void sent_info(){        
	               MPI_Send( &data[nxs+x],    1, printtype, 0, 10, cartcomm );
	}
	
    /**
    * Set distance between grid points in printing sections.
    */
	void  set_distace(double const x1,double const y1){
        hx=x1;
        hy=y1;
    }
    
    /**
    * Set offset  in printing sections.
    */
	void  set_offset(int const x1,int const y1){
        offsetx=x1;
        offsety=y1; 
    }
	
    size_t size()const{return n;}
	size_t offsetx__()const{return offsetx;}
	size_t offsety__()const{return offsety;}
	double hx__() const{return hx;}
	double hy__() const{return hy;}
	size_t nx__() const{return nx;}
	size_t ny__() const{return ny;}
	size_t nxs__()const{return nxs;}
	size_t nys__()const{return nys;}
	size_t x__()  const{return x;}
	size_t y__()  const{return y;}
	//end vector members
    
    //template <class A, class B> friend class Add ;//: public Expr<Add<A,B>>
    
    public:
    // initialization order
    //:n(n_),x(x_),y(y_),nx(nx_),ny(ny_),nxs(nxs_),nys(nys_),nxi(nxi_),nyi(nyi_)
	
    //number of grid points
	size_t n;
    //top-left grid point that are inner point not from goast layer.
    size_t x;
	size_t y;
    //inner point + boundaries 
	size_t nx;
    size_t ny;
    //inner point + boundaries + layers
	size_t nxs;
    size_t nys;
    // inner points (are updated)
	size_t nxi;
    size_t nyi;
    
    double *data;
    //distance between grid points
	double hx;
    double hy;
	size_t offsetx;
    size_t offsety;
	MPI_Datatype columntype; 
	MPI_Datatype printtype;
};

//vector print operators

/**
* Print for testing purpose.
*/
std::ostream & operator<<( std::ostream & os, const Vector & v )
{

	 for( size_t i=0; i < v.nys__(); ++i )
		{
			for( size_t j=0; j < v.nxs__(); ++j )
			{
				os << v[i*v.nxs__()+j] << " ";        
			}  
		    os<<"\n";
		}
	
	os<<"\n";
    return os;
}

/**
* Print for solution.
* print inner points + boundaries.
*/
std::ostream & operator>>( std::ostream & os, const Vector & v )
{
	
	 for( size_t i=0; i < v.ny__(); ++i )
		{
			for( size_t j=0; j < v.nx__(); ++j )
			{  
				os<<(j+v.offsetx__())*v.hx__()<<' '<<(i+v.offsety__())*v.hy__()<<' ' << v[(v.y__()+i)*v.nxs__()+j+v.x__()] << "\n";        
			}  
		   
		}
	
	os<<"\n";
    return os;
}

/**
* Print for solution.
* print inner points + boundaries.
*/
std::ostream & operator&( std::ostream & os, const Vector & v )
{
	os<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n <UnstructuredGrid>\n";
	os<<"<Piece NumberOfPoints=\""+std::to_string(v.ny__()*v.nx__())+"\" NumberOfCells=\""+std::to_string(v.ny__()*v.nx__())+"\">\n";
	os<<"<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	std::cout<<v.hx__()<<"--------------\n";
	 for( size_t i=0; i < v.ny__(); ++i )
		{
			for( size_t j=0; j < v.nx__(); ++j )
			{  
				os<<(j+v.offsetx__())*v.hx__()<<' '<<(i+v.offsety__())*v.hy__()<<' ' <<"0.00\n"; //<< v[(v.y__()+i)*v.nxs__()+j+v.x__()] << "\n";        
			}
		   
		}

	os<<"</DataArray>\n </Points>\n<Cells>\n<DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\">\n";	
	int ii = 0;
	for( size_t i=0; i < v.ny__(); ++i )
			for( size_t j=0; j < v.nx__(); ++j )os<<ii++<<" ";
	os<<"\n</DataArray>\n<DataArray type=\"UInt32\" Name=\"offsets\" format=\"ascii\">\n";	
	ii = 0;
	for( size_t i=0; i < v.ny__(); ++i )
			for( size_t j=0; j < v.nx__(); ++j )os<<++ii<<" ";		
	os<<"\n</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	for( size_t i=0; i < v.ny__(); ++i )
			for( size_t j=0; j < v.nx__(); ++j )os<<"1 ";
	os<<"\n</DataArray>\n</Cells>\n<PointData>\n<DataArray type=\"Float64\" Name=\"Temperature\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	for( size_t i=0; i < v.ny__(); ++i )
		{
			for( size_t j=0; j < v.nx__(); ++j )
			{  
				os<< v[(v.y__()+i)*v.nxs__()+j+v.x__()] << " ";        
			}  
		   os<<"\n";
		}	
	os<<" </DataArray>\n</PointData>\n  </Piece>\n </UnstructuredGrid>\n</VTKFile>";		
    return os;
}

//vector end



class Stencil : public Expr<Stencil > {


public :
    /**
    * Construct the Matrix based on Stencil.
    */
	Stencil(int nx1,int ny1, int * nbrs_,double tau_,double k_,double alfa_)
	:nx_(nx1+1),ny_(ny1+1),pi(3.141592653589793),tau(tau_), k(k_), alfa(alfa_),nbrs(nbrs_){
	        double hx_ = 1.0/nx1;
            double hy_ = 1.0/ny1;
            xst =  (tau*k)/(hx_*hx_);
            yst =  (tau*k)/(hy_*hy_);
            mst = 2.0/(hx_*hx_)+2.0/(hy_*hy_);          
            pg = (nx1+1)*(ny1+1);  
	}
	~Stencil(){}	
	/*
	Vector operator - (Vector& u){
		//std::cout<<"--";
		 Vector r(pg,u.nx__(),u.ny__());
		 for(int i=1;i<ny_-1;++i)
					for(int j=1;j<nx_-1;++j)
                    	r[i*nx_+j] = mst*u[i*nx_+j] - xst*(u[i*nx_+j+1]+u[i*nx_+j-1])-yst*(u[j+(i+1)*nx_]+u[j+(i-1)*nx_]);
		 return r;
	}*/
	
	
	
	int nx() const{return nx_;}
	int ny() const{return ny_;}
	
	int left() const{return nbrs[LEFT];}
	int right()const{return nbrs[RIGHT];}
	int up()   const{return nbrs[UP];}
	int down() const{return nbrs[DOWN];}
	
	double yst_()const{return yst;}
	double xst_()const{return xst;}
	double mst_()const{return mst;}

	double tau__()const{return tau;}
	double k__()const{return k;}
	double alfa__()const{return alfa;}
	private :
        // grid size.
		int nx_,ny_,pg;
		
        // coefficients for grid point update.
		double yst, xst, mst;
		const double pi;
		double tau, k, alfa;
		int * nbrs;
};//stencil end

//operation begin

template <class A, class B>
	class Add : public Expr<Add<A,B>> {
	const A& a_;
	const B& b_;
	public :
        /**
        * Vector Adition.
        */
		Add(const A& a,const B& b): a_(a),b_(b){}
		double operator[](int i)const{
		return (a_[i] + b_[i]);
		}
    
		int x__()const{return b_.x__();}
		int y__()const{return a_.y__();}
		int nxs__()const{return b_.nxs__();}
	};

template <class A, class B>
	class Minus : public Expr<Minus<A,B>> {
	const A& a_;
	const B& b_;
	public :
        /**
        * Vector Substraction.
        */
		Minus(const A& a,const B& b): a_(a),b_(b){}
		
		int x__()const{return b_.x__();}
		int y__()const{return a_.y__();}
		int nxs__()const{return b_.nxs__();}
		
		double operator[](int i)const{
		return (a_[i] - b_[i]);
		}
		
	};

template <class A>
	class Add<A,double> : public Expr<Add<A,double>> {
	const A& a_;
	const double& b_;
	public :
        /**
        * Vector scalart multiplication.
        */
		Add(const A& a,const double& b): a_(a),b_(b){}
		double operator[](int i)const{
		return (a_[i] * b_);
		}		
		int x__()const{return a_.x__();}
		int y__()const{return a_.y__();}
		int nxs__()const{return a_.nxs__();}
	};

template <>
	class Add<Stencil,Vector> : public Expr<Add<Stencil,Vector>> {
	const Vector& a_;
	const Stencil& b_;
    MPI_Request reqs[8];
    MPI_Status stats[8];	
	public :
        /**
        * Matrix(Stencil) - Vector multiplication.
        */
		Add(const Vector& a,const Stencil& b): a_(a),b_(b){
            //columntype - nxi blocks, of 1 elem, with nxs distance between blocks.
            
            //         begin adress,                   #elem   #datatype     #rank     tag envirom.  request
			MPI_Isend( &a.data[1+a.nxs],                a.nxi, MPI_DOUBLE,   b.up()  ,  0, cartcomm, &reqs[0]   ); // first row
			MPI_Isend( &a.data[1+a.nxs*(a.nys-2)],      a.nxi, MPI_DOUBLE,   b.down(),  1, cartcomm, &reqs[1]   ); //last row
			MPI_Isend( &a.data[1+a.nxs],                1,     a.columntype, b.left() , 2, cartcomm, &reqs[2]   ); //first column
			MPI_Isend( &a.data[(a.nxs*2)-2],            1,     a.columntype, b.right(), 3, cartcomm, &reqs[3]   );//last column

			MPI_Irecv( &a.data[1],                      a.nxi, MPI_DOUBLE,   b.up(),    1, cartcomm, &reqs[4]   );//first row
			MPI_Irecv( &a.data[1+a.nxs*(a.nys-1)],      a.nxi, MPI_DOUBLE,   b.down(),  0, cartcomm, &reqs[5]   );//last row
			MPI_Irecv( &a.data[-1+(a.nxs*2)],            1,     a.columntype, b.right(), 2, cartcomm, &reqs[6]   );//last column
			MPI_Irecv( &a.data[a.nxs],                  1,     a.columntype, b.left(),  3, cartcomm, &reqs[7]   );//first column

            // Wait till all the messages are sent and receive.
            MPI_Waitall( 8, reqs, stats );
		}

		bool diagonal(int i)const{
			if(a_.offsetx__() + (i % a_.nxs__()) - 1 == a_.offsety__() + (i / a_.nxs__()))
				return true;
			return false;
		}
        
        /**
        * Update formula.
        */
		double operator[](int i)const{
		return (b_.mst_()*a_[i] - b_.xst_()*(a_[i+1]+a_[i-1])-b_.yst_()*(a_[i+a_.nxs__()]+a_[i-a_.nxs__()]));			
		}
        
		int x__()const{return a_.x__();}
		int y__()const{return a_.y__();}
		int nxs__()const{return a_.nxs__();}
	};



template <class A, class B>
	inline Add<A,B> operator+ (const Expr<A>& a, const Expr<B>& b){
	return Add<A,B>(a,b);
	}

template <class A, class B>
	inline Minus<A,B> operator- (const Expr<A>& a, const Expr<B>& b){
	return Minus<A,B>(a,b);
	}

template <class A>
	inline Add<A,double> operator* (const Expr<A>& a, const double& b){
	return Add<A,double>(a,b);
	}

template <class A>
	inline Add<Stencil,A> operator* ( const Stencil& b, const Expr<A>& a){
	return Add<Stencil,A>(a,b);
	}
//operation end

/**
* Get the matrix of length as the one from process with rank = rank
*
* nx - Number grid points on x axe.
* ny - Number grid points on y axe.
* px - Cart's Width
* py - Cart's Height
*
* return Vector with dimensions as on process rank.
*/
Vector *get_vector(int rank,int nx,int ny,int px = 2,int py = 2){
	    //get the coorddonate of the process.
        // coords[0] - y, coords[1] - x.
	    int coords[2] = {0,0};
        MPI_Cart_coords(cartcomm, rank, 2, coords);
	    std::swap(coords[0],coords[1]);
	//number of inner points in current process on x and y axes.
	    int nx_ = (nx+1)/px , ny_ = (ny+1)/py; 
    //number of inner points for stencil calculation(point that we will update)
		int nxi = nx_, nyi = ny_;
    //total size of a block (+ goast layers)
		int nxs = nx_ , nys = ny_;
     //first inner point coordonates
		int x=0,y=0;
    
		int offsetx=coords[1]*nx_, offsety =coords[0]*ny_;	

		if(coords[0]==0 || coords[0]==py-1){//if we are last or first row
			if(coords[0]==py-1)
				ny_ += (ny+1)%py; // add the remaining rows
            if(py != 1)
			    nys = ny_+1; // we have just one layer because we are at bondary
		}else nys +=2; //add the layers

		if(coords[1]==0 || coords[1]==px-1){//if we are last or first column
			if(coords[1]==px-1)
				nx_ += (nx+1)%px; // add the remaining columns
            if(px != 1)
			    nxs = nx_+1; // we have just one layer because we are at bondary
		}else nxs +=2; //add the layers


		if(coords[0]==0 && coords[1]==0){//we are the process on the top-left
			x=0;y=0; //we do not have layer on top and on left
		}else if(coords[0]==0){//first column
		   x=1;y=0; // we dont have layers on top
		}else if(coords[1]==0)//first row
		{
		x=0; y=1; // we do not have layer on left
		}else{
			x=1;y=1;//we have layer on top and on left
		}//others

		nxi = nx_;nyi = ny_;
		if(coords[0]==0 || coords[0]==py-1){//first row or last row
			--nyi; // we are at the boundary	
            if(py == 1)--nyi;
		}
		if(coords[1]==0 || coords[1]==px-1){//first column
			--nxi; // we are at the boundary
            if(px == 1) --nxi;
		}
    //number of grid points
        int pg = nxs*nys;
    //create an matrix full with 0
		Vector * r = new Vector(pg,nx_,ny_,nxi,nyi,nxs,nys,x,y);	
    // get the info from process rank
		r->set_offset(offsetx,offsety);
		return r;
}

int primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89};

/**
* Compute the height of the cart.
*
* size - number of processes
*
* return - Cart's height.
*/
int getHeight(int size){
    int height = 1;
    bool use = true;
    for(int i = 20;i >= 0;--i){
        while(size%primes[i] == 0 && size != 1){
            size /= primes[i];
            use = !use;
            if(use) {
                height *= primes[i];
            }
        }
    }
    return height;
}

/**
* Run the Conjugent Gradient methog.
* Make the cart, initiate the variable based on the rank, run CG, print solution.
*
* nx - Number grid points on x axe.
* ny - Number grid points on y axe.
* c - number of iterations.
* eps - error.
* rank - rank of process.
* nrpr - Number of processes.
*
* return residual.
*/
void Expr_CG(int nx,int ny,int c,double eps, int py, int px, int nbrs[], int coords[], int cartrank, int tau, int k, double alfa1, Vector & uOld){
	//MPI topology      
	  /*int py = getHeight(nrpr);
      int px = nrpr/py;
      if(rank == 0){
        std::cout<<"Height="<<py<<", Width="<<px<<";\n";
      }
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
	  std::swap(coords[0],coords[1]);	*/
	
	//initialization	
    /*
    * If an process has a neighbour then it will have an goast layer in that drections,
    * otherwise the goast layer will be the boundary of the grid, so we need to find the begining of our grid
    * it can be (1,1) for the process that have top and left neighbours, (1,0) it we don't have the process 
    * from top, and (0,1) if we don't have process from left, and (0,0) if we are the top-left process.
    */
    //return;
    //#inner points on axes in each process(need for f initialization)
	int nx_ = (nx+1)/px , ny_ = (ny+1)/py; 
    //#inner points on axes for stencil calculation(need for u initialization)
	int nxi = nx_, nyi = ny_;
    //total size of a block
	int nxs = nx_ , nys = ny_;
    //first inner point coordonates
	int x=0,y=0; 
	int offsetx=coords[1]*nx_, offsety =coords[0]*ny_;
	
	if(coords[0]==0 || coords[0]==py-1){//if we are last or first row
	    if(coords[0]==py-1)
			ny_ += (ny+1)%py; // add the remaining rows
		if(py-1 != 0){//if  height of our cart is not 1
            nys = ny_+1;
        }
	}else nys +=2; //add the layers
	
	if(coords[1]==0 || coords[1]==px-1){//if we are last or first column
	    if(coords[1]==px-1)
			nx_ += (nx+1)%px; // add the remaining columns
        if(px != 1)// if width of our cart is not 1 
		   nxs = nx_+1;
	}else nxs +=2; //add the layers
	
	//first point that are initiated in f matrix
	if(coords[0]==0 && coords[1]==0){//first point
	    x = 0;
        y = 0;		
	}else if(coords[0]==0){//first column
	   x=1;y=0;
	}else if(coords[1]==0)//first row
	{
	   x = 0;
       y = 1;
	}else{
		x = 1;
        y = 1;
	}//others
	
    //points that are updated
	nxi = nx_;nyi = ny_;
	if(coords[0]==0 || coords[0]==py-1){//first row or last row
	    --nyi;	
        if(py == 1)--nyi;
	}
	if(coords[1]==0 || coords[1]==px-1){//first column
	    --nxi;
        if(px == 1)--nxi;
	}
	
	
	//number of inner grid points
	int pg = nxs*nys;	
	//int n_,int nx_,int ny_,int  nxi,int nyi,int nxs_,int nys_,int x_,int y_,size_t offsetx,size_t offsety
    //initializate the vector with 0.
	Vector  r(pg,nx_,ny_,nxi,nyi,nxs,nys,x,y),
			d(pg,nx_,ny_,nxi,nyi,nxs,nys,x,y),
			z(pg,nx_,ny_,nxi,nyi,nxs,nys,x,y); 
	//double pi = 3.141592653589793;
    //distance between grid points
	//double hx_ = 1.0/nx;
    //double hy_ = 1.0/ny;
    //coefficients for function initialization 
	//double C = 4*pi*pi;   
	//double freqx = 2*pi*hx_;   
	//double freqy = 2*pi*hy_;     
	
    /*
    * For the processes from the border top and left we initiate base on the fist inner point.
    */
	//Vector f(pg, nx_, ny_, nxi, nyi, nxs, nys , x, y,
	//		 [C,freqx,freqy,offsetx,offsety](int x1,int y1)->double{return C*sin(freqx*(x1+offsetx))*sinh(freqy*(y1+offsety));});	//4π^2 sin(2πx) sinh(2πy)
    
    //double SINH = sinh(2*pi); 	
	//sin(2πx) sinh(2πy)
	/*Vector *ut;
	if(coords[0]==py-1){//we are at border's block that are not 0
	     ut= new Vector(pg,nx_,ny_,nxi,nyi,nxs,nys,x,y, [SINH,freqx,offsetx](int x1)->double{return sin((x1+offsetx)*freqx) * SINH;});
    }
	else{
        //set vector with 0.
        ut = new Vector(pg,nx_,ny_,nxi,nyi,nxs,nys,x,y);
    }*/
    /*
    * When we update our points, we always begin with the point (1,1). Because we took care about right f initialization.
    */
	Vector u(pg,nx_,ny_,nxi,nyi,nxs,nys,x,y);
	//u(1,1)=1;
	Vector f(pg,nx_,ny_,nxi,nyi,nxs,nys,x,y);
	f = uOld;
	u.set_offset(offsetx,offsety);
	d.set_offset(offsetx,offsety);
	//u = *ut;
	
	Stencil A(nx, ny, nbrs, tau, k, alfa1);
	double delta0 = 0, delta1 = 0, beta = 0, alfa = 0 ,scalar = 0;
	//initialization end
	double eps2 = eps * eps;
	int nrIteration = 0;

	
	//CG
	r = f - A * u;
	delta0 = r^r;
	d = r;
	if((delta0 / ((nx-1)*(ny-1)))  > eps2 || eps <= 0){
	  for(int i = 0 ;i < c; ++i){
          ++nrIteration;
		  z = A*d;		  
		  scalar = (d ^ z);
		  alfa = delta0 / scalar;
		  u = u + d * alfa;
		  r = r - z * alfa;
		  delta1 = r^r;
		  beta = delta1 / delta0;
		  delta0 = delta1;
		  if(eps > 0 && (delta1 / ((nx-1)*(ny-1))) < eps2)break;		  
		  d = r + d * beta;		  
	  }
    }
    delta0 = r.LNorm() / ((nx-1)*(ny-1));
    if(cartrank == 0){
       std::cout<<"Number of iteration:"<<nrIteration<<";\n"
       			<<"Residual:"<<delta0<<"\n";
    }
    uOld = u;
	//return u;
	//print
	/*
	std::ofstream out("solution.txt");
    if(out.is_open() && nx != 10000 && nx != 10000) {
        if(cartrank == 0) {
            std::cout<<"writing output\n";
        }
        for( int i = 0; i < px * py; ++i )
              {
                 if( cartrank == 0 )
                 {
                    if( i == 0 )
                    {//if we are proccess 0 we write our solution.
                       u.set_distace(hx_,hy_);
                       u.set_offset(offsetx,offsety);	
                       out&u;
                    }
                    else
                    { //we receive data from other processes.
                       //out<<i<<'\n';	
                       Vector * pr = get_vector(i,nx,ny,px,py);	
                       pr->get_info(i);	
                       pr->set_distace(hx_,hy_);		
                       out&(*pr);	


                    }
                 }
                 else if( cartrank == i )
                 { // we are not process 0, then we sent data to process 0.
                       u.sent_info();	
                 }
              }
    }
    else {
      if(cartrank == 0) {    
        std::cout<<"No file writing.\n";
      }
    }
	*/	
}
