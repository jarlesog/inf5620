/*
  Auther: Haakon Osterbo og Jarle Sogn


*/



#include <cmath>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>


using namespace std;
class waveFunctions
{
        public:
                waveFunctions(double Lx, double Ly, double dx, double dy, double dt);
                double c(double x, double y);
                double f(double x, double y);
                double I(double x, double y);
                double V(double x, double y);
                double getCFL();
                bool isStabile();
        private:
                double Lx; double Ly;
                double dx; double dy; double dt;
                double max_c;
};

ofstream ofile;

void interate_v(int,int, double, double, double, double, waveFunctions, double *,double *,double *);
void create_initial_v(int, int, double, double, double, double, waveFunctions, double *, double *);
void neuman_boundary_cond(int, int, double, double, double, double, waveFunctions, double*, double*, double*);
void printToFile(int, int, char *, double *);


double c(double, double);
double f(double, double);
double I(double, double, double, double);
double V(double, double);

int main(int argc, char* argv[])
{ 
  //Input variabels(read form commandline)
  int Nx; int Ny; int M; double T; double Lx; double Ly;

  if( argc <= 6 ){
    cout << "Bad Usage: " << argv[0] << "Read also in: Nx, Ny, M, T, Lx, Ly on same line" << endl;
    exit(1);
  }
  else{
    Nx = atoi(argv[1]); Ny = atoi(argv[2]);  M = atoi(argv[3]); 
    T = atof(argv[4]); Lx = atof(argv[5]); Ly = atof(argv[6]);
  }
  
  
  //Reserving space for my vactor(matrises)
  double *v_prev = new double [(Nx+1)*(Ny+1)];
  double *v_now  = new double [(Nx+1)*(Ny+1)];
  double *v_next = new double [(Nx+1)*(Ny+1)];
  
  //#Alternative method
  //#double** matrix = new double[N];
  //# for(i = 0; i<N;++1)
  //#matrix[i]= buf + i*N //buf er en vector.
  
  //Creating contans
  double dx = Lx/Nx;
  double dy = Ly/Ny;
  double dt= T/M;
  
  waveFunctions w(Lx, Ly, dx, dy, dt);
  
  if(!w.isStabile())
  {
          cout << "WARNING: CFL condition is not satisfied. c_max*dt*sqrt(1/(dx*dx) + 1/(dy*dy)) = " << w.getCFL() << endl;
  }

  double b = 1.0;
  double *temp_pointer;
  char outfilename[60]; 
  int ff   = 1; //Frame frevense
  if(M >=(int)T*24){ff = M/((int)(T*24));}
  

  //Create the initial condition
  create_initial_v(Nx, Ny, dx, dy, Lx, Ly, w, v_now, v_prev);

  //Write the IC to a file 
  sprintf(outfilename, "wave_squar_2D_Nx%d_Ny%d_M%d_t%4.2f.dat", Nx, Ny, M, 0*dt);
  printToFile(Nx,Ny, outfilename, v_now);

  //Main Loop:
  //For each interation it move one timestep dt forward
  for (int i=1; i <= M;i++){
    interate_v(Nx, Ny, dx, dy, dt, b,w, v_next,v_now,v_prev);
    neuman_boundary_cond(Nx, Ny, dx, dy, dt, b,w, v_next, v_now, v_prev);
    //updateing pointers
    temp_pointer = v_prev;
    v_prev = v_now;
    v_now  = v_next;
    v_next = temp_pointer;
    
    
    //#2byte per tall som blir laget

    //I don't need to make plot of all the interations
    //so i use 24 frames per seconds(to save run-time 
    //and space.
    if(i%ff == 0){
      sprintf(outfilename, "wave_squar_2D_Nx%d_Ny%d_M%d_t%4.2f.dat", Nx, Ny,M, i*dt);
      printToFile(Nx, Ny, outfilename, v_now);}
  }

  //Clean up
  delete [] v_prev;
  delete [] v_now ;
  delete [] v_next;
  cout << "TEST ------ MAIN ALL DONE!!!"<< endl;
  return 0;
}



//Moves one time step forward
void interate_v(int Nx, int Ny, double dx, double dy, double dt, double b, waveFunctions w, double* v_next, double* v_now, double* v_prev)
{
  double *temp_pointer;
  //Div constants to save flops
  double cx_tmp = 2*dt*dt/((2+b*dt)*dx*dx);
  double cy_tmp = 2*dt*dt/((2+b*dt)*dy*dy);
  double cf_tmp = (2 + b*dt)/(2*dt*dt);
  double c_damp = -2/(2+b*dt);
  double c_prev = b*dt/(2+b*dt);
  double temp0, temp1, temp2;
  //Create/fill the v_new vector/matrix 
  //Denne for-loopen gaar kun igjennom de indre punktene
  for(int i = 1; i < Ny; i++){// i is y axis
    for(int j = 1; j < Nx; j++){//j is x axis
      temp0 = cx_tmp*(w.c(j*dx+dx/2,i*dy)*(v_now[(i+1)*(Nx+1)+j] - v_now[i*(Nx+1)+j]) - w.c(j*dx-dx/2,i*dy)*(v_now[i*(Nx+1)+j] - v_now[(i-1)*(Nx+1)+j]));
      temp1 = cy_tmp*(w.c(j*dx,i*dy+dy/2)*(v_now[i*(Nx+1)+j+1]   - v_now[i*(Nx+1)+j]) - w.c(j*dx,i*dy-dy/2)*(v_now[i*(Nx+1)+j] - v_now[i*(Nx+1)+j-1]))  ;
      temp2 = cf_tmp*w.f(j*dx,i*dy) + c_prev*v_prev[i*(Nx+1)+j] + c_damp*(v_prev[i*(Nx+1)+j] - 2*v_now[i*(Nx+1)+j]);
	v_next[i*(Nx+1)+j] = temp0+temp1+temp2;
      }
  }
  //Updateing the vectors/matrises(Change the pointer)
  //temp_pointer = v_prev;
  //v_prev = v_now;
  //v_now  = v_next;
  //v_next = temp_pointer;
}

//Sets Neuman boundary condition
void neuman_boundary_cond(int Nx, int Ny, double dx, double dy, double dt, double b,waveFunctions w, double* v_next, double* v_now, double* v_prev)
{
        double cx_tmp = 2*dt*dt/((2+b*dt)*dx*dx);
        double cy_tmp = 2*dt*dt/((2+b*dt)*dy*dy);
        double cf_tmp = (2 + b*dt)/(2*dt*dt);
        double c_damp = -2/(2+b*dt);
        double c_prev = b*dt/(2+b*dt);
	double temp0, temp1, temp2;
        
        //Boundary conditions for y
        for(int j = 1; j<Nx; j++){
	  //y = 0 boundary;
	  temp0 = cy_tmp*(v_now[1*(Nx+1)+j]-v_now[0*(Nx+1)+j])*(w.c(j*dx,.5*dy) + w.c(j*dx,-.5*dy));
	  temp1 = cx_tmp*(w.c(j*dx+0.5*dx,0)*(v_now[0*(Nx+1)+j+1]-v_now[0*(Nx+1)+j]) - w.c(j*dx-0.5*dx,0)*(v_now[0*(Nx+1)+j]-v_now[0*(Nx+1)+j-1]));
	  temp2 = cf_tmp*w.f(j*dx,0) + c_prev*v_prev[0*(Nx+1)+j] + c_damp*(v_prev[0*(Nx+1)+j] - 2*v_now[0*(Nx+1)+j]);
	  v_next[0*(Nx+1)+j] = temp0+temp1+temp2;
	  
	  //y = Ly boundary
	  temp0 = cy_tmp*(v_now[(Nx-1)*(Nx+1)+j]-v_now[Nx*(Nx+1)+j])*(w.c(j*dx,Ny*dy + 0.5*dy) + w.c(j*dx,Ny*dy - 0.5*dy));
	  temp1 = cx_tmp*(w.c(j*dx+.5*dx,Ny*dy)*(v_now[Nx*(Nx+1)+j+1]-v_now[Nx*(Nx+1)+j]) - w.c(j*dx-.5*dx,Ny*dy)*(v_now[Nx*(Nx+1)+j]-v_now[Nx*(Nx+1)+j-1]));
	  temp2 = cf_tmp*w.f(j*dx,Ny*dy) + c_prev*v_prev[Nx*(Nx+1)+j] + c_damp*(v_prev[Nx*(Nx+1)+j] - 2*v_now[Nx*(Nx+1)+j]);
	  v_next[Nx*(Nx+1)+j] = temp0+temp1+temp2;
        }

        //Boundary conditions for x
        for(int i = 1; i<Ny; i++){
	  //x = 0 boundary
	  temp0 = cy_tmp*(w.c(0,i*dy+.5*dy)*(v_now[(i+1)*(Nx+1)+0]-v_now[i*(Nx+1)+0]) - w.c(0,i*dy-.5*dy)*(v_now[i*(Nx+1)+0]-v_now[(i-1)*(Nx+1)+0]));
	  temp1 = cx_tmp*(v_now[i*(Nx+1)+1]-v_now[i*(Nx+1)+0])*(w.c(.5*dx,i*dy) + w.c(-.5*dx,i*dy));
	  temp2 = cf_tmp*w.f(0,i*dy) + c_prev*v_prev[i*(Nx+1)+0] + c_damp*(v_prev[i*(Nx+1)+0] - 2*v_now[i*(Nx+1)+0]);
	  v_next[i*(Nx+1)+0] = temp0+temp1+temp2;
	  
	  //x = Lx boundary
	  temp0 = cy_tmp*(w.c(Nx*dx,i*dy+.5*dy)*(v_now[(i+1)*(Nx+1)+Ny]-v_now[i*(Nx+1)+Ny]) - w.c(Nx*dx,i*dy-.5*dy)*(v_now[i*(Nx+1)+Ny]-v_now[(i-1)*(Nx+1)+Ny]));
	  temp1 = cx_tmp*(v_now[i*(Nx+1)+Ny-1]-v_now[i*(Nx+1)+Ny])*(w.c(Nx*dx+.5*dx,i*dy) + w.c(Nx*dx-.5*dx,i*dy));
	  temp2 = cf_tmp*w.f(Nx*dx,i*dy) + c_prev*v_prev[i*(Nx+1)+Ny] + c_damp*(v_prev[i*(Nx+1)+Ny] - 2*v_now[i*(Nx+1)+Ny]);
	  v_next[i*(Nx+1)+Ny] = temp0+temp1+temp2;
        }
	//disse er implimentert rast, mulig det er noen feil her
	//Manualy taking the corner x = 0, y = 0
	temp0 = cy_tmp*(v_now[1*(Nx+1)+0]-v_now[0*(Nx+1)+0])*(w.c(0,.5*dy) + w.c(0,-.5*dy));
	temp1 = cx_tmp*(v_now[0*(Nx+1)+1]-v_now[0*(Nx+1)+0])*(w.c(.5*dx,0) + w.c(-.5*dx,0));
	temp2 = cf_tmp*w.f(0,0) + c_prev*v_prev[0*(Nx+1)+0] + c_damp*(v_prev[0*(Nx+1)+0] - 2*v_now[0*(Nx+1)+0]);
	v_next[0*(Nx+1)+0] = temp0+temp1+temp2;
	//Manualy taking the corner x = Nx, y = Ny
	temp0 = cy_tmp*(v_now[(Ny-1)*(Nx+1)+Nx]-v_now[Ny*(Nx+1)+NX])*(w.c(Nx,Ny+.5*dy) + w.c(Nx,-.Ny+5*dy));
	temp1 = cx_tmp*(v_now[Ny*(Nx+1)+Nx-1]-v_now[Ny*(Nx+1)+Nx])*(w.c(Nx+.5*dx,Ny) + w.c(Nx-.5*dx,Ny));
	temp2 = cf_tmp*w.f(Nx,Ny) + c_prev*v_prev[Ny*(Nx+1)+Nx] + c_damp*(v_prev[Ny*(Nx+1)+Nx] - 2*v_now[Ny*(Nx+1)+Nx]);
	v_next[Ny*(Nx+1)+Nx] = temp0+temp1+temp2;
	//Manualy taking the corner x = Nx, y = 0
	temp0 = cy_tmp*(v_now[1*(Nx+1)+Nx]-v_now[0*(Nx+1)+Nx])*(w.c(Nx,.5*dy) + w.c(Nx,-.5*dy));
	temp1 = cx_tmp*(v_now[0*(Nx+1)+Nx-1]-v_now[0*(Nx+1)+Nx])*(w.c(Nx+.5*dx,0) + w.c(Nx-.5*dx,0));
	temp2 = cf_tmp*w.f(Nx,0) + c_prev*v_prev[0*(Nx+1)+Nx] + c_damp*(v_prev[0*(Nx+1)+Nx] - 2*v_now[0*(Nx+1)+Nx]);
	v_next[0*(Nx+1)+Nx] = temp0+temp1+temp2;	
	//Manualy taking the corner x = 0, y = Ny
	temp0 = cy_tmp*(v_now[(Ny-1)*(Nx+1)+0]-v_now[Ny*(Nx+1)+0])*(w.c(0,Ny+.5*dy) + w.c(0,Ny-.5*dy));
	temp1 = cx_tmp*(v_now[Ny*(Nx+1)+1]-v_now[Ny*(Nx+1)+0])*(w.c(.5*dx,Ny) + w.c(-.5*dx,Ny));
	temp2 = cf_tmp*w.f(0,Ny) + c_prev*v_prev[Ny*(Nx+1)+0] + c_damp*(v_prev[Ny*(Nx+1)+0] - 2*v_now[Ny*(Nx+1)+0]);
	v_next[Ny*(Nx+1)+0] = temp0+temp1+temp2;

	
}

//Creats the initial condition
void create_initial_v(int Nx, int Ny, double dx, double dy, double Lx, double Ly,waveFunctions w, double *v_now, double *v_prev)
{
  //u(x,y,t=0), n = 0
  for(int i = 0; i < Ny+1; i++){
    for(int j = 0; j < Nx+1; j++){
      v_now[i*(Nx+1)+j] = w.I(dx*j,dy*i);
      v_prev[i*(Nx+1)+j] = v_now[i*(Nx+1)+j] - w.V(j*dx,i*dy);//Backward Euler
    }}
}


//Prints data to file
void printToFile(int Nx, int Ny, char *outfilename, double *v)
{
  ofile.open(outfilename); //Open outputfile;
  //Prints as matrix A_{i,j}
  for(int i = 0; i <= Ny; i++){
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    for(int j = 0; j <= Nx; j++){
      ofile << setw(15) << setprecision(8) << v[i*(Nx+1)+j];
    }
    ofile << endl;
  }
  ofile.close(); //Close outupfile
}


waveFunctions::waveFunctions(double Lxin, double Lyin, double dxin, double dyin, double dtin)
{
        Lx = Lxin; Ly = Lyin; 
        dx = dxin; dy = dyin; dt = dtin;
        max_c = 1;
}

double waveFunctions::c(double x, double y )
{
        return 1.0;
}

double waveFunctions::f(double x, double y )
{
        return 0.0;
}

double waveFunctions::I(double x, double y)
{
        double a = 20;
        return exp(-a*((x-0.5*Lx)*(x-0.5*Lx) + (y-0.5*Ly)*(y-0.5*Ly)));
}

double waveFunctions::V(double x, double y)
{
        return 0.0;
}

double waveFunctions::getCFL()
{
        return max_c*dt*sqrt(1/(dx*dx) + 1/(dy*dy));
}

bool waveFunctions::isStabile()
{
        return getCFL() <= 1;
}



