/*
  Auther: Haakon Osterbo og Jarle Sogn


*/



#include <cmath>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>


using namespace std;

ofstream ofile;

void interate_v(int,int, double, double, double, double, double *,double *,double *);
void create_initial_v(int, int, double, double, double *, double *);
void printToFile(int, char *, double *);


double c(double, double);
double f(double, double);
double I(double, double);
double V(double, double);

int main(int argc, char* argv[])
{ 
  //Input variabels(read form commandline)
  int Nx; int Ny; int M; double T; double c;

  if( argc <= 5 ){
    cout << "Bad Usage: " << argv[0] << 
      "Read also in: Nx, Ny, M, T, and c, on same line" << endl;
    exit(1);
  }
  else{
    Nx = atoi(argv[1]); Ny = atoi(argv[2]);  M = atoi(argv[3]); 
    T = atof(argv[4]) ; c  = atof(argv[5]);
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
  double dx = 1.0/Nx;
  double dy = 1.0/Ny;
  double dt= T/M;
  double *temp_pointer;
  char outfilename[60]; 
  int ff   = 1; //Frame frevense
  if(M >=(int)T*24){ff = M/((int)T*24);}
  

  //Create the initial condition
  create_initial_v(Nx, Ny, dx, dy, v_now, V_prev);

  //Write the IC to a file 
  sprintf(outfilename, "wave_squar_2D_N%d_M%d_t%4.2f.dat", N, M, 0*dt);
  printToFile(N, outfilename, v_now);

  //Main Loop:
  //For each interation it move one timestep dt forward
  for (int i=1; i <= M;i++){
    interate_v(Nx, Ny, dx, dy, dt, b, v_next,v_now,v_prev);
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
      sprintf(outfilename, "wave_squar_2D_N%d_M%d_t%4.2f.dat", N, M, i*dt);
      printToFile(N, outfilename, v_now);}
  }

  //Clean up
  delete [] v_prev;
  delete [] v_now ;
  delete [] v_next;
  cout << "TEST ------ MAIN ALL DONE!!!"<< endl;
  return 0;
}



//Moves one time step forward
void interate_v(int Nx, ,int Ny, double dx, double dy, double dt, double b,double* v_next, double* v_now, double* v_prev)
{
  double *temp_pointer;
  //Div constants to save flops
  double cx_tmp = 2*dt*dt/((2+b*dt)*dx*dx);
  double cx_tmp = 2*dt*dt/((2+b*dt)*dy*dy);
  double cf_tmp = (2 + b*dt)/(2*dt*dt);
  double c_damp = -2/(2+b*dt);
  double c_prev = b*dt/(2+b*dt);

  //Create/fill the v_new vector/matrix 
  //Denne for-loopen gaar kun igjennom de indre punktene
  for(int i = 1; i < Ny; i++){
    for(int j = 1; j < Nx; j++){
      double temp0 = cx_tmp*(c(i*dx+dx/2,j*dy)*(v_now[(i+1)*(Nx+1)+j] - v_now[i*(Nx+1)+j]) - c(i*dx-dx/2,j*dy)*(v_now[i*(Nx+1)+j] - v_now[(i-1)*(Nx+1)+j]));
      double temp1 = cy_tmp*(c(i*dx,j*dy+dy/2)*(v_now[i*(Nx+1)+j+1]   - v_now[i*(Nx+1)+j]) - c(i*dx,j*dy-dy/2)*(v_now[i*(Nx+1)+j] - v_now[i*(Nx+1)+j-1]))  ;
      double temp3 = cf_tmp*f(i*dx,j*dy) + c_prev*v_prev[i*(N+1)+j] + c_damp*(v_prev[i*(N+1)+j] - 2*v_now[i*(N+1)+j]);
	v_next[i*(Nx+1)+j] = temp0+temp1+temp3;
      }
  }
  //Updateing the vectors/matrises(Change the pointer)

  temp_pointer = v_prev;
  v_prev = v_now;
  v_now  = v_next;
  v_next = temp_pointer;
}

//Sets Neuman boundary condition
void neuman_boundary_cond(int Nx, int Ny, double dx, double dy, double b, double* v_next, double* v_now, double* v_prev)
{
        double cx_tmp = 2*dt*dt/((2+b*dt)*dx*dx);
        double cy_tmp = 2*dt*dt/((2+b*dt)*dy*dy);
        double cf_tmp = (2 + b*dt)/(2*dt*dt);
        double c_damp = -2/(2+b*dt);
        double c_prev = b*dt/(2+b*dt);

        //Boundary conditions for x
        for(int j = 1; j<Ny; j++){
                //x = 0 boundary
                double temp0 = cx_tmp*(v_now[1*(Nx+1)+j]-v_now[0*(Nx+1)+j])*(c(0.5*dx,j*dy) + c(-0.5*dx,j*dy));
                double temp1 = cy_tmp*(c(0,j*dy+0.5*dy)*(v_now[0*(Nx+1)+j+1]-v_now[0*(Nx+1)+j]) - c(0,j*dy-0.5*dy)*(v_now[0*(Nx+1)+j]-v_now[0*(Nx+1)+j-1]));
                double temp3 = cf_tmp*f(0,j*dy) + c_prev*v_prev[0*(N+1)+j] + c_damp*(v_prev[0*(N+1)+j] - 2*v_now[0*(N+1)+j]);
                v_next[0*(Nx+1)+j] = temp0+temp1+temp3;

                //x = Lx boundary
                double temp0 = cx_tmp*(v_now[(Nx-1)*(Nx+1)+j]-v_now[Nx*(Nx+1)+j])*(c(Nx*dx + 0.5*dx,j*dy) + c(Nx*dx - 0.5*dx,j*dy));
                double temp1 = cy_tmp*(c(Nx*dx,j*dy+0.5*dy)*(v_now[Nx*(Nx+1)+j+1]-v_now[Nx*(Nx+1)+j]) - c(Nx*dx,j*dy-0.5*dy)*(v_now[Nx*(Nx+1)+j]-v_now[Nx*(Nx+1)+j-1]));
                double temp3 = cf_tmp*f(Nx*dx,j*dy) + c_prev*v_prev[Nx*(N+1)+j] + c_damp*(v_prev[Nx*(N+1)+j] - 2*v_now[Nx*(N+1)+j]);
                v_next[Nx*(Nx+1)+j] = temp0+temp1+temp3;
        }

        for(int i = 1; i<Nx; i++){
                //y = 0 boundary
                double temp1 = cx_tmp*(c(i*dx+.5*dx,0)*(v_now[(i+1)*(Nx+1)+0]-v_now[i*(Nx+1)+0]) - c(i*dx-.5*dx,0)*(v_now[i*(Nx+1)+0]-v_now[(i-1)*(Nx+1)+0]));
                double temp0 = cy_tmp*(v_now[i*(Nx+1)+1]-v_now[i*(Nx+1)+0])*(c(i*dx,.5*dy) + c(i*dx,-.5*dy));
                double temp3 = cf_tmp*f(i*dx,0) + c_prev*v_prev[i*(N+1)+0] + c_damp*(v_prev[i*(N+1)+0] - 2*v_now[i*(N+1)+0]);
                v_next[i*(Nx+1)+0] = temp0+temp1+temp3;

                //y = Ly boundary
                double temp1 = cx_tmp*(c(i*dx+.5*dx,Ny*dy)*(v_now[(i+1)*(Nx+1)+Ny]-v_now[i*(Nx+1)+Ny]) - c(i*dx-.5*dx,Ny*dy)*(v_now[i*(Nx+1)+Ny]-v_now[(i-1)*(Nx+1)+Ny]));
                double temp0 = cy_tmp*(v_now[i*(Nx+1)+Ny-1]-v_now[i*(Nx+1)+Ny])*(c(i*dx,Ny*dy+.5*dy) + c(i*dx,Ny*dy-.5*dy));
                double temp3 = cf_tmp*f(i*dx,Ny*dy) + c_prev*v_prev[i*(N+1)+Ny] + c_damp*(v_prev[i*(N+1)+Ny] - 2*v_now[i*(N+1)+Ny]);
                v_next[i*(Nx+1)+Ny] = temp0+temp1+temp3;
        }
}

//Creats the initial condition
void create_initial_v(int Nx, int Ny, double dx, double dy,double *v_now, double *v_prev)
{
  //u(x,y,t=0), n = 0
  for(int i = 0; i < Ny+1; i++){
    for(int j = 0; j < Nx+1; j++){
      v_prev[i*(Nx+1)+j] = I(dx*i,dy*j);
    }}
  //u(x,y,t=dt), n = 1
  for(int i = 0; i < Ny+1; i++){
    for(int j = 0; j < Nx+1; j++){
      v_now[i*(Nx+1)+j] = 0.0;
    }}
}


//Prints data to file
void printToFile(int N, char *outfilename, double *v)
{
  ofile.open(outfilename); //Open outputfile;
  for(int i = 1; i < N-1; i++){
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    for(int j = 1; j < N-1; j++){
      ofile << setw(15) << setprecision(8) << v[i*(N)+j];
    }
    ofile << endl;
  }
  ofile.close(); //Close outupfile
}



double c(double x, double y );
{
return 1.0
}


double f(double x, double y );
{
return 0.0
}


double I(double, double);
{
return 1.0
}

double V(double, double);
{
return 0.0
}
