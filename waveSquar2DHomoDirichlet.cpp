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
void create_initial_v(int, double *);;
void printToFile(int, char *, double *);

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
  create_initial_v(N, v_prev);
  create_initial_v(N, v_now);

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
  double c_damp = 2/(2+b*dt);
  double c_prev = b*dt/(2+b*dt);

  //Create/fill the v_new vector/matrix 
  //Denne for-loopen gaar kun igjennom de indre punktene
  for(int i = 1; i < Ny; i++){
    for(int j = 1; j < Nx; j++){
      //double temp0 = C2*(v_now[i*(N)+j+1] + v_now[i*(N)+j-1] + v_now[(i-1)*(N)+j] + v_now[(i+1)*(N)+j] - 4*v_now[i*(N)+j]);
      //double temp1 = 2*v_now[i*(N)+j] - v_prev[i*(N)+j]; 
      v_next[i*(Nx+1)+j] = temp0+temp1;
      }
  }
  //Updateing the vectors/matrises(Change the pointer)						\
  temp_pointer = v_prev;
  v_prev = v_now;
  v_now  = v_next;
  v_next = temp_pointer;
}



//Creats the initial condition
void create_initial_v(int N, double *v_now)
{
  double h = 1.0/(N-1);
  for(int i = 0; i < N; i++){
    for(int j = 0; j < N; j++){
      if(i == 0 || i == N || j == 0 || j == N){
	v_now[i*(N)+j] = 0;
      }
      v_now[i*(N)+j] = exp(-40*(pow((2*(i*h)-1.4),2) + pow((2*(j*h)-1),2)));
    }
  }

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
