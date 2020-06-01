//A celluar automaton model to simulate tumor growth in d-dimensions
//Using a sequential implementation, with sophiscated data structure to maximize efficiency
//Complementary to the paralell version by Jana

//author: Yang JIAO, yjiao@princeton.edu
//started: Oct. 27, 2010



/* Input file required for algorithm: input.dat							*/
/*  - This file has parameters for algorithm and input files 					*/
/*      and specify the name of the other three input files and output file name                */
/*  - Input file 1: inputPoints: takes list of cell centers sorted from furthest from center 	*/
/*      to closest										*/
/*  - Input file 2: inputArea: area of each automaton cell					*/
/*  - Input file 3: inputNeighbors: neighbors of each automaton cell				*/



//Novel Features of this method:
//(1) We use a grain-corasening idea to check cells on the boundary, i.e., first divide the whole 
//    domain into small grids, then only make a list of grid that contains boundary cells;
//    when computing distance between bd cells and any particular cell, first compute the 
//    distance between the two grids, when the closest grid is found, then compute the bd cells.
//
//(2) For the tumor boundary, it is dynamic, using grid requries frequent update, not efficient
//    Instead, we use a local-rule, i.e., each proliferative has a closest-edge cell, if this 
//    cell becomes tumor cell, all its neighbors are searched for new healthy cell or bd cell
//    The neighbor search can be 2nd order or even higher order in depth, meaning searching 
//    the neighbors of a neighbor cell. But for now, we just search the first neighbors.
//  !!!!!- Note: this neighbor-search-method requires update quiscent cells first...
//  !!!!!- It turns out that (2) would lead to shape not closely follow the bd...NOT USING THIS...

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

using namespace std;

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <map>
#include <string>
#include <time.h>
#include <vector>

#define NAME_LEN 256 //for the name of the input and output files....
//#define RAND_MAX 10000000 //denominator for maximum number
#define TOL 0.00000001 //a numerical tolerance...

int n_dim; //the spatial dimension
int n_grid; //the linear size of grid along each direction
int n_neig; //this is not for read in, computed as 3^n_dim, not a necessary num.


//Tumor algorithm variables
double a, b, p0, Rmax; //the four minimalist parameters Rmax is not defined here 
double length; //the length of a square domain 
int iterations, t1, t2, t3, t4, t5; //total number of iterations and time points for output
int bd_type; //specify the boundary shape and type (whether can provide nuritiion or not)


//names of input files
char dataPoints[NAME_LEN];	
char dataVol[NAME_LEN];	
char dataNeighbors[NAME_LEN];
//name of the output file for tumor radius and volume
//char outRadVol[NAME_LEN];

//information for cells
vector< vector<double> > cell_coord; //the coordinates of the cells
vector<int> cell_type; //-1 healthy outside bd; 0 - healthy within bd; 1 - proliferative; 2 - quiscent; 3 - necrotic
vector<int> cell_bd_flag; //whehter this cell is on the boundary
vector<int> cell_edge_flag; //whether this cell is on the tumor edge
vector<int> cell_close_edge_index; //the index of the closest tumor edge cells to the current cell, obtained by search neighors...

vector<double> cell_vol; //volume of each cell
vector< vector<int> > cell_neig; //the neighbor list of each cell
vector<int> cell_neig_ct; //the counter for the neighbors

int n_cell; //the total number of cells (including those within and outside the growth region)
double cell_len; //a characteristic length of cells


//informtion of for the tumor
vector<double> tumor_ct_coord; //the coordinates of the tumor center

vector<int> necrotic_list; //the list of cells associated with tumors, constantly updated
vector<int> quiscent_list;
vector<int> prolif_list;

int n_nec, n_qui, n_pro, n_tumor; //numbers
double vol_nec, vol_qui, vol_pro; //volumes

//vector<int> healthy_edge_list; //the list of healthy cells at tumor edge, constantly updated
vector<int> healthy_bd_list; //the list of healthy cells at the boundary, static

//information of the grids, this is for the quick search of boundary cells...
int  n_tot_grid;

vector< vector<int> > bd_grid_list; //a list of grid in which each one has a list of index of boundary cells
vector<int> bd_grid_index; //the index of grid that contain boundary cells
vector< vector<double> > bd_grid_coord; //the coordinates of the grid center...
vector<int> bd_grid_ct; //the counter
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Init_Data()
{
  //read in the basic parameters...
  ifstream infile;
  infile.open("input.dat", ios::in);
  if(!infile)
    { 
      cout << "Can't open input.dat\n";
      exit(1);
    }

  infile >> n_dim; //the dimension
  infile >> n_grid; //num.of grid, not in the original input file
  n_neig = (int)pow(3.0, n_dim); //the maximum number of neighbors at this dimension 

  for(int i=0; i<n_dim; i++)
    {
      tumor_ct_coord.push_back(-1); //the order can be messed up using push_back for new data...
    }
  for(int i=0; i<n_dim; i++)
    {
      infile >> tumor_ct_coord[i]; //this is a saver way of doing this
      //at this stage, it is not the real tumor center, just a specified one...
    }
  

  infile >> a;  infile >> b;  infile >> p0;  infile >> Rmax;
  infile >> length;  infile >> iterations;
  infile >> t1; infile >> t2; infile >> t3;  infile >> t4; infile >> t5; 
  infile >> bd_type;
  infile >> dataPoints; infile >> dataVol;
  infile >> dataNeighbors; 
  //infile >> outRadVol;
  infile.close();

  //double check for correct input...
  cout<< ". n_dim = "<< n_dim <<endl;
  cout<< ". n_grd = "<<n_grid <<endl;
  //cout<<" x0, y0, z0"<<endl;
  cout << ". a = " << a << endl;
  cout << ". b = " << b << endl;
  cout << ". p0 = " << p0 << endl;
  cout << ". Rmax = " << Rmax << endl;
  cout << ". length = " << length << endl; //length of square simulation domain, containing the boundary
  cout << ". iterations = " << iterations << endl;
  cout << ". t1 = " << t1 << endl;
  cout << ". t2 = " << t2 << endl;
  cout << ". t3 = " << t3 << endl;
  cout << ". t4 = " << t4 << endl;
  cout << ". t5 = " << t5 << endl;
  cout << ". type = " << bd_type << endl; //specify the type of the boundary
  cout << ". dataPoints = " << dataPoints << endl;
  cout << ". dataArea = " << dataVol << endl;
  cout << ". dataNeighbors = " << dataNeighbors << endl;
  //cout << ". outRadArea = " << outRadVol << endl;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //now read in the confiugration data for the cells...
  ifstream inputPoints(dataPoints);
  if(!inputPoints)
    {
      cout << "Cannot open file " << dataPoints << "\n"; 
      exit(1);
    }
  ifstream inputVol(dataVol);
  if(!inputVol)
    {
      cout << "Cannot open file " << dataVol << "\n"; 
      exit(1);
    }
  ifstream inputNeighbors(dataNeighbors);	
  if(!inputNeighbors)
    {
      cout << "Cannot open file " << dataNeighbors << "\n"; 
      exit(1);
    }


  inputPoints >> n_cell; //read in cell numbers
  inputPoints >> cell_len; //read in cell length

  for(int i=0; i<n_cell; i++)
    {
      cell_coord.push_back(vector<double>(n_dim, -1.0));
      cell_type.push_back(0);
      cell_bd_flag.push_back(0); //on bd, 1, else 0
      cell_edge_flag.push_back(0);
      
      cell_close_edge_index.push_back(-1); //the index of the closest edge cell for each cell

      cell_vol.push_back(-1.0);
      
      cell_neig.push_back(vector<int>(n_neig, -1));
      cell_neig_ct.push_back(0);
    }

  //read in the coordinates and cell volumes...
  int temp;
  for(int i=0; i<n_cell; i++)
    {
      for(int j=0; j<n_dim; j++)
	inputPoints >> cell_coord[i][j];

      inputVol >> temp >> cell_vol[i];
      if(temp!=i)
	{
	  cout << "Problem reading area.dat file\n"; 
	  exit(1);
	}
    }
  inputPoints.close();
  inputVol.close();

  //read in the neighbor list
  int temp_num;
  for(int i=0; i<n_cell; i++)
    {
      inputNeighbors >> temp >> temp_num;

      if(i!=temp)
	{
	  cout << "problem reading neighbors file: since i = " << i << ", temp = " << temp << "\n"; 
	  exit(1); 
	}
      else if(temp_num>n_neig)
	{
	  cout<<" the number of neighbors is greater than the max, i.e., "<<temp_num<<" > "<<n_neig<<endl;
	  exit(1);
	}

      cell_neig_ct[i] = temp_num;
      for(int j=0; j<temp_num; j++)
	inputNeighbors >> cell_neig[i][j];
    }
  inputNeighbors.close();


  //now we initialize the grid...~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  n_tot_grid = (int)pow(n_grid, n_dim);
  int n_cell_per_grid = (int)floor((double)2*n_cell/(double)n_tot_grid);

  for(int i=0; i<n_tot_grid; i++)
    {
      bd_grid_list.push_back(vector<int> (n_cell_per_grid, -1)); //the index of bd cells
      bd_grid_ct.push_back(0); //number of bd_cells in each grid
      bd_grid_coord.push_back(vector<double> (n_dim, -1.0)); //the coordinates of the grids

    
    }

  for(int i=0; i<n_tot_grid; i++)
    {
      //now we compute the coordinates of the grid...
      if(n_dim == 2)
	{
	  double x = i%n_grid;
	  double y = floor((double)i/(double)n_grid);
	  
	  bd_grid_coord[i][0] = x;
	  bd_grid_coord[i][1] = y;
	}
      else if(n_dim == 3)
	{
	  double x = (i%(n_grid*n_grid))%n_grid;
	  double y = (int)floor((double)(i%(n_grid*n_grid))/(double)n_grid);
	  double z = (int)floor((double)i/((double)n_grid*n_grid));
	  
	  bd_grid_coord[i][0] = x;
	  bd_grid_coord[i][1] = y;
	  bd_grid_coord[i][2] = z;
	}
       else 
	{
	  cout<<" dim = "<<n_dim<<" is not programmed!"<<endl;
	  exit(1);
	}
    }
  
}


//get the shape and the growth permitting domain, 
//get the boundary 
//also the grid containing the list of boundary cells
void Init_Boundary(int type)
{
  //first get the boundary~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cout<<" Generating the growth permitting domain..."<<endl;

  if(n_dim == 2)
    {
      if(type==1)
	{ //complex environment is an ellipse: (x-x0)^2/ra^2 + (y-y0)^2/rb^2 = 1
	  double x0 = 0.5*length;
	  double y0 = 0.5*length;
	  double ra = 0.48;
	  double rb = 0.12;
	  for(int i = 0; i<n_cell; i++)
	    {
	      double xDist = (cell_coord[i][0]-x0)/ra;
	      double yDist = (cell_coord[i][1]-y0)/rb;
	      double curDist = (xDist*xDist)+(yDist*yDist);
	      if(curDist>=1)
		cell_type[i] = -1;
	    }		
	}
      else if(type==2)
	{ //complex environment is 3 of 4 quadrants of the unit square, minus a small rim about the square boundary
	  double midpt = 0.5*length;

	  for(int i = 0; i<n_cell; i++)
	    {
	      if((cell_coord[i][0]>=midpt)&&(cell_coord[i][1]>=midpt)) cell_type[i] = -1;	
	      else if(cell_coord[i][0]<0.1*midpt) cell_type[i] = -1;
	      else if(cell_coord[i][0]>1-0.1*midpt) cell_type[i] = -1;
	      else if(cell_coord[i][1]<0.1*midpt) cell_type[i] = -1;
	      else if(cell_coord[i][1]>1-0.1*midpt) cell_type[i] = -1;		
	    }
	}
      else
	{
	  cout<<" type == "<<type<<" hasn't programmed yet..."<<endl;
	  exit(1);
	}
    }
  else if(n_dim == 3)
    {
      if(type==1)
	{ //complex environment is an ellipsoid: (x-x0)^2/ra^2 + (y-y0)^2/rb^2 + (z-z0)^2/rc^2= 1
	  double x0 = 0.5*length;
	  double y0 = 0.5*length;
	  double z0 = 0.5*length;
	  double ra = 0.48;
	  double rb = 0.3184;
	  double rc = 0.28;
	  for(int i = 0; i<n_cell; i++)
	    {
	      double xDist = (cell_coord[i][0]-x0)/ra;
	      double yDist = (cell_coord[i][1]-y0)/rb;
	      double zDist = (cell_coord[i][2]-z0)/rc;
	      double curDist = (xDist*xDist)+(yDist*yDist)+(zDist*zDist);
	      if(curDist>=1)
		cell_type[i] = -1;
	    }		
	}
      else
	{
	  cout<<" type == "<<type<<" hasn't programmed yet..."<<endl;
	  exit(1);
	}
    }
  else
    {
      cout<<" boundary for dim = "<<n_dim<<" hasn't programmed yet..."<<endl;
      exit(1);
    }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //now we have the cell_type, we can find cell_bd_flag...
  //loop over all cells with cell_type == -1, check for neighbor with cell_type = 0
  cout<<" Search for boundary cells with cell_type = -1..."<<endl;

  for(int i=0; i<n_cell; i++)
    {
      if(cell_type[i] == -1)
	{
	  for(int j=0; j<cell_neig_ct[i]; j++)
	    {
	      if(cell_type[cell_neig[i][j]] == 0)
		{
		  cell_bd_flag[i] = 1;
		  break;
		}
	    }
	}
    }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //make the grid list for the boundary cells....
  //loop over cell_bd_flag for boundary cells...
  cout<<" generating bd_cell_index ..."<<endl;
  
  for(int i=0; i<n_cell; i++)
    {
      if(cell_bd_flag[i] == 1)
	{
	  //the list containing all healthy cells on the boundary
	  healthy_bd_list.push_back(i);

	  int temp_bd_grid_index;

	  if(n_dim == 2)
	    {
	      int xt = (int)floor(cell_coord[i][0]*n_grid/length);
	      int yt = (int)floor(cell_coord[i][1]*n_grid/length);
	      
	      temp_bd_grid_index = yt*n_grid + xt; 
	    }
	  else if(n_dim == 3)
	    {
	      int xt = (int)floor(cell_coord[i][0]*n_grid/length);
	      int yt = (int)floor(cell_coord[i][1]*n_grid/length);
	      int zt = (int)floor(cell_coord[i][2]*n_grid/length);
	      
	      temp_bd_grid_index = zt*n_grid*n_grid + yt*n_grid + xt; 
	    }

	  if(temp_bd_grid_index > n_tot_grid)
	    {
	      cout<<" the grid index is not correct...! recheck!"<<endl;
	      exit(1);
	    }

	  //found the grid_index, push back this cell in the grid list...
	  bd_grid_list[temp_bd_grid_index][bd_grid_ct[temp_bd_grid_index]] = i;
	  bd_grid_ct[temp_bd_grid_index]++;
	  
	}
    }
  
  //a list of grid which contains boundary cells...
  for(int i=0; i<n_tot_grid; i++)
    {
      if(bd_grid_ct[i]>0)
	{
	  bd_grid_index.push_back(i);
	 
	}
    }

}

//return the distance between the cells with indexI and indexJ
double dist_cell(int indexI, int indexJ)
{
  double dist = 0.0;

  for(int i=0; i<n_dim; i++)
    dist += (cell_coord[indexI][i]-cell_coord[indexJ][i])*(cell_coord[indexI][i]-cell_coord[indexJ][i]);

  return sqrt(dist);
}

//return the distance between the grids with indexI and indexJ
double dist_grid(int indexI, int indexJ)
{
  double dist = 0.0;

  for(int i=0; i<n_dim; i++)
    dist += (bd_grid_coord[indexI][i]-bd_grid_coord[indexJ][i])*(bd_grid_coord[indexI][i]-bd_grid_coord[indexJ][i]);

  return sqrt(dist);
}

//return the distance between a cell and the center of the tumor
double radial_dist(int index)
{
  double dist = 0.0;

  for(int i=0; i<n_dim; i++)
    dist += (cell_coord[index][i] - tumor_ct_coord[i])*(cell_coord[index][i] - tumor_ct_coord[i]);

  return sqrt(dist); 
}


//search the neighbors of edge cell index and return the closest healthy cell 
//also compute the dist Le
//if no healthy cell are found from the neighbors of the edge cell, return -1
int close_edge_cell(int edge_index, int cell_index, double& Le_dist)
{
  double mini_dist = 10000000.0;
  double dist;
  int new_edge_index = -1;

  //search the neighbors of the edge cell...
  for(int i=0; i<cell_neig_ct[edge_index]; i++)
    {
      int temp_index = cell_neig[edge_index][i];

      //make sure it's a healthy cell..
      //################# THIS can be modified to allow boundary nutrition feeding...
      if(cell_type[temp_index] == 0 || cell_type[temp_index] == -1)
	{
	  
	  dist = dist_cell(temp_index, cell_index);
	  
	  if(dist < mini_dist)
	    {
	      mini_dist = dist;
	      new_edge_index = temp_index;
	    }
	}
    }
  //we have to perform a higher order neighor search... see whether it works,
  //if still doesn't work, we need to deal with edge_list rigorously...
  if(new_edge_index != -1)
    {
      Le_dist = mini_dist;
  
      //cout<<"new_edge_index = "<<new_edge_index<<endl;
      return new_edge_index;
    }
  else
    {
      //we search for the second order neighbors...
      mini_dist = 10000000.0;
      
      for(int i=0; i<cell_neig_ct[edge_index]; i++)
	{
	  int temp_index1 = cell_neig[edge_index][i];

	  for(int j=0; j<cell_neig_ct[temp_index1]; j++)
	    {
	      int temp_index2 = cell_neig[temp_index1][j];

	      //make sure it's a healthy cell..
	      //################# THIS can be modified to allow boundary nutrition feeding...
	      if(cell_type[temp_index2] == 0 || cell_type[temp_index2] == -1)
		{
		  
		  dist = dist_cell(temp_index2, cell_index);
		  
		  if(dist < mini_dist)
		    {
		      mini_dist = dist;
		      new_edge_index = temp_index2;
		    }
		}
	    }
	}
      if(new_edge_index != -1)
	{
	  Le_dist = mini_dist;
	  
	  //cout<<"new_edge_index = "<<new_edge_index<<endl;
	  return new_edge_index;
	}
      else
	{
	  //this is reserved for a 3rd order neighbor search, if necessary...
	  mini_dist = 10000000.0;
      
	  for(int i=0; i<cell_neig_ct[edge_index]; i++)
	    {
	      int temp_index1 = cell_neig[edge_index][i];
	      
	      for(int j=0; j<cell_neig_ct[temp_index1]; j++)
		{
		  int temp_index2 = cell_neig[temp_index1][j];
		  
		  for(int k=0; k<cell_neig_ct[temp_index2]; k++)
		    {
		      int temp_index3 = cell_neig[temp_index2][k];

		      //make sure it's a healthy cell..
		      //################# THIS can be modified to allow boundary nutrition feeding...
		      if(cell_type[temp_index3] == 0 || cell_type[temp_index3] == -1)
			{
			  
			  dist = dist_cell(temp_index3, cell_index);
			  
			  if(dist < mini_dist)
			    {
			      mini_dist = dist;
			      new_edge_index = temp_index3;
			    }
			}
		    }
		  
		  
		}
	    }
	  
	  
	  //I won't go to 4th order neighbors...
	  Le_dist = mini_dist;
	  
	  //cout<<"new_edge_index = "<<new_edge_index<<endl;
	  return new_edge_index;
	}
    }
}

double inner_product(double* vect1, double* vect2)
{
  double sum = 0.0;
  double sum1 = 0.0;
  double sum2 = 0.0;

  for(int i=0; i<n_dim; i++)
    {
      sum += vect1[i]*vect2[i];
      sum1 += vect1[i]*vect1[i];
      sum2 += vect2[i]*vect2[i];
    }

  return sum/(sqrt(sum1)*sqrt(sum2));
}

//directly search all the boundary points...
int close_bd_cell2(int index, int edge_index)
{
  int bd_cell_index;

  double* vect1 = new double[n_dim];
  double* vect2 = new double[n_dim];


  //this local rule is important to determine tumor shape
  for(int i=0; i<n_dim; i++)
    {
      vect1[i] = cell_coord[edge_index][i]-tumor_ct_coord[i];
    }
  
  double max_ang = -10000000.0;
  double tmp_ang;

  for(int i=0; i<healthy_bd_list.size(); i++)
    {
      int temp_bd_index = healthy_bd_list[i];
      
      for(int k=0; k<n_dim; k++)
	{
	  vect2[k] = cell_coord[temp_bd_index][k] - tumor_ct_coord[k];
	}

      tmp_ang = inner_product(vect1, vect2);

      if(tmp_ang > max_ang)
	{
	  max_ang = tmp_ang;

	  bd_cell_index = temp_bd_index;
	}
    }
  
  return bd_cell_index;
}

//return the closest cells on the boundary...
//this is fast, but cannot correclty reflect the boundaries ...
/*
int close_bd_cell(int index, int edge_index)
{
  int bd_cell_index;

  double* vect1 = new double[n_dim];
  double* vect2 = new double[n_dim];
  //double* vect3 = new double[n_dim];
  

  //this local rule is important to determine tumor shape
  for(int i=0; i<n_dim; i++)
    {
      vect1[i] = cell_coord[edge_index][i]-tumor_ct_coord[i];
    }

  //cout<<"here1"<<endl;

  //using the box...
  if(n_dim == 2)
    {
      //int xt = (int)floor(cell_coord[index][0]*n_grid/length);
      //int yt = (int)floor(cell_coord[index][1]*n_grid/length);

      int xt = (int)floor(tumor_ct_coord[0]*n_grid/length);
      int yt = (int)floor(tumor_ct_coord[1]*n_grid/length);
      
      int temp_grid = yt*n_grid + xt;
      int bd_grid;

      //find the grid that closest to the current one...
      for(int i=0; i<bd_grid_index.size(); i++)
	{
	  for(int k=0; k<n_dim; k++)
	    vect2[k] = bd_grid_coord[bd_grid_index[i]][k] - bd_grid_coord[temp_grid][k];

	  double max_ang = -100000000.0;
	  double tmp_ang = inner_product(vect1, vect2);

	  //cout<<"temp_ang = "<<tmp_ang<<endl;

	  if(tmp_ang>max_ang && tmp_ang <(1- TOL))
	    {
	      max_ang = tmp_ang;
	      bd_grid = bd_grid_index[i];
	    }
	}

      //cout<<"here2"<<endl;

      //cout<<"bd_grid"<<bd_grid<<endl;

      //loop over the cells within bd_grid...
      for(int i=0; i<bd_grid_ct[bd_grid]; i++)
	{

	  //cout<<"here3"<<endl;
	  
	  int temp_cell_index = bd_grid_list[bd_grid][i];

	   
	  
	  //for(int k=0; k<n_dim; k++)
	    //vect2[k] = cell_coord[temp_cell_index][k] - cell_coord[index][k];
	  


	  for(int k=0; k<n_dim; k++)
	    vect2[k] = cell_coord[temp_cell_index][k] - tumor_ct_coord[k]   ;

	  double max_ang = -100000000.0;
	  double tmp_ang = inner_product(vect1, vect2);

	  if(tmp_ang>max_ang && tmp_ang <(1- TOL))
	    {
	      max_ang = tmp_ang;
	      bd_cell_index = temp_cell_index;
	      
	    }
	}

     
    }
  else if(n_dim == 3)
    {
      int xt = (int)floor(cell_coord[index][0]*n_grid/length);
      int yt = (int)floor(cell_coord[index][1]*n_grid/length);
      int zt = (int)floor(cell_coord[index][2]*n_grid/length);
      
      int temp_grid = zt*n_grid*n_grid + yt*n_grid + xt;
      int bd_grid;

      //find the grid that closest to the current one...
      for(int i=0; i<bd_grid_index.size(); i++)
	{
	  for(int k=0; k<n_dim; k++)
	    vect2[k] = bd_grid_coord[bd_grid_index[i]][k] - bd_grid_coord[temp_grid][k];

	  double max_ang = -100000000.0;
	  double tmp_ang = inner_product(vect1, vect2);

	  if(tmp_ang>max_ang && tmp_ang <(1- TOL))
	    {
	      max_ang = tmp_ang;
	      bd_grid = bd_grid_index[i];
	    }
	}

      //loop over the cells within bd_grid...
      for(int i=0; i<bd_grid_ct[bd_grid]; i++)
	{
	  int temp_cell_index = bd_grid_list[bd_grid][i];

	  for(int k=0; k<n_dim; k++)
	    vect2[k] = cell_coord[temp_cell_index][k] - cell_coord[index][k];

	  double max_ang = -100000000.0;
	  double tmp_ang = inner_product(vect1, vect2);

	  if(tmp_ang>max_ang && tmp_ang <(1- TOL))
	    {
	      max_ang = tmp_ang;
	      bd_cell_index = temp_cell_index;
	      
	    }
	}
    }

  return bd_cell_index;
}
*/

//provide the coords of the proliferative cells
//initialize all cell type, all the list, number and vol for each type
//initialize the tumor edge
void Init_Tumor()
{
  //search all cells to find the one closest to the specifed values...
  cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
  cout<<"Initializing the tumor center..."<<endl;
  double temp_mini_dist2 = 1000000.0;
  double temp_dist2;
  int ct_index;

  for(int i=0; i<n_cell; i++)
    {
      temp_dist2 = 0;

      for(int j=0; j<n_dim; j++)
	temp_dist2 += (cell_coord[i][j] - tumor_ct_coord[j])*(cell_coord[i][j] - tumor_ct_coord[j]);

      if(temp_dist2 < temp_mini_dist2)
	{
	  temp_mini_dist2 = temp_dist2;
	  ct_index = i;
	}
    }

  //now we find the real center, give new values to tumor_ct_coord
  for(int i=0; i<n_dim; i++)
    tumor_ct_coord[i] = cell_coord[ct_index][i];

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //now we initialize all cell types, including both the list and the number...
  cout<<"Center cell index = "<<ct_index<<endl;
  cout<<"Initializing the prolif, qui and necrotic cells..."<<endl;

  cell_type[ct_index] = 1;

  n_pro = 1; n_tumor = 1;
  n_nec = n_qui = 0;

  prolif_list.push_back(ct_index);
  cell_close_edge_index[ct_index] = close_edge_cell(ct_index, ct_index, temp_dist2);

  vol_nec = vol_qui = 0.0;
  vol_pro = cell_vol[ct_index];

  //############### Can not think of anything else...
}


//for each proliferative cell, check to see whether it will divide or not
//- if divide, add new cells using push_back
//- if not, turn the cell to quiscent, remove the cell from prolif_list...
//the check for nearest edge cell is based on local rules...
void Update_Proliferative()
{

  int n_pro_old = n_pro;

  //loop over the list containing proliferative cells...
  for(int i=0; i<n_pro_old; i++)
    {
      double Le, Lt, R, Rmaxt; //the important distances
      double p_d, p_temp; //the dividing probability
      int edge_index; //index of the closest healthy cell on tumor edge
      int bd_index; //index of the closest boundary cells
      
      int cell_index = prolif_list[i];
      //cout<<" ***** Working with cell = "<<cell_index<<endl;

      if(cell_type[cell_close_edge_index[cell_index]] == 0)
	{
	  edge_index = cell_close_edge_index[cell_index];
	  Le = dist_cell(cell_index, edge_index); //this computation is duplicate...
	  //cout<<"edge_index0 = "<<edge_index<<endl;
	}
      else
	{
	  edge_index = close_edge_cell(cell_close_edge_index[cell_index], cell_index, Le);
	  //cout<<"edge_index1 = "<<edge_index<<endl;
	}
     
      
      if(edge_index != -1)
	{
	  //update the closest edge cell associated with the prolif cell
	  cell_close_edge_index[cell_index] = edge_index;

	  //the cell might divide
	  //check whether within the distance using Lt, Le 
	  Lt = radial_dist(edge_index); //the distance of the edge cell to the tumor center
	  
	  //using multiplication instead of pow, which is more costy
	  double Delta_p = 1;
	  double B_d = 1;
	  double L_d1 = 1;
	  for(int j=0; j<(n_dim-1); j++)
	    {
	      Delta_p = Delta_p*Le;
	      B_d = B_d*b;
	      L_d1 = L_d1*Lt;
	    }
	  Delta_p = Delta_p*Le;
	  B_d = B_d*b;
	  L_d1 = L_d1*B_d;
	  
	  //cout<<"Delta_p = "<<Delta_p<<"   L_d1 = "<<L_d1<<endl; 
	  
	  if(Delta_p<=L_d1)
	    {
	      //- if so, check whether divide using p0, using R, Rmax
	      bd_index = close_bd_cell2(cell_index, edge_index);
	      
	      

	      R = radial_dist(cell_index);
	      Rmaxt = radial_dist(bd_index);

	      //cout<<"here2"<<endl;
	      //double ratio = R/Rmaxt;
	      //if(ratio>=1) ratio = 1;

	      p_d = p0*(1-R/Rmaxt); //probability of division
	      p_temp = (double)(rand()%RAND_MAX)/(double)RAND_MAX;

	      //cout<<"here3"<<endl;

	      //the cell will divide, otherwise nothing is changed...
	      if(p_temp < p_d && cell_type[edge_index] != -1)
		{
		  //cout<<"A new prolif cell is found ..."<<endl;

		  cell_type[edge_index] = 1;

		  n_pro++;
		  //update the n_pro and n_pro_old, the new prolif cell is treated next round...
		  //this makes sure we don't  go into a dead loop 

		  prolif_list.push_back(edge_index);

		  vol_pro += cell_vol[edge_index];

		  //now we need to assing a neighbor to this new cell...
		  double temp_d;
		  int tmp_neig_new_cell = close_edge_cell(edge_index, edge_index, temp_d);
		  
		  if(tmp_neig_new_cell != -1)
		    {
		      cell_close_edge_index[edge_index] = tmp_neig_new_cell;
		    }
		  else
		    {
		      //it means this cell is now on the tumor boundary...
		      cell_close_edge_index[edge_index] = edge_index;
		    }
		    

		  //only proliferative cell contribute to the total number of tumor cells
		  //a new cell is added, need to update the tumor center...
		  for(int k=0; k<n_dim; k++)
		    {
		      tumor_ct_coord[k] = n_tumor*tumor_ct_coord[k] + cell_coord[edge_index][k];
		      tumor_ct_coord[k] = tumor_ct_coord[k]/(double)(n_tumor + 1);
		    }

		  //cout<<"ct = "<<tumor_ct_coord[0]<<" "<<tumor_ct_coord[1]<<endl;
		    
		  n_tumor++;

		  //cout<<"The new prolif cell is added ..."<<endl;
		}
	      /*
	      else
		{
		  cout<<"The cell is not dividing now..."<<endl;
		}
	      */
	      //if not dividing, nothing needs to be done...
	    }
	  else
	    { //- if not, turn into quiscent cell....
	      //cout<<"turning the into quiscent 1..."<<endl;

	      cell_type[cell_index] = 2; 
	      n_qui ++;
	      quiscent_list.push_back(cell_index);
	      vol_qui += cell_vol[cell_index];

	      //erase the prolif cell from the prolif list...
	      n_pro --;
	      n_pro_old --; //this is needed in both case...
	      
	      vol_pro = vol_pro - cell_vol[cell_index];

	      prolif_list.erase(prolif_list.begin()+i);

	      //################## just put a checker here...
	      if(n_pro != prolif_list.size())
		{
		  cout<<"Counting is not correct for n_pro and prolif_list.size()! Abort! "<<endl;
		  exit(1);
		}
	    }
	 
	 
	  //update the list, num., and vol.,
	}
      else
	{
	  //################## There is a tricky thing about the closest distance... 
	  //################## We could go second round check 
	  //no healthy cells are found, the pro is turned into quiscent, 
	  //update both list, vol. num.,

	  //cout<<"turning the into quiscent 2..."<<endl;

	  cell_type[cell_index] = 2; 
	  n_qui ++;
	  quiscent_list.push_back(cell_index);
	  vol_qui += cell_vol[cell_index];
	  
	  //erase the prolif cell from the prolif list...
	  n_pro --;
	  n_pro_old --; //this is needed in both case...
	  
	  vol_pro = vol_pro - cell_vol[cell_index];
	  
	  prolif_list.erase(prolif_list.begin()+i);
	  
	  //################# just put a checker here...
	  if(n_pro != prolif_list.size())
	    {
	      cout<<"Counting is not correct for n_pro and prolif_list.size()! Abort! "<<endl;
	      exit(1);
	    }
	  
	}
      
    }

  //cout<<"return to the main iteration..."<<endl;
}

//our current way of searching closest tumor edges requires updateing Quiscent cells first
void Update_Quiscent()
{
  //loop over all quiscent cells, 
  //check whether they will turn into necrotic cells...
  for(int i=0; i<n_qui; i++)
    {
      int cell_index = quiscent_list[i];
      int edge_index; //index of the closest healthy cell on tumor edge

      double Lt, Le;
     
      if(cell_type[cell_close_edge_index[cell_index]] == 0)
	{
	  edge_index = cell_close_edge_index[cell_index];
	  Le = dist_cell(cell_index, edge_index); //this computation is duplicate...
	}
      else
	{
	  edge_index = close_edge_cell(cell_close_edge_index[cell_index], cell_index, Le);
	}

      
      if(edge_index != -1)
	{
	  cell_close_edge_index[cell_index] = edge_index; //update the closest tumor edge cell...

	  //check if the distance is OK for quiscent to live
	  //otherwise turn into necrotic
	  Lt = radial_dist(edge_index);

	  //using multiplication instead of pow, which is more costy
	  double Delta_p = 1;
	  double A_d = 1;
	  double L_d1 = 1;
	  for(int j=0; j<(n_dim-1); j++)
	    {
	      Delta_p = Delta_p*Le;
	      A_d = A_d*a;
	      L_d1 = L_d1*Lt;
	    }
	  Delta_p = Delta_p*Le;
	  A_d = A_d*a;
	  L_d1 = L_d1*A_d;


	  if(Delta_p > L_d1)
	    { // turn into necrotic cell....

	      //cout<<"turning the into nectrotic 1..."<<endl;
	      
	      cell_type[cell_index] = 3; 
	      n_nec++;
	      necrotic_list.push_back(cell_index);
	      vol_nec += cell_vol[cell_index];

	      //erase the quiscent cell from the quiscent list...
	      n_qui --;
	      
	      vol_qui = vol_qui - cell_vol[cell_index];

	      quiscent_list.erase(quiscent_list.begin()+i);

	      //################## just put a checker here...
	      if(n_qui != quiscent_list.size())
		{
		  cout<<"Counting is not correct for n_pro and prolif_list.size()! Abort! "<<endl;
		  exit(1);
		}
	    }
	  
	}
      else
	{
	  //################## again, a second around search might needed here...
	  //no healthy cells close enough, turn into necrotic 


	  //cout<<"turning the into nectrotic 2..."<<endl;

	  cell_type[cell_index] = 3; 
	  n_nec++;
	  necrotic_list.push_back(cell_index);
	  vol_nec += cell_vol[cell_index];
	  
	  //erase the quiscent cell from the quiscent list...
	  n_qui --;
	  
	  vol_qui = vol_qui - cell_vol[cell_index];
	  
	  quiscent_list.erase(quiscent_list.begin()+i);
	  
	  //################## just put a checker here...
	  if(n_qui != quiscent_list.size())
	    {
	      cout<<"Counting is not correct for n_qui and quiscent_list.size()! Abort! "<<endl;
	      exit(1);
	    }
	}
      
    }
}


//Get the tumor radius...
//this is a complete update each time the function is called...
//search for all proliferative cells, find those on the boundary...
double tumor_radius()
{
  int edge_ct = 0; //the total number 

  double edge_dist = 0.0; //the distance...

  for(int i=0; i<n_pro; i++)
    {
      int pro_index = prolif_list[i];

      for(int j=0; j<cell_neig_ct[pro_index]; j++)
	{
	  int neig_index = cell_neig[pro_index][j];
	  
	  if(cell_type[neig_index] == 0)
	    {
	      double tmp_dist = 0.0;
	      
	      for(int k=0; k<n_dim; k++)
		tmp_dist += (cell_coord[pro_index][k]-tumor_ct_coord[k])*(cell_coord[pro_index][k]-tumor_ct_coord[k]);

	      edge_dist += sqrt(tmp_dist);
	      edge_ct ++;

	      continue;
	    }
	}
    }

  return edge_dist/(double)edge_ct;
}

//return the volume of the tumor...
double tumor_vol()
{
  return (vol_pro+vol_qui+vol_nec);
}

//we print out the tumor at specified time point...
//this is not very concise in terms of codeing ....
void Print_Tumor(int t)
{
  if(t==t1)
    {
      ofstream outfile0;
      outfile0.open("Tumor_t1.dat");
      
      outfile0<<"# t1 = "<<t1<<endl;
      outfile0 <<"# n_pro = "<<n_pro<<"   vol_pro = "<<vol_pro<<endl;
      outfile0 <<"# n_qui = "<<n_qui<<"   vol_qui = "<<vol_qui<<endl;
      outfile0 <<"# n_nec = "<<n_nec<<"   vol_nec = "<<vol_nec<<endl;
      

      int index;

      for(int i=0; i<n_pro; i++)
	{
	  index = prolif_list[i];
	  
	  for(int j=0; j<n_dim; j++)
	    outfile0 <<cell_coord[index][j] << "  ";
	  outfile0 <<endl;
	}	 
      outfile0<<"################################"<<endl<<endl;

      for(int i=0; i<n_qui; i++)
	{
	  index = quiscent_list[i];
	  
	  for(int j=0; j<n_dim; j++)
	    outfile0 <<cell_coord[index][j] << "  ";
	  outfile0 <<endl;
	}	 
      outfile0<<"################################"<<endl<<endl;

      for(int i=0; i<n_nec; i++)
	{
	  index = necrotic_list[i];
	  
	  for(int j=0; j<n_dim; j++)
	    outfile0 <<cell_coord[index][j] << "  ";
	  outfile0 <<endl;
	}	 
      outfile0<<"################################"<<endl<<endl;
      
      outfile0.close();
    }
  if(t==t2)
    {
      ofstream outfile0;
      outfile0.open("Tumor_t2.dat");
      
      outfile0<<"# t2 = "<<t2<<endl;
      outfile0 <<"# n_pro = "<<n_pro<<"   vol_pro = "<<vol_pro<<endl;
      outfile0 <<"# n_qui = "<<n_qui<<"   vol_qui = "<<vol_qui<<endl;
      outfile0 <<"# n_nec = "<<n_nec<<"   vol_nec = "<<vol_nec<<endl;
      

      int index;

      for(int i=0; i<n_pro; i++)
	{
	  index = prolif_list[i];
	  
	  for(int j=0; j<n_dim; j++)
	    outfile0 <<cell_coord[index][j] << "  ";
	  outfile0 <<endl;
	}	 
      outfile0<<"################################"<<endl<<endl;

      for(int i=0; i<n_qui; i++)
	{
	  index = quiscent_list[i];
	  
	  for(int j=0; j<n_dim; j++)
	    outfile0 <<cell_coord[index][j] << "  ";
	  outfile0 <<endl;
	}	 
      outfile0<<"################################"<<endl<<endl;

      for(int i=0; i<n_nec; i++)
	{
	  index = necrotic_list[i];
	  
	  for(int j=0; j<n_dim; j++)
	    outfile0 <<cell_coord[index][j] << "  ";
	  outfile0 <<endl;
	}	 
      outfile0<<"################################"<<endl<<endl;

      outfile0.close();
      
    }
  if(t==t3)
    {
      ofstream outfile0;
      outfile0.open("Tumor_t3.dat");
      
      outfile0<<"# t3 = "<<t3<<endl;
      outfile0 <<"# n_pro = "<<n_pro<<"   vol_pro = "<<vol_pro<<endl;
      outfile0 <<"# n_qui = "<<n_qui<<"   vol_qui = "<<vol_qui<<endl;
      outfile0 <<"# n_nec = "<<n_nec<<"   vol_nec = "<<vol_nec<<endl;
      

      int index;

      for(int i=0; i<n_pro; i++)
	{
	  index = prolif_list[i];
	  
	  for(int j=0; j<n_dim; j++)
	    outfile0 <<cell_coord[index][j] << "  ";
	  outfile0 <<endl;
	}	 
      outfile0<<"################################"<<endl<<endl;

      for(int i=0; i<n_qui; i++)
	{
	  index = quiscent_list[i];
	  
	  for(int j=0; j<n_dim; j++)
	    outfile0 <<cell_coord[index][j] << "  ";
	  outfile0 <<endl;
	}	 
      outfile0<<"################################"<<endl<<endl;

      for(int i=0; i<n_nec; i++)
	{
	  index = necrotic_list[i];
	  
	  for(int j=0; j<n_dim; j++)
	    outfile0 <<cell_coord[index][j] << "  ";
	  outfile0 <<endl;
	}	 
      outfile0<<"################################"<<endl<<endl;
      
    }
  if(t==t4)
    {
      ofstream outfile0;
      outfile0.open("Tumor_t4.dat");
      
      outfile0<<"# t4 = "<<t4<<endl;
      outfile0 <<"# n_pro = "<<n_pro<<"   vol_pro = "<<vol_pro<<endl;
      outfile0 <<"# n_qui = "<<n_qui<<"   vol_qui = "<<vol_qui<<endl;
      outfile0 <<"# n_nec = "<<n_nec<<"   vol_nec = "<<vol_nec<<endl;
      

      int index;

      for(int i=0; i<n_pro; i++)
	{
	  index = prolif_list[i];
	  
	  for(int j=0; j<n_dim; j++)
	    outfile0 <<cell_coord[index][j] << "  ";
	  outfile0 <<endl;
	}	 
      outfile0<<"################################"<<endl<<endl;

      for(int i=0; i<n_qui; i++)
	{
	  index = quiscent_list[i];
	  
	  for(int j=0; j<n_dim; j++)
	    outfile0 <<cell_coord[index][j] << "  ";
	  outfile0 <<endl;
	}	 
      outfile0<<"################################"<<endl<<endl;

      for(int i=0; i<n_nec; i++)
	{
	  index = necrotic_list[i];
	  
	  for(int j=0; j<n_dim; j++)
	    outfile0 <<cell_coord[index][j] << "  ";
	  outfile0 <<endl;
	}	 
      outfile0<<"################################"<<endl<<endl;
      
    }
  if(t==t5)
    {
      ofstream outfile0;
      outfile0.open("Tumor_t5.dat");
      
      outfile0<<"# t5 = "<<t5<<endl;
      outfile0 <<"# n_pro = "<<n_pro<<"   vol_pro = "<<vol_pro<<endl;
      outfile0 <<"# n_qui = "<<n_qui<<"   vol_qui = "<<vol_qui<<endl;
      outfile0 <<"# n_nec = "<<n_nec<<"   vol_nec = "<<vol_nec<<endl;
      

      int index;

      for(int i=0; i<n_pro; i++)
	{
	  index = prolif_list[i];
	  
	  for(int j=0; j<n_dim; j++)
	    outfile0 <<cell_coord[index][j] << "  ";
	  outfile0 <<endl;
	}	 
      outfile0<<"################################"<<endl<<endl;

      for(int i=0; i<n_qui; i++)
	{
	  index = quiscent_list[i];
	  
	  for(int j=0; j<n_dim; j++)
	    outfile0 <<cell_coord[index][j] << "  ";
	  outfile0 <<endl;
	}	 
      outfile0<<"################################"<<endl<<endl;

      for(int i=0; i<n_nec; i++)
	{
	  index = necrotic_list[i];
	  
	  for(int j=0; j<n_dim; j++)
	    outfile0 <<cell_coord[index][j] << "  ";
	  outfile0 <<endl;
	}	 
      outfile0<<"################################"<<endl<<endl;
      
    }
}


//rt...
void Print_Boundary()
{
  //check whether the boundary cells are correct....
  ofstream fout;
  fout.open("boundary.txt");

  int temp_size = healthy_bd_list.size();
  for(int i=0; i<temp_size; i++)
    {
      for(int j=0; j<n_dim; j++)
	fout<<cell_coord[healthy_bd_list[i]][j]<<"   ";
      fout<<endl;
    }

  fout.close();

  /*
  //check whether the grids work
  fout.open("bd_grid.txt");
  for(int i=0; i<bd_grid_index.size(); i++)
    {
      int temp_grid_index = bd_grid_index[i];

      //for the bd cells in each grid...
      for(int j=0; j<bd_grid_ct[temp_grid_index]; j++)
	{
	  int temp_cell_index = bd_grid_list[temp_grid_index][j];

	  for(int k=0; k<n_dim; k++)
	    fout<<cell_coord[temp_cell_index][k]<<"   ";
	  fout<<endl;
	}

      fout<<"###########"<<endl<<endl;
    }
  fout.close();
  */
}


main()
{

  srand(time(NULL));
  
  Init_Data();

  Init_Boundary(bd_type);
  Print_Boundary();

  //~~~~~~~~~~~~~~~~~~~~~~~~~
  Init_Tumor();

  ofstream outRad;
  outRad.open("tumor_rad.txt");
  
  ofstream outVol;
  outVol.open("tumor_vol.txt");

  for(int t=0; t<iterations; t++)
    {
      cout<<"At T = "<<t<<"..."<<endl;
      
      Update_Quiscent();

      //cout<<"here1"<<endl;
      
      Update_Proliferative();

      //cout<<"here2"<<endl;

      Print_Tumor(t);

      cout<<"n_tumor = "<<n_tumor<<endl;
      cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

      outRad << t <<"    "<<tumor_radius()<<endl;
      outVol << t <<"    "<<tumor_vol()<<endl;
      
    }


  outRad.close();
  outVol.close();
 
}
