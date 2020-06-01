//A celluar automaton model to simulate tumor growth in d-dimensions
//Using a sequential implementation, with sophiscated data structure to maximize efficiency
//Complementary to the paralell version by Jana


//##################################################
//This version uses the latest CA model to study the mechanism of dormancy caused by immunosurveillance

//author: Originall written by Yang Jiao, modified by Duyu Chen, duyu@princeton.edu
//started: Jan. 22, 2013



/* Input file required for algorithm: input.dat							*/
/*  - This file has parameters for algorithm and input files 					*/
/*      and specify the name of the other three input files and output file name                */
/*  - Input file 1: inputPoints: takes list of cell centers sorted from furthest from center 	*/
/*      to closest										*/
/*  - Input file 2: inputArea: area of each automaton cell					*/
/*  - Input file 3: inputNeighbors: neighbors of each automaton cell				*/



//Major Features of this Invasion CA Model:
//(1) (We use a grain-coarsening idea to check cells on the boundary, i.e., first divide the whole 
//    domain into small grids, then only make a list of grid that contains boundary cells;
//    when computing distance between bd cells and any particular cell, first compute the 
//    distance between the two grids, when the closest grid is found, then compute the bd cells.)
//
//(2) For the tumor edge, it is dynamic, using grid requries frequent update, not efficient
//    Instead, we use a local-rule, i.e., each proliferative has a closest-edge cell, if this 
//    cell becomes tumor cell, all its neighbors are searched for new healthy cell or bd cell
//    The neighbor search can be 2nd order or even higher order in depth, meaning searching 
//    the neighbors of a neighbor cell. But for now, we just search the first neighbors.
//  !!!!!- Note: this neighbor-search-method requires update quiscent cells first...
//  !!!!!- Note: to take into account the effects of collagen density, the search return 
//         a cell_index that is mimium of dist+cell_density
//
//(3) For invasive cells,
//   (a) at each interaction (stage), when updating proliferation cells, a fraction of new 
//       cells are selected to be possible invasive, with a intrisic degradation ability, i.e., rho_degrad
//       and cell motility, i.e., the largest steps it can move in each update
//       If rho_degrad is smaller than the smallest rho of its neighbors, this cell is turned back into proliferative
//       Else, it is inserted into the invasive cell list, and make the movement
//   (b) After updating the non-invasive cells, update the invasive list with only the 
//       old invasive cells (the new cells are updated when they are generated)
//       including making the movement, but not proliferating 
//
//(4) Modeling both long-range and short-range elastic interactions
//    (a) long-range: through the shape of boundary, by modifying the probability of proliferation
//    (b) short-range: through the local ECM density, by modifying the probability of proliferation
//        i.e.. p_prof_ECM = 1 - rho_ECM
//
//(5a) We print-out the tumor configuration every 5 itetrations to make movies of the growth



//Novel features of the Pressure CA Model
//(1) In the previous model, it is considered ECM is completely degraded before a tumor 
//    cell can take up the ECM. Here we consider only a portion \chi of ECM are 
//    degraded and the rest are pushed away, which in turn can exert a pressure on the tumor
//
//(2) The density of ECM is considered to be a function of time, since ECM is deformed.
//    In particular, rho_ECM depends on the current tumor radius, number of tumor cells
//    and number of ECM cells. Please see the associated notes for the exact math formula.
//
//(3) The probability of division is purely a function of rho_ECM now, i.e., 
//    pt = max{0, p0*(1-rho_ECM)}
//
//(4) We may also consider the daughter cells take up the automaton cell 
//    further into ECM will get growth advantage, due stress concentration,
//    which is proportional to its scaled local successive distance (see the implementation for details)
//
//(5) Also need to update the rules for invasive cells...
//    In the worest case, can still program pressure induced invasion
//
//(6) Add function to compute specific surface and asphericity


//
//(7) considering environment triggered invasion by explicitly considering 
//    the effects of cell-cell adhesion: this is implemented as a pre-check for 
//    invasive cells: only those with a number of contact neighbors smaller 
//    than a prescribed value.

//version: 03/26/2013
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
#include <sstream>

#define NAME_LEN 256 //for the name of the input and output files....
//#define RAND_MAX 10000000 //denominator for maximum number
#define TOL 0.00000001 //a numerical tolerance...

int n_dim; //the spatial dimension
//int n_grid; //the linear size of grid along each direction
int n_neig; //this is not for read in, computed as 3^n_dim, not a necessary num.
double vol_growth_domain = 3.14159265/4.0;

double rho0 = 0.30; //this specifies the intitial density for homogeneous and random ECM

double wp = 2*rho0;//parameter measure division prob. reduced by pressure
//Tumor algorithm variables
double p_t = 0.005; //transition rate
double frac_gentle = 0.90;
double frac_mal = 0.10;
double frac_kill = 0.20;
double r_gentle = 0.75;
double r_mal = 0.15;
double a, b, p0, Rmax; //the four minimalist parameters  
//Rmax is the characteristic para. used for mixing the effect of maximizing oxygen, nutrient supply and degradation degree
double p_inva; //the fraction of invasive cells on the tumor edge over total number of cells on the tumor edge

double degradation_max; //the maximum possible level of degradation for invasive cells
double chi; //the average degradation fraction of the ECM
            //(in the current implementation, we consider chi = 0.2*degradation_max)
            //i.e., the degradation of ECM by noninvasive cell is small compared to that of invasive cells

int motility_max; //the maximum possible trials for invasive cell, but each one should have its own motiltiy
double length; //the length of a square domain 
int iterations, t1; //total number of iterations and time points for output
int bd_type; //specify the boundary shape and type (whether can provide nuritiion or not)
int density_type; //specify the distribution of collagen density
                  // 1 - random distribution

//for the ECM
double den_ave_ECM; //the scaled time-dependent average ECM density 
double den_ave0; //the initial average density ...
double sum_den = 0.0;// the sum of all ECm density, for computing den_ave_ECM at each day
double sum_den_degrad = 0.0; //the sum of degraded ECM density, for computing den_ave_ECM at each day
int adhesion = 0; //the minimum number of neighbors an invasive has to keep attached to the central tumor
                    //an invasive cell will never invade if adhesion = 0.


//names of input files
char dataPoints[NAME_LEN];	
char dataVol[NAME_LEN];	
char dataNeighbors[NAME_LEN];
//name of the output file for tumor radius and volume
//char outRadVol[NAME_LEN];

//information for cells
vector< vector<double> > cell_coord; //the coordinates of the cells
vector<int> cell_type; //-1 healthy outside bd; 0 - healthy within bd; 1 - proliferative; 2 - quiscent; 3 - necrotic; 100 - invasive
vector<int> cell_status;//0 dormany; 1 gentle; 2 malignant
vector<double> cell_density; //in [0, 1] the mass density of ECM within that CA cell
vector<int> cell_bd_flag; //whehter this cell is on the boundary
vector<int> cell_edge_flag; //whether this cell is on the tumor edge
vector<int> cell_close_edge_index; //the index of the closest tumor edge cells to the current cell

vector<double> cell_vol; //volume of each cell
vector< vector<int> > cell_neig; //the neighbor list of each cell
vector<int> cell_neig_ct; //the counter for the neighbors

int n_cell; //the total number of cells (including those within and outside the growth region)
double cell_len; //a characteristic length of cells,not used in the current model


//informtion of for the tumor
vector<double> tumor_ct_coord; //the coordinates of the tumor center

vector<int> necrotic_list; //the list of cells associated with tumors, constantly updated
vector<int> quiscent_list;
vector<int> prolif_list;
vector<int> invasive_list; //the list of invasive cells, 
                           //when updating the list, don't delete entry, just change the cell index to take into account the movement
vector<int> inva_motility; //the motility of each invasive cells
vector<double> inva_degradation; //the degradation ability of each invasive cell

int n_nec, n_qui, n_pro, n_pro_acti, n_inva, n_tumor; //numbers
double vol_nec, vol_qui, vol_pro, vol_pro_acti, vol_inva; //volumes

double Rad_T; //tumor radius
double Rad_Tmax; //max temp tumor radius for asphericty
double Rad_Tmin; //min temp tumor radius for asphericty

//asphericy
double alpha_T; 

//specific surface
double ss_T; 

//vector<int> healthy_edge_list; //the list of healthy cells at tumor edge, constantly updated
vector<int> healthy_bd_list; 
//the list of healthy cells at the boundary, static!

//information of the grids, this is for the quick search of boundary cells...
//int  n_tot_grid;

//vector< vector<int> > bd_grid_list; //a list of grid in which each one has a list of index of boundary cells
//vector<int> bd_grid_index; //the index of grid that contain boundary cells
//vector< vector<double> > bd_grid_coord; //the coordinates of the grid center...
//vector<int> bd_grid_ct; //the counter
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
  //infile >> n_grid; //num.of grid, not in the original input file
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
  infile >> p_inva; //not in the original input...
  infile >> degradation_max; //the max is 1, should always be less than this
  infile >> chi; //the average degradation ability of the nonivasive cells
  infile >> motility_max; //we choose 5, meaning 5 invasive trails at max
  infile >> length;  infile >> iterations;
  infile >> t1; //infile >> t2; infile >> t3;  infile >> t4; infile >> t5; 
  infile >> bd_type; infile >> density_type;
  infile >> adhesion;
  infile >> dataPoints; infile >> dataVol;
  infile >> dataNeighbors; 
  //infile >> outRadVol;
  infile.close();

  //double check for correct input...
  cout<< ". n_dim = "<< n_dim <<endl;
  //cout<< ". n_grd = "<<n_grid <<endl;
  //cout<<" x0, y0, z0"<<endl; //not output, but do read in: the position of initial tumor cell
  cout << ". a = " << a << endl;
  cout << ". b = " << b << endl;
  cout << ". p0 = " << p0 << endl;
  cout << ". Rmax = " << Rmax << endl;
  cout << ". p_inva = " << p_inva << endl;
  cout << ". degradation_max = " << degradation_max << endl;
  cout << ". chi = " << chi << endl;  //the average degradation ability of noninvasive cells
  cout << ". motility_max = " << motility_max << endl;
  cout << ". length = " << length << endl; //length of square simulation domain, containing the boundary
  cout << ". iterations = " << iterations << endl;
  cout << ". t1 = " << t1 << endl; //time interval between successive configurations printing
  
  cout << ". bd_type = " << bd_type << endl; //specify the type of the boundary
  cout << ". density_type = " << density_type << endl; //specify the distribution of collagen density
  cout << ". adhesion = " << adhesion << endl;
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
	  cell_status.push_back(0);
      cell_bd_flag.push_back(0); //on bd, 1, else 0
      cell_edge_flag.push_back(0);
      cell_density.push_back(0.0); //the density of collagen 
      
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
  /*n_tot_grid = (int)pow(n_grid, n_dim);
  int n_cell_per_grid = (int)floor((double)2*n_cell/(double)n_tot_grid);
  */
  /*for(int i=0; i<n_tot_grid; i++)
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
*/ 
}


//now we initialize the collagen density
void Init_CellDensity(int density_type)
{
  if(density_type == 0) //uniform...
    {
      srand(time(NULL));

      double rho = rho0;

      for(int i=0; i<n_cell; i++)
	{ 
	  cell_density[i] = rho;
	}
      
      sum_den = n_cell*rho;
      
    }
  else if(density_type == 1)
    {
      srand(time(NULL));

      double rho;

      sum_den = 0.0;

      for(int i=0; i<n_cell; i++)
	{
	  rho = (double)(rand()%RAND_MAX)/(double)RAND_MAX;
	  
	  cell_density[i] = rho0*rho;

	  sum_den += rho0*rho;
	}
      
    }
  else if(density_type == 2)
    {
      //for the sine density profile
      double lambda = length/6.0;
      double rho;
      double pi = 3.1415926;


      sum_den = 0.0;
     
      //the following is dimensional indepedent
      for(int i=0; i<n_cell; i++)
	{
	  rho = 1.0;
	  
	  for(int j=0; j<n_dim; j++)
	    rho = rho*(sin(2*pi*cell_coord[i][j]/lambda)+1)*0.5;

	  //make sure the density is in [0, 1]

	  cell_density[i] = rho;

	  sum_den += rho;
	}
      
    }
  else if(density_type == 3)
    {
      //this is reserved for the Guassian model
      cout<<"Gaussian model is not programmed yet!"<<endl;
      exit(1);
    }
  else
    {
      cout<<"collagen density distribution for density_type = "<<density_type<<" hasn't programmed yet...!"<<endl;
      exit(1);
    }
}



//get the shape and the growth permitting domain, 
//get the boundary 
//(also the grid containing the list of boundary cells)
void Init_Boundary(int type)
{
  //first get the boundary~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cout<<" Generating the growth permitting domain..."<<endl;

  if(n_dim == 2)
    {
      if(type==0)
	{ //complex environment is an circle: (x-x0)^2/r^2 + (y-y0)^2/r^2 = 1
	  double x0 = 0.5*length;
	  double y0 = 0.5*length;
	  double ra = 0.5*length;
	  double rb = 0.5*length;
	  for(int i = 0; i<n_cell; i++)
	    {
	      double xDist = (cell_coord[i][0]-x0)/ra;
	      double yDist = (cell_coord[i][1]-y0)/rb;
	      double curDist = (xDist*xDist)+(yDist*yDist);
	      if(curDist>=1)
		cell_type[i] = -1;
	    }		
	}
      else if(type==1)
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
      if(type==0)
	{ //complex environment is an ellipsoid: (x-x0)^2/ra^2 + (y-y0)^2/rb^2 + (z-z0)^2/rc^2= 1
	  double x0 = 0.5*length;
	  double y0 = 0.5*length;
	  double z0 = 0.5*length;
	  double ra = 0.48;
	  double rb = 0.48;
	  double rc = 0.48;
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
		  healthy_bd_list.push_back(i);
		  break;
		}
	    }
	}
    }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //make the grid list for the boundary cells....
  //loop over cell_bd_flag for boundary cells...
  /*
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
*/
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
/*double dist_grid(int indexI, int indexJ)
{
  double dist = 0.0;

  for(int i=0; i<n_dim; i++)
    dist += (bd_grid_coord[indexI][i]-bd_grid_coord[indexJ][i])*(bd_grid_coord[indexI][i]-bd_grid_coord[indexJ][i]);

  return sqrt(dist);
}*/

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
		      
		      for(int h=0; h<cell_neig_ct[temp_index3]; h++)
			{

			  int temp_index4 = cell_neig[temp_index3][h];
			  //make sure it's a healthy cell..
			  //################# THIS can be modified to allow boundary nutrition feeding...
			  if(cell_type[temp_index4] == 0 || cell_type[temp_index4] == -1)
			  {
			    
			    dist = dist_cell(temp_index4, cell_index);
			    
			    
			    if(dist < mini_dist)
			      {
				mini_dist = dist;
				new_edge_index = temp_index4;
			      }
			  }
			}
		    }
		  
		  
		}
	    }
	  
	  
	  //I won't go to 5th order neighbors...
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

//directly search all the boundary points...which is closet along the center-edge line
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


//get and update average ECM density
double Get_den_ave_ECM()
{
  double vol_T = vol_nec + vol_qui + vol_pro;//total volume of noninvasive cells

  double vol_D = length*length;//volume of the square domain
  
  double temp_rho = (sum_den*vol_D/n_cell-sum_den_degrad*chi*vol_T/n_tumor) /(vol_D - vol_T);
  //the fluctuation of cell volume is small; not consider the effect of degradation by invasive cell

  return temp_rho;
}



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
  cell_density[ct_index] = 0.0;

  n_pro = 1; n_tumor = 1;
  n_nec = n_qui = 0;
  n_inva = 0;
  n_pro_acti = 0;

  prolif_list.push_back(ct_index);
  cell_close_edge_index[ct_index] = close_edge_cell(ct_index, ct_index, temp_dist2);

  vol_nec = vol_qui = 0.0; vol_inva = 0.0; vol_pro_acti = 0.0;
  vol_pro = cell_vol[ct_index];

  //for the den_ave_ECM
  sum_den_degrad += cell_density[ct_index]; //this should be updated for each proliferative cell and invasive cell
  double temp_neigh_dist[28];//not a good way, could be modified
  for(int i=0; i<cell_neig_ct[ct_index]; i++)
  {
	  cell_type[cell_neig[ct_index][i]] = 1;
	  cell_density[cell_neig[ct_index][i]] = 0.0;
	  n_pro++;
	  n_tumor++;
	  prolif_list.push_back(cell_neig[ct_index][i]);
	  cell_close_edge_index[cell_neig[ct_index][i]] = close_edge_cell(cell_neig[ct_index][i], cell_neig[ct_index][i], temp_neigh_dist[i]);
	  vol_pro += cell_vol[cell_neig[ct_index][i]];
	  sum_den_degrad += cell_density[cell_neig[ct_index][i]];
  }
  
  den_ave0 = Get_den_ave_ECM();
  den_ave_ECM = den_ave0;


  Rad_T = 1.0;
 

  //############### Can not think of anything else...
}


//given the index of the invasive cell, and its motility
//make the movement if it is allowed (based on degradation and rho)
//return the cell index it moved to...
//
//It is easy to modify this if there are other source that bias the motion of the invasive cells
int cell_move(int cell_index, double degradation, int motility)
{
  
  int temp_index = cell_index;

  //implement the trial moves
  //this could be an actual move, or just a degradation of the ECM
  for(int i=0; i<=motility; i++)
    {
      //check all neighbors
      double degrad_level_min = 10000000.0;
      double dlevel;
      double dist; //dist to the tumor center...
      double mix_max = -10000000.0;
      double mix;
      int neig_index;

    
      for(int j=0; j<cell_neig_ct[temp_index]; j++)
	{
	  if(cell_type[cell_neig[temp_index][j]] == 0)
	    {

	      dlevel = cell_density[cell_neig[temp_index][j]]*den_ave_ECM/den_ave0 - degradation;
	      dist = radial_dist(cell_neig[temp_index][j]);

	      if(dlevel<=0)
		{
		  mix = fabs(dlevel) + dist/Rmax;//mix the influence of optimum degradation direction and optimum direction for oxygen and nutrient supply

		  if(mix > mix_max)
		    {
		      mix_max = mix;
		      neig_index = cell_neig[temp_index][j];
		    }

		  if(dlevel<degrad_level_min) degrad_level_min = dlevel;
		}
	      
	      if(dlevel<0) cell_density[cell_neig[temp_index][j]] = 0;
	      else cell_density[cell_neig[temp_index][j]] = dlevel;
	      
	    }
	}

      //if the degradation is sufficent, the cell moved to neig_index
      if(degrad_level_min<=0) 
	{
	  cell_type[temp_index] = 0; //the original place turns back into a normal tissue cell with rho = 0;
	  temp_index = neig_index;
	  cell_type[temp_index] = 100;
	  cell_density[temp_index] = 0.0;
	}
      
      
      
    }
  
  return temp_index;
  
}


double Get_Max(double a, double b)
{
  if(a>=b) return a;
  else return b;
}

double Get_Min(double a, double b)
{
  if(a>=b) return b;
  else return a;
}


//get the growth advantage by comparing the local proliferative neighbors
double Get_delta_Lt(int edge_index)
{
  double L0 = length/pow(n_cell, 1.0/(double)n_dim); //characteristics length of a single cell
  //for normalization purpose

  double sum_Lt = 0.0;
  int counter = 0;

  for(int i=0; i<cell_neig_ct[edge_index]; i++)
    {
      int temp_index = cell_neig[edge_index][i];

      if(cell_type[temp_index] == 1)
	{
	  sum_Lt += radial_dist(temp_index);
	  
	  counter++;
	}
    }

  double delta = radial_dist(edge_index) - sum_Lt/(double)counter;

  return delta/L0;

  
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
	  int sign = 0;
      //cout<<" ***** Working with cell = "<<cell_index<<endl;
	  if(cell_status[cell_index]==0)//check mutation of dormant cell
	  {
		  double pt_temp = (double)(rand()%RAND_MAX)/(double)RAND_MAX;
		  if(pt_temp <= p_t)//determine whether trasition of status occurs
		  {
			  double frac_temp = (double)(rand()%RAND_MAX)/(double)RAND_MAX;
			  if(frac_temp<=frac_gentle)//gentle transition
				  {
					  cell_status[cell_index]=1;
				  }
			  else //malignant transition
				  {
					  cell_status[cell_index]=2;
				  }
			  n_pro_acti ++;
			  vol_pro_acti += cell_vol[cell_index];
		  }
	  }

	
      if(cell_type[cell_close_edge_index[cell_index]] == 0)
	{
	  edge_index = cell_close_edge_index[cell_index];
	  Le = dist_cell(cell_index, edge_index); //
	  //cout<<"edge_index0 = "<<edge_index<<endl;
	  //cout<<"Le = "<<Le<<endl;
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

	  //cout<<"Lt = "<<Lt<<endl;
	  
	  //using multiplication instead of pow, which is more costy


	  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	  //only when the tumor is sufficiently large, this cretia works;
	  //for very small tumors, every cell gets sufficient nutrients.
	  //check if there is space to place a daughter cell
	  double Delta_p = 1;
	  double B_d = 1;
	  double L_d1 = 1;
	  double b1;

	  if(n_pro_old>50)
	    {
	     
	      for(int j=0; j<(n_dim-1); j++)
		{
		  b1 = b*(2.0-1.0*(double)n_pro_acti/(double)n_pro);//set b as a function of the fraction of actively divided cells
		  Delta_p = Delta_p*Le;
		  B_d = B_d*b1;
		  L_d1 = L_d1*Lt;
		}
	      Delta_p = Delta_p*Le;
	      B_d = B_d*b1;
	      L_d1 = L_d1*B_d;
	      
	      //cout<<"Delta_p = "<<Delta_p<<"   L_d1 = "<<L_d1<<endl; 
	    }
	  else
	    {
	      Delta_p = Lt;
	      L_d1 = 2.0*Lt;
	    }
	  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	  if(Delta_p<=L_d1)
	    {
	      //- if so, check whether divide using p0, using R, Rmax
	      //bd_index = close_bd_cell2(cell_index, edge_index);
	      //this is for the growth in confined container, not used for the effects of pressure...
	      
	      

	      R = radial_dist(cell_index);
	      //Rmaxt = radial_dist(bd_index);

	      //cout<<"here2"<<endl;
	      //double ratio = R/Rmaxt;
	      //if(ratio>=1) ratio = 1;

	      //#######################################################
	      
		  double dp = wp*(den_ave_ECM/den_ave0-1);
	      
	      

	      double delta_Lt = Get_delta_Lt(edge_index)*den_ave_ECM/den_ave0;//this accounts for growth advantage due to surface fluctuations

	      p_d = Get_Max(0, p0*((1-dp-cell_density[edge_index]*den_ave_ECM/den_ave0)+delta_Lt)); //probability of division
	      //the first term takes into account long range elastic interaction, i.e.,through bd
	      //the second term takes into accout short range interaction, i.e., through ECM
	      //the rule for enhanced growth has to be local...

	      p_temp = (double)(rand()%RAND_MAX)/(double)RAND_MAX;

	      //cout<<"here3"<<endl;

	      //the cell will divide, otherwise nothing is changed...
	      //!!!!!!!!!
	      //if it divides, whether it divides into an invasive cell or not
	      if(p_temp < p_d && cell_type[edge_index] != -1&&cell_status[cell_index]!=0)
		{

		  double p_temp2 = (double)(rand()%RAND_MAX)/(double)RAND_MAX;

		  double degrad_flag = 1;

		  //double temp_p_inva = Get_Min(p_inva*(den_ave_ECM/den_ave0)*(den_ave_ECM/den_ave0), 1.0);
		  //################################################
		  //everything associated with invasion is proportional to the square of the pressure


		  //first check if the number of local number satisfy the dettaching condition
		  int dettach_flag = 0;

		  for(int l=0; l<cell_neig_ct[edge_index]; l++)
		    {
		      if(cell_type[cell_neig[edge_index][l]] == 1)
			{
			  dettach_flag++;
			}
		    }
		  //the detach flag should smaller than adhesion for an invasive to truely invade...
		  //also this is enhanced by the growth tip

		  //we could have a invasive cell....
		  if(p_temp2 < p_inva && dettach_flag < adhesion*(1+delta_Lt))
		    {
		      //first, get the invasive phenotype and see whether it can leave the tumor mass...
		      //if not, this invasive cell will just behave like another proliferative cell
		      int rho_neig_index;
			  
			  //now get the phenotype...
		      double temp_dg = (double)(rand()%RAND_MAX)/(double)RAND_MAX;
		      //temp_dg = temp_dg*degradation_max;
		      //here we could potential make it depndent on pressue,
		      //but the comprison is still wrt the intial cell_density, this favors invasions at high pressure..., i.e.,
		      temp_dg = temp_dg*degradation_max*(den_ave_ECM/den_ave0)*(den_ave_ECM/den_ave0);
              double temp_mo1 = (double)(rand()%RAND_MAX)/(double)RAND_MAX;
		      int temp_mo = (int)floor(temp_mo1*motility_max*(den_ave_ECM/den_ave0)*(den_ave_ECM/den_ave0));
			  double degrad_level_min = 10000000.0;
              double dlevel;
              double dist_center; //dist to the tumor center...
              double mix_max = -10000000.0;
              double mix;
              for(int l=0; l<cell_neig_ct[edge_index]; l++)
				{ 
				if(cell_type[cell_neig[edge_index][l]] == 0)
				{

	            dlevel = cell_density[cell_neig[edge_index][l]]*den_ave_ECM/den_ave0 - temp_dg;
	            dist_center = radial_dist(cell_neig[edge_index][l]);

	             if(dlevel<=0)
		        {
		         mix = fabs(dlevel) + dist_center/Rmax;//mix the influence of optimum degradation direction and optimum direction for oxygen and nutrient supply

		         if(mix > mix_max)
		        {
		          mix_max = mix;
		          rho_neig_index = cell_neig[edge_index][l];
		         }

		         if(dlevel<degrad_level_min) degrad_level_min = dlevel;
		         }
	      
	           }
	           }

      
              if(degrad_level_min>0) degrad_flag = 0;//the invasive cell doesn't gain enough level of degradation...
		      //this pre-selection based on cell-cell adhsion is Extremely Important
		     
		      else
			{
			  //this means it is really an invasive cell...
			  degrad_flag = 1;

			  sum_den_degrad += cell_density[edge_index]*den_ave_ECM/den_ave0;

			  //the following steps are taken care of in cell_move
			  cell_type[edge_index] = 100;
			  cell_density[edge_index] = 0.0; //eat up the ECM

			  //make the move...
			  //the first time, make the largest jump
			  int temp_inva_index = cell_move(edge_index, temp_dg, temp_mo);

			  
			  
			  n_inva++;
			  //update the n_pro and n_pro_old, the new prolif cell is treated next round...
			  //this makes sure we don't  go into a dead loop 
			  
			  //insert this cell to the list, 
			  invasive_list.push_back(temp_inva_index);

			  inva_motility.push_back(temp_mo);
			  inva_degradation.push_back(temp_dg);
			  
			  vol_inva += cell_vol[temp_inva_index];
			  
			  
			}




		    }
		  if(p_temp2 >= p_inva || degrad_flag == 0)
		    {
		      //cout<<"A new prolif cell is found ..."<<endl;
		      sum_den_degrad += cell_density[edge_index]*den_ave_ECM/den_ave0;
		      //keep track of this for computing the average ECM density

		      cell_type[edge_index] = 1;

		      cell_density[edge_index] = 0.0; //eat up the ECM
		      
		      n_pro++;
			  n_pro_acti++;
		      //update the n_pro and n_pro_old, the new prolif cell is treated next round...
		      //this makes sure we don't  go into a dead loop 
		      
		      prolif_list.push_back(edge_index);
		      
		      vol_pro += cell_vol[edge_index];
		      vol_pro_acti += cell_vol[edge_index];
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
		      if(cell_status[cell_index] == 1)
			  {
				  double pt_temp3 = (double)(rand()%RAND_MAX)/(double)RAND_MAX;
				  if(pt_temp3 <= frac_gentle)
				  {
					  cell_status[edge_index] = 1;
				  }
				  else
				  {
					  cell_status[edge_index] = 2;
				  }
			  }
			  if(cell_status[cell_index] == 2)
				  cell_status[edge_index] = 2;
		      //cout<<"The new prolif cell is added ..."<<endl;
		    }
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
		  sign = 1;//the prolif cell should be erased
	      //erase the prolif cell from the prolif list...
	      n_pro --;
		  
	      n_pro_old --; //this is needed in both case...
	      
	      vol_pro = vol_pro - cell_vol[cell_index];
		  if(cell_status[cell_index]!=0)
		  {
			  n_pro_acti --;
			  vol_pro_acti = vol_pro_acti - cell_vol[cell_index];
		  }

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
      sign = 1;
	  cell_type[cell_index] = 2; 
	  n_qui ++;
	  quiscent_list.push_back(cell_index);
	  vol_qui += cell_vol[cell_index];
	  
	  //erase the prolif cell from the prolif list...
	  n_pro --;
	  n_pro_old --; //this is needed in both case...
	  
	  vol_pro = vol_pro - cell_vol[cell_index];
	  if(cell_status[cell_index]!=0)//proliferative cells could turn into quiscent before they actively divide
		  {
			  n_pro_acti --;
			  vol_pro_acti = vol_pro_acti - cell_vol[cell_index];
		  }
	  prolif_list.erase(prolif_list.begin()+i);
	  
	  //################# just put a checker here...
	  if(n_pro != prolif_list.size())
	    {
	      cout<<"Counting is not correct for n_pro and prolif_list.size()! Abort! "<<endl;
	      exit(1);
	    }
	  
	}
	
	if(cell_status[cell_index]==1&&sign==0)//response of immune system to gentle mutation
	  {
		  double res_temp1=(double)(rand()%RAND_MAX)/(double)RAND_MAX;
		  if(res_temp1 <= r_gentle)
			  {
				  double kill_temp = (double)(rand()%RAND_MAX)/(double)RAND_MAX;
				  if(kill_temp <= frac_kill)//kill the actively dividing cell
				  {
					   cell_type[cell_index] = 0;
					   cell_density[cell_index] = 0.0;
					   for(int k=0; k<n_dim; k++)
						{
							tumor_ct_coord[k] = n_tumor*tumor_ct_coord[k] - cell_coord[cell_index][k];
							tumor_ct_coord[k] = tumor_ct_coord[k]/(double)(n_tumor - 1);
						}
					   n_pro --;
					   n_pro_old --;
					   n_pro_acti --;
					   vol_pro = vol_pro - cell_vol[cell_index];
					   vol_pro_acti = vol_pro_acti - cell_vol[cell_index];
					   n_tumor --;
					   prolif_list.erase(prolif_list.begin()+i);
					   if(n_pro != prolif_list.size())
						{
						 cout<<"Counting is not correct for n_pro and prolif_list.size()! Abort! "<<endl;
						 exit(1);
						}
				  }
				  else
				  {
					cell_status[cell_index] = 0;
					n_pro_acti --;
					vol_pro_acti = vol_pro_acti - cell_vol[cell_index];
				  }
		       }
	  }
   if(cell_status[cell_index]==2&&sign==0)//response of immune system to malignant mutation
	  {
		  double res_temp2=(double)(rand()%RAND_MAX)/(double)RAND_MAX;
		  if(res_temp2 <= r_mal)
			  {
				  double kill_temp2 = (double)(rand()%RAND_MAX)/(double)RAND_MAX;
				  if(kill_temp2 <= frac_kill)//kill the actively dividing cell
				  {
					   cell_type[cell_index] = 0;
					   cell_density[cell_index] = 0.0;
					   for(int k=0; k<n_dim; k++)
						{
							tumor_ct_coord[k] = n_tumor*tumor_ct_coord[k] - cell_coord[cell_index][k];
							tumor_ct_coord[k] = tumor_ct_coord[k]/(double)(n_tumor - 1);
						}
					   n_pro --;
					   n_pro_old --;
					   n_pro_acti --;
					   vol_pro = vol_pro - cell_vol[cell_index];
					   vol_pro_acti = vol_pro_acti - cell_vol[cell_index];
					   n_tumor --;
					   prolif_list.erase(prolif_list.begin()+i);
					   if(n_pro != prolif_list.size())
						{
						 cout<<"Counting is not correct for n_pro and prolif_list.size()! Abort! "<<endl;
						 exit(1);
						}
				  }
				  else
				  {
					cell_status[cell_index] = 0;
					n_pro_acti --;
					vol_pro_acti = vol_pro_acti - cell_vol[cell_index];
				  }
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
	  double a1;
	  for(int j=0; j<(n_dim-1); j++)
	    {
	      a1 = a*(1.6-0.6*(double)n_pro_acti/(double)n_pro);//set a as a function of the fraction of actively divided cells
		  Delta_p = Delta_p*Le;
	      A_d = A_d*a1;
	      L_d1 = L_d1*Lt;
	    }
	  Delta_p = Delta_p*Le;
	  A_d = A_d*a1;
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


//these are the most aggressive phenotypes, update first, i.e., before update_proliferative
//at this stage, we consider the invasive cells do not divide, could program this at a later time
void Update_Invasive()
{
  //loop over all invasive cells, 
  //make the move ...
  //!!!! The phenotype could also be updated, but now leave it alone
  for(int i=0; i<n_inva; i++)
    {
      int temp_index = invasive_list[i];

      //search its neighbors to see whether it is movable...
      //if so, keep it invasive, otherwise, turn it to proliferative
      int inva_flag = 0;
      for(int j=0; j<cell_neig_ct[temp_index]; j++)
	{
	  if(cell_type[cell_neig[temp_index][j]] == 0 || cell_type[cell_neig[temp_index][j]] == -1)
	    {
	      //either at tumor edge or the boundary of growth permitable regions
	      inva_flag = 1;
	      break;
	    }
	}
      
      //have room to move, make the movement
      if(inva_flag == 1)
	{
	  int inva_index = cell_move(invasive_list[i], inva_degradation[i], inva_motility[i]);
      
	  invasive_list[i] = inva_index;
	  
	  //update the invasive cells volume ... this is not very important, but just make consistent...
	  vol_inva = vol_inva - cell_vol[temp_index] + cell_vol[inva_index];
	}
      else //turn it back to proliferative,
	{
	  cell_type[temp_index] = 1;
	  cell_density[temp_index] = 0.0; //eat up the ECM		  
	  
	  n_pro++;
	  //update the n_pro and n_pro_old, the new prolif cell is treated next round...
	  //this makes sure we don't  go into a dead loop 
	  
	  prolif_list.push_back(temp_index);
	  
	  vol_pro += cell_vol[temp_index];
		      
	  //now we need to assing a neighbor to this new cell...
	  double temp_d;
	  int tmp_neig_new_cell = close_edge_cell(temp_index, temp_index, temp_d);
	  
	  if(tmp_neig_new_cell != -1)
	    {
	      cell_close_edge_index[temp_index] = tmp_neig_new_cell;
	    }
	  else
	    {
	      //it means this cell is now on the tumor boundary...
	      cell_close_edge_index[temp_index] = temp_index;
	    }
	  	 

	  //now remove it from the invasive list...

	  invasive_list.erase(invasive_list.begin()+i);
	  n_inva--;
	  vol_inva -= cell_vol[temp_index];
	  
	  //################# just put a checker here...
	  if(n_inva != invasive_list.size())
	    {
	      cout<<"Counting is not correct for n_inva and invasive_list.size()! Abort! "<<endl;
	      exit(1);
	    }
	  
	  
	 
	}
    }
}


//Get the tumor radius, asphreicity and specific surface
//this is a complete update each time the function is called...
//search for all proliferative cells, find those on the boundary...
void tumor_stat()
{
 
  //for asphericity
  Rad_Tmax = 0.0;
  Rad_Tmin = 1000000.0;


  //for tumor radius
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

	      tmp_dist = sqrt(tmp_dist);

	      if(tmp_dist<Rad_Tmin) Rad_Tmin = tmp_dist;
	      if(tmp_dist>Rad_Tmax) Rad_Tmax = tmp_dist;

	      edge_dist += tmp_dist;
	      edge_ct ++;

	      continue;
	    }
	}
    }

  //for tumor specific surface
  double L0 = length/pow(n_cell, 1.0/(double)n_dim); //characteristics length of a single cell

  //tumor radius
  Rad_T = edge_dist/(double)edge_ct;

  //asphericy
  alpha_T = Rad_Tmax/Get_Max(Rad_T, L0);

 
  double s0 = (double)n_dim/Rad_T; //s for circles, for normalization purpose
  
  //specific surface
  ss_T = ((double)edge_ct/((n_nec+n_qui+n_pro)*L0))/s0;
}


//return the volume of the tumor...
double tumor_vol()
{
  return (vol_pro+vol_qui+vol_nec+vol_inva);
}


//print the 
void Print_Degrad(int t)
{
  //ofstream fout;
  //fout.open("D_Frames.dat");
  
  if(t%t1==0)
    {
      int findex = t/t1;
      string s;
      stringstream out;
      out << findex;
      s = out.str();
      
      string sf = "D_Frames.";
      sf.append(s);

      int sizes = sf.size();
      char* filename = new char[sizes];
      for(int i=0; i<sizes; i++)
	filename[i] = sf[i];
      
      ofstream outfile0;
      outfile0.open(filename);

      outfile0<<"#################################################################################"<<endl;
      outfile0<<"#################################################################################"<<endl;
      outfile0<<"# t="<<t<<endl;
      

      for(int i=0; i<n_cell; i++)
	{
	  if(cell_density[i]*den_ave_ECM/den_ave0<0.00001)
	    {
	      for(int j=0; j<n_dim; j++)
		outfile0 <<cell_coord[i][j] << "  ";
	      outfile0 <<endl;
	    }
	}

      outfile0<<"#################################################################################"<<endl;
      outfile0<<"#################################################################################"<<endl;
      outfile0<<"#################################################################################"<<endl;
      outfile0<<"#################################################################################"<<endl;
      outfile0<<"#################################################################################"<<endl;
      outfile0<<"#################################################################################"<<endl;

      outfile0.close();
    }

}

//we print out the tumor at specified time point...
//this is not very concise in terms of codeing ....
void Print_Tumor(int t)
{
  //ofstream fout;
  //fout.open("T_Frames.dat");
  

  if(t%t1==0)
    {
      int findex = t/t1;
      string s;
      stringstream out;
      out << findex;
      s = out.str();

      string sf = "T_Frames.";
      sf.append(s);
      
      int sizes = sf.size();
      char* filename = new char[sizes];
      for(int i=0; i<sizes; i++)
	filename[i] = sf[i];


      ofstream outfile0;
      outfile0.open(filename);

      outfile0<<"#################################################################################"<<endl;
      outfile0<<"#################################################################################"<<endl;
      outfile0<<"#################################################################################"<<endl;
      outfile0<<"#################################################################################"<<endl;
      outfile0<<"#################################################################################"<<endl;
      outfile0<<"#################################################################################"<<endl;
      
      outfile0<<"# t = "<<t<<endl;
      outfile0<<"#tumor center = "<<tumor_ct_coord[0]<<" "<<tumor_ct_coord[1]<<" "<<tumor_ct_coord[2]<<endl;
      outfile0 <<"# n_pro_acti = "<<n_pro_acti<<"   vol_pro_acti = "<<vol_pro_acti<<endl;
	  outfile0 <<"# n_pro_inac = "<<n_pro-n_pro_acti<<"   vol_pro_inac = "<<vol_pro-vol_pro_acti<<endl;
      outfile0 <<"# n_qui = "<<n_qui<<"   vol_qui = "<<vol_qui<<endl;
      outfile0 <<"# n_nec = "<<n_nec<<"   vol_nec = "<<vol_nec<<endl;
      outfile0 <<"# n_inv = "<<n_inva<<"   vol_inv = "<<vol_inva<<endl;
      

      int index;

      for(int i=0; i<n_pro; i++)//print actively dividing proliferative cells
	{
	  index = prolif_list[i];
	  if(cell_status[index]!=0)
	  {
	  for(int j=0; j<n_dim; j++)
	    outfile0 <<cell_coord[index][j] << "  ";
	  outfile0 <<endl;
	  }
	}	 
      outfile0<<"################################"<<endl<<endl;

	   for(int i=0; i<n_pro; i++)//print dormant proliferative cells
	{
	  index = prolif_list[i];
	  if(cell_status[index]==0)
	  {
	  for(int j=0; j<n_dim; j++)
	    outfile0 <<cell_coord[index][j] << "  ";
	  outfile0 <<endl;
	  }
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

      for(int i=0; i<n_inva; i++)
	{
	  index = invasive_list[i];
	  
	  for(int j=0; j<n_dim; j++)
	    outfile0 <<cell_coord[index][j] << "  ";
	  outfile0 <<endl;
	}	 
      outfile0<<"################################"<<endl<<endl;
      
      outfile0.close();

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

  Init_CellDensity(density_type);
  //~~~~~~~~~~~~~~~~~~~~~~~~~
  Init_Tumor();

  ofstream outRad;
  outRad.open("tumor_rad.txt");  
  ofstream outVol;
  outVol.open("tumor_vol.txt");
  ofstream outVol_type;
  outVol_type.open("tumor_vol_type.txt");
  ofstream outDen;
  outDen.open("den_ave.txt");

  ofstream outSS;
  outSS.open("spe_surface.txt");
  ofstream outAlpha;
  outAlpha.open("alpha.txt");

  for(int t=0; t<iterations; t++)
    {
      cout<<"At T = "<<t<<"..."<<endl;

      Update_Invasive();
      
      Update_Quiscent();

      //cout<<"here1"<<endl;
      
      Update_Proliferative();

      //cout<<"here2"<<endl;

      den_ave_ECM = Get_den_ave_ECM();

      Print_Tumor(t);
      Print_Degrad(t);

      cout<<"n_tumor = "<<n_tumor<<endl;
      cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

      //compute tumor stat
      tumor_stat();

      outRad << t <<"    "<<Rad_T<<endl;
      outVol << t <<"    "<<tumor_vol()/vol_growth_domain<<endl;
	  outVol_type << t <<"    "<<vol_pro_acti/vol_growth_domain<<"    "<<(vol_pro-vol_pro_acti)/vol_growth_domain<<"    "<<vol_qui/vol_growth_domain<<"    "<<vol_nec/vol_growth_domain<<"    "<<vol_inva/vol_growth_domain<<endl;
      outDen << t <<"    "<<den_ave_ECM/den_ave0<<endl;
      outSS << t <<"    "<<ss_T<<endl;
      outAlpha << t <<"    "<<alpha_T<<endl;
      
    }


  outRad.close();
  outVol.close();
  outVol_type.close();
  outDen.close();
}
