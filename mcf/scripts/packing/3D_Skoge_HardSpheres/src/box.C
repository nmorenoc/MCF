/* 
   Packing of hard spheres via molecular dynamics
   Developed by Monica Skoge, 2006, Princeton University
   Contact: Aleksandar Donev (adonev@math.princeton.edu) with questions
   This code may be used, modified and distributed freely.
   Please cite:
   
   "Packing Hyperspheres in High-Dimensional Euclidean Spaces"
   	M. Skoge, A. Donev, F. H. Stillinger and S. Torquato, 2006
	
   if you use these codes.	
*/

#include "box.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iomanip>

//==============================================================
//==============================================================
//  Class Box: Fills box with hardspheres to given packing fraction
//  and evolves spheres using molecular dynamics!
//==============================================================
//==============================================================


//==============================================================
// Constructor
//==============================================================
box::box(int N_i, double r_i, double growthrate_i, double maxpf_i):
  r(r_i),          
  N(N_i),
  growthrate(growthrate_i),
  h(N_i+1),
  maxpf(maxpf_i)
{
  ngrids = Optimalngrids(maxpf);
  cells.set_size(ngrids);
 
  s = new sphere[N];
  binlist = new int[N];
  x = new vector<DIM>[N];        
  h.s = s;

  gtime = 0.;
  rtime = 0.;
  ncollisions = 0;
  ntransfers = 0;
  nchecks = 0;
  neventstot = 0;
  ncycles = 0;
  xmomentum = 0.; 
  pressure = 0.;

  cells.set_size(ngrids);
  cells.initialize(-1);      // initialize cells to -1
  srandom(::time(0));        // initialize the random number generator
  for (int i=0; i<N; i++)    // initialize binlist to -1
    binlist[i] = -1;
  
  time(&start);
}


//==============================================================
// Destructor
//==============================================================
box::~box() 
{
  delete[] s;
  delete[] binlist;
  delete[] x;
}


//==============================================================
// ReadFile
//==============================================================
void box::ReadPositions(const char* filename)
{  
  // open file to read in arrays
  std::ifstream infile(filename);
  
  infile.ignore(256, '\n');  // ignore the dim line
  infile.ignore(256, '\n');  // ignore the #sphere 1 line
  infile.ignore(256, '\n');  // ignore the #sphere line
  infile.ignore(256, '\n');  // ignore the diameter line
  infile.ignore(1000, '\n'); // ignore the 100 010 001 line
  infile.ignore(256, '\n');  // ignore the T T T line

  for (int i=0; i<N; i++)
    for (int k=0; k<DIM; k++)
      infile >> s[i].x[k];

  infile.close();
}


//==============================================================
// Recreates all N spheres at random positions
//==============================================================
void box::RecreateSpheres(const char* filename, double temp)
{
  ReadPositions(filename);  // reads in positions of spheres
  VelocityGiver(temp);      // gives spheres initial velocities
  AssignCells();            // assigns spheres to cells
  SetInitialEvents();
}


//==============================================================
// Creates all N spheres at random positions
//==============================================================
void box::CreateSpheres(double temp)
{
  int Ncurrent = 0;
  for(int i=0; i<N; i++)
    {
      CreateSphere(Ncurrent);
      Ncurrent++;
    }
  if (Ncurrent != N)
    std::cout << "problem! only made " << Ncurrent << " out of " << N << " desired spheres" << std::endl;
  
  VelocityGiver(temp);
  SetInitialEvents();
}

   
//==============================================================
// Creates a sphere of radius r at a random unoccupied position 
//==============================================================
void box::CreateSphere(int Ncurrent)
{
  int keeper;    // boolean variable: 1 means ok, 0 means sphere already there
  int counter = 0;   // counts how many times sphere already exists
  vector<DIM> xrand;  // random new position vector
  double d = 0.;

  while (counter<1000)
    {
      keeper = 1;
      
      for(int k=0; k<DIM; k++) 
	xrand[k] = ((double)random()/(double)RAND_MAX)*SIZE;
      
      for (int i=0; i<Ncurrent; i++)  // need to check nearest image!
	{
	  d=0.;
	  for (int k=0; k<DIM; k++)
	    {
	      if ((xrand[k] - s[i].x[k])*(xrand[k] - s[i].x[k]) > SIZE*SIZE/4.)
		{
		  if (xrand[k] > SIZE/2.)  // look at right image
		    d += (xrand[k] - (s[i].x[k]+SIZE))*
		      (xrand[k] - (s[i].x[k]+SIZE));
		  else                     // look at left image
		    d += (xrand[k] - (s[i].x[k]-SIZE))*
		      (xrand[k] - (s[i].x[k]-SIZE));
		}
	      else
		d += (xrand[k] - s[i].x[k])*(xrand[k] - s[i].x[k]);
	    }
	  	   
	  if (d <= 4*r*r)
	    {
	      keeper = 0;
	      counter++;
	      break;
	    }
	}
      if (keeper == 1)
	break;
    }
  if (counter >= 1000)
    {
      std::cout << "counter >= 1000" << std::endl;
      exit(-1);
    }

  // now convert xrand into index vector for cells
  vector<DIM,int> cell;
  cell = vector<DIM>::integer(xrand*((double)(ngrids))/SIZE);
  
  s[Ncurrent] = sphere(Ncurrent, xrand, cell, gtime);
  
  //first check to see if entry at cell
  if (cells.get(cell) == -1) //if yes, add Ncurrent to cells gridfield
    cells.get(cell) = Ncurrent;

  else  // if no, add i to right place in binlist
    {  
      int iterater = cells.get(cell); // now iterate through to end and add Ncurrent
      int pointer = iterater;
      while (iterater != -1)
	{
	  pointer = iterater;
	  iterater = binlist[iterater];
	}
      binlist[pointer] = Ncurrent;
    }

  Ncurrent++;
}


//==============================================================
// Assign cells to spheres read in from existing configuration
//==============================================================
void box::AssignCells()
{
  for (int i=0; i<N; i++)
    {
      // now convert x into index vector for cells
      vector<DIM,int> cell;
      cell = vector<DIM>::integer(s[i].x*((double)(ngrids))/SIZE);
      s[i].cell = cell;
      
      //first check to see if entry at cell
      if (cells.get(cell) == -1) //if yes, add Ncurrent to cells gridfield
	cells.get(cell) = i;
      
      else  // if no, add i to right place in binlist
	{  
	  int iterater = cells.get(cell); // now iterate through to end and add Ncurrent
	  int pointer = iterater;
	  while (iterater != -1)
	    {
	      pointer = iterater;
	      iterater = binlist[iterater];
	    }
	  binlist[pointer] = i;
	}
    }
}

	
//==============================================================
// Velocity Giver, assigns initial velocities from Max/Boltz dist.
//==============================================================
void box::VelocityGiver(double T)
{
  for (int i=0; i<N; i++)
    {
      for (int k=0; k<DIM; k++)
	{
	  if (T==0.)
	    s[i].v[k] = 0.;
	  else
	    s[i].v[k] = Velocity(T);
	}
    }
}


//==============================================================
// Velocity, gives a single velocity from Max/Boltz dist.
//==============================================================
double box::Velocity(double T)
{
  double rand;                       // random number between -0.5 and 0.5
  double sigmasquared = T;    // Assumes M = mass of sphere = 1
  double sigma = sqrt(sigmasquared); // variance of Gaussian
  double stepsize = 1000.;           // stepsize for discretization of integral
  double vel = 0.0;                  // velocity
  double dv=sigma/stepsize;
  double p=0.0;
  
  rand = (double)random() / (double)RAND_MAX - 0.5;
  if(rand < 0) 
    {
      rand = -rand;
      dv = -dv;
    }
  
  while(fabs(p) < rand) // integrate until the integral equals rand
    {
      p += dv * 0.39894228 * exp(-vel*vel/(2.*sigmasquared))/sigma;
      vel += dv;
    }
  return vel;
}


//==============================================================
// Finds next events for all spheres..do this once at beginning
//==============================================================
void box::SetInitialEvents()
{
  for (int i=0; i<N; i++)  // set all events to checks
    {
      event e(gtime, i, INF); 
      s[i].nextevent = e;
      h.insert(i);
    }
}


//==============================================================
// Finds next event for sphere i 
//==============================================================
event box::FindNextEvent(int i)
{
  event t = FindNextTransfer(i);
  event c = FindNextCollision(i);

  if ((c.time < t.time)&&(c.j == INF)) // next event is check at DBL infinity
    return c;
  else if (c.time < t.time) // next event is collision!
    {
      CollisionChecker(c); 
      return c; 
    }
  else // next event is transfer!
    return t;  
} 


//==============================================================
// Checks events of predicted collision partner to keep collisions
// symmetric
//==============================================================
void box::CollisionChecker(event c)
{
  int i = c.i;
  int j = c.j;
  event cj(c.time,j,i,c.v*(-1));

  // j should have NO event before collision with i!
  if (!(c.time  < s[j].nextevent.time))
    std::cout << i << " " <<  j << " error collchecker, s[j].nextevent.time= " << s[j].nextevent.time << " " << s[j].nextevent.j << ", c.time= " << c.time << std::endl;
  
  int k = s[j].nextevent.j; 
  if ((k < N) && (k!=i)) // j's next event was collision so give k a check
    s[k].nextevent.j = INF;
  
  // give collision cj to j
  s[j].nextevent = cj;
  h.upheap(h.index[j]);
}


//==============================================================
// Find next collision for sphere i 
//==============================================================
event box::FindNextTransfer(int i)
{
  double ttime = dblINF;  
  int wallindex = INF;   // -(k+1) left wall, (k+1) right wall

  vector<DIM> xi = s[i].x + s[i].v*(gtime - s[i].lutime);
  vector<DIM> vi = s[i].v;

  for (int k=0; k<DIM; k++)
    {
      double newtime;
      if (vi[k]==0.) 
	newtime= dblINF;
      else if (vi[k]>0)  // will hit right wall
	{
	  newtime = ((double)(s[i].cell[k]+1)*SIZE/((double)(ngrids))
		     - xi[k])/(vi[k]);
	  if (newtime < 0)
	    std::cout << "error in FindNextTransfer right newtime < 0 " << k << std::endl;
	  if (newtime<ttime)
	    {
	      wallindex = k+1;
	      ttime = newtime;
	    }
	}
      else if (vi[k]<0)  // will hit left wall
	{
	  newtime = ((double)(s[i].cell[k])*SIZE/((double)(ngrids)) 
		     - xi[k])/(vi[k]);
	  if (newtime < 0)
	    {
	      if (newtime > -10.*DBL_EPSILON) // this should happen only when reading in a configuration and spheres is on left boundary moving left
		newtime = 0.;
	      else
		std::cout << "error in FindNextTransfer left newtime < 0 " << k << std::endl;
	    }
	  if (newtime<ttime)
	    {
	      wallindex = -(k+1);
	      ttime = newtime;
	    }
	}
    }

  if (ttime < 0)
    {
      std::cout << "error in FindNextTransfer ttime < 0" << std::endl;
      std::cout << i << std::endl;
      std::cout << xi << " " << s[i].x << std::endl;
      std::cout << vi << " " << s[i].v << std::endl;
      std::cout << s[i].cell << std::endl;
      exit(-1);
    }
  // make the event and return it
  event e = event(ttime+gtime,i,wallindex+DIM+N+1);
  return e;
}


//==============================================================
// Find next collision for sphere i 
//==============================================================
void box::ForAllNeighbors(int i, vector<DIM,int> vl, vector<DIM,int> vr,
			  neighbor& operation)
{
  vector<DIM,int> cell = s[i].cell;

  // now iterate through nearest neighbors
  vector<DIM, int> offset;          // nonnegative neighbor offset
  vector<DIM, int> pboffset;        // nearest image offset

   vector<DIM,int> grid;

   int ii ;

   grid=vl;
   while(1)
   {
     //if (vr[0] > 1)
     //std::cout << grid << "..." << cell+grid << "\n";
     for(int k=0; k<DIM; k++)
     {
        offset[k]=grid[k]+ngrids;  // do this so no negatives 	
        if (cell[k]+grid[k]<0) //out of bounds to left
          pboffset[k] = -1;
	else if (cell[k]+grid[k]>=ngrids) // out of bounds to right
	  pboffset[k] = 1;
        else
          pboffset[k] = 0;
     }     
     int j = cells.get((cell+offset)%ngrids);
     while(j!=-1)
       {
	 operation.Operation(j,pboffset);
	 j = binlist[j];
       }

     // A. Donev:     
     // This code makes this loop dimension-independent
     // It is basically a flattened-out loop nest of depth DIM
     for(ii=0;ii<DIM;ii++)
     {
       grid[ii] += 1;
       if(grid[ii]<=vr[ii]) break;
       grid[ii]=vl[ii];
     }
     if(ii>=DIM) break;
   }  
}


//==============================================================
// PredictCollision
//==============================================================
void box::PredictCollision(int i, int j, vector<DIM, int> pboffset, 
			     double& ctime, int& cpartner, 
			     vector<DIM, int>& cpartnerpboffset)
{
  double ctimej;
  
  if (i!=j)
    {	 
      ctimej = CalculateCollision(i,j,pboffset.Double())+gtime;
      
      if (ctimej < gtime)
	std::cout << "error in find collision ctimej < 0" << std::endl;
      
      if ((ctimej < ctime)&&(ctimej < s[j].nextevent.time))
	{
	  ctime = ctimej;
	  cpartner = j;
	  cpartnerpboffset = pboffset;
	}	
    }
}


//==============================================================
// Find next collision
//==============================================================
event box::FindNextCollision(int i)
{
  collision cc(i, this);
  
  vector<DIM, int> vl, vr;

  for (int k=0; k<DIM; k++)  // check all nearest neighbors
    {
      vl[k] = -1;
      vr[k] = 1;
    }
  
  ForAllNeighbors(i,vl,vr,cc);

  event e;
  if (cc.cpartner == i)  // found no collisions in neighboring cells
    {
      if (cc.ctime != dblINF)
	std::cout << "ctime != dblINF" << std::endl;
      e = event(dblINF,i,INF);  // give check at double INF
    }
  else
    e = event(cc.ctime,i,cc.cpartner,cc.cpartnerpboffset);

  return e;
}


//==============================================================
// Calculates collision time between i and image of j using quadratic formula
//==============================================================
double box::CalculateCollision(int i, int j, vector<DIM> pboffset)
{
// calculate updated position and velocity of i and j
  vector<DIM> xi = s[i].x + s[i].v*(gtime - s[i].lutime);
  vector<DIM> vi = s[i].v;
  vector<DIM> xj = s[j].x + pboffset*SIZE + s[j].v*(gtime - s[j].lutime);
  vector<DIM> vj = s[j].v;
  
  double r_now = r + gtime*growthrate;
  
  double A,B,C;
  A = vector<DIM>::norm_squared(vi - vj) - 4*growthrate*growthrate;
  B = vector<DIM>::dot(xi - xj, vi - vj) - 4*r_now*growthrate;
  C = vector<DIM>::norm_squared(xi - xj) - 4*r_now*r_now;

  if (C < -1E-12*2.*r_now)
    {
      std::cout << "error, " << i << " and " << j << " are overlapping at time "<< gtime << " and A, B, C = "  << A << " " << " " << B << " " << " " << C <<  std::endl;
      std::cout << "velocity i=  " << s[i].v << ", velocity j= " << s[j].v << ", gtime= " << gtime << ", det= " << B*B - A*C << std::endl;
      exit(-1);
    }
      
  return QuadraticFormula(A, B, C);
}


//==============================================================
// Quadratic Formula ax^2 + bx + c = 0
//==============================================================
 double box::QuadraticFormula(double a, double b, double c)
{
  double x = dblINF;
  double xpos;
  double xneg;
  double det = b*b - a*c;

  if (c <= 0.)
    {
      if(b < 0.) // spheres already overlapping and approaching
	{
	  //std::cout << "spheres overlapping and approaching" << std::endl;
	  //std::cout << "# events= " << neventstot << std::endl;
	  x = 0.;	
	}
    }
  else if (det > -10.*DBL_EPSILON)
    {
      if (det < 0.)  // determinant can be very small for double roots
	det = 0.;    
      if (b < 0.)
	x = c/(-b + sqrt(det));
      else if ((a < 0.)&&(b > 0.))
	x = -(b + sqrt(det))/a;
      else
	x = dblINF;
    }
  return x;
}


//==============================================================
// Returns first event
//==============================================================
void box::ProcessEvent()
{  
  neventstot++;
  // Extract first event from heap
  int i = h.extractmax();   
  event e = s[i].nextevent; // current event
  event f;                  // replacement event

  if ((e.j>=0)&&(e.j<N))  // collision!
    {
      ncollisions++;
      //std::cout << "collision between " << e.i << " and " << e.j << " at time " << e.time << std::endl;
      Collision(e);
      f = FindNextEvent(i);
      s[i].nextevent = f;
      h.downheap(1);
      if (f.time < e.time)
	{
	  std::cout << "error, replacing event with < time" << std::endl;
	  exit(-1);
	}
      
      /*
      if (f.time == e.time)
	{
	  std::cout << "replacing event with = time (it's ok)" << std::endl;
	  std::cout << "# events= " << neventstot << std::endl;
	}
      */

      // make sure collision was symmetric and give j a check
      if ((s[e.j].nextevent.j != i)||(s[e.j].nextevent.time != gtime))
	{
	  std::cout << "error collisions not symmetric" << std::endl;
	  std::cout << "collision between " << e.i << " and " << e.j << " at time " << e.time << std::endl;
	  std::cout << "but " << e.j << " thinks it has " << s[e.j].nextevent.j<< " "  << s[e.j].nextevent.time << std::endl;
	  exit(-1);
	}
      else  // give j a check
	s[e.j].nextevent.j = INF;
    }
  else if (e.j==INF)      // check!  
    {
      nchecks++;
      //std::cout << "check for " << e.i << " at time " << e.time << std::endl;
      f = FindNextEvent(i);
      s[i].nextevent = f;
      h.downheap(1);
    }
  else                    // transfer!
    {
      ntransfers++;
      //std::cout << "transfer for " << e.i << " at time " << e.time << std::endl;
      Transfer(e);
      f = FindNextEvent(i);
      s[i].nextevent = f;
      h.downheap(1);
      //r = FindNextEvent(i, e.j-N-DIM-1);
      if (f.time <= e.time)
	{
	  std::cout << "error after transfer, replacing new event with <= time" << " " << std::endl;
	  std::cout << "e.time= " << e.time << ", f.time= " << f.time << ", f.i= " << f.i << ", f.j= " << f.j << std::endl;
	  std::cout << "difference= " << e.time - f.time << std::endl;
	  exit(-1);
	}
    }
}


//==============================================================
// Processes a collision
//=============================================================
void box::Collision(event e)
{
  double ctime = e.time;
  int i = e.i;
  int j = e.j;
  vector<DIM,int> v = e.v;  // virtual image
  gtime = ctime;

  // Update positions and cells of i and j to ctime
  s[i].x += s[i].v*(gtime-s[i].lutime);
  s[j].x += s[j].v*(gtime-s[j].lutime);

  // Check to see if a diameter apart
  double r_now = r + gtime*growthrate;
  double distance = vector<DIM>::norm_squared(s[i].x - s[j].x- v.Double()*SIZE) - 4*r_now*r_now;
  if (distance*distance > 10.*DBL_EPSILON)
    std::cout << "overlap " << distance << std::endl;

  s[i].lutime = gtime;
  s[j].lutime = gtime;
  
  vector<DIM,double> vipar;          // parallel comp. vi
  vector<DIM,double> vjpar;          // parallel comp. vj
  vector<DIM,double> viperp;         // perpendicular comp. vi
  vector<DIM,double> vjperp;         // perpendicular comp. vj

  // make unit vector out of displacement vector
  vector<DIM,double> dhat;
  dhat = s[i].x - s[j].x - v.Double()*SIZE;  // using image of j!!
  double dhatmagnitude = sqrt(dhat.norm_squared());
  dhat /= dhatmagnitude;

  vipar = dhat*vector<DIM>::dot(s[i].v, dhat);
  vjpar = dhat*vector<DIM>::dot(s[j].v, dhat);
  viperp = s[i].v - vipar;
  vjperp = s[j].v - vjpar;

  s[i].v = vjpar + dhat*2.*growthrate + viperp;
  s[j].v = vipar - dhat*2.*growthrate + vjperp;

  // momentum exchange
  double xvelocity;   // exchanged velocity
  xvelocity = vector<DIM>::dot(s[i].v - s[j].v, dhat) - 4.*growthrate;
  xmomentum += M*xvelocity*dhatmagnitude;
}


//==============================================================
// Transfer, takes care of boundary events too
//=============================================================
void box::Transfer(event e)
{
  gtime = e.time;
  int i = e.i;
  int j = e.j;
  int k=0;           // dimension perpendicular to wall it crosses
 
  // update position and lutime (velocity doesn't change)
  s[i].x += s[i].v*(gtime-s[i].lutime);
  s[i].lutime = gtime;

  vector<DIM,int> celli;  // new cell for i
  celli = s[i].cell;  // this is not redundant
  
  // update cell
  if (j>N+DIM+1)  // right wall
    {
      k = j-N-DIM-2;
      celli[k] = s[i].cell[k] + 1;

      // if in right-most cell, translate x and cell
      if (s[i].cell[k] == ngrids - 1)
	{
	  s[i].x[k] -= SIZE;
	  celli[k] -= ngrids;
	}
    }
  else if (j<N+DIM+1)  // left wall
    {
      k = -j+N+DIM;
      celli[k] = s[i].cell[k] - 1;

      // if in left-most cell, translate x and cell
      if (s[i].cell[k] == 0)
	{
	  s[i].x[k] += SIZE;
	  celli[k] += ngrids;
	}
    }
  else
    std::cout << "error in Transfer" << std::endl;
      
  UpdateCell(i, celli); 
}


//==============================================================
// Updates cell of a sphere to time
//=============================================================
void box::UpdateCell(int i, vector<DIM,int>& celli)
{
  if(ngrids==1) return;
   
  if (celli == s[i].cell)
    std::cout << "error in update cell..shouldn't be the same" << std::endl;
  
  // delete i from cell array at cell

  else if (cells.get(s[i].cell) == i) 
    {
      if (binlist[i] == -1)
	cells.get(s[i].cell) = -1;
      else
	{
	  cells.get(s[i].cell) = binlist[i];
	  binlist[i] = -1;
	}
    }

  else if (cells.get(s[i].cell) == -1)
    {
      std::cout << "error " << i << " not in claimed cell UpdateCell" << std::endl;
      OutputCells();
    }

  else  // if no, find i in binlist
    {  
      int iterater = cells.get(s[i].cell);
      int pointer = iterater;
      while ((iterater != i)&&(iterater != -1))
	{
	  pointer = iterater;
	  iterater = binlist[iterater];
	}
      if (iterater == -1)  // got to end of list without finding i
	{
	  std::cout << "problem " << i << " wasn't in claimed, cell iterater = -1" << std::endl;
	  OutputCells();
	}
      else  // we found i!
	{
	  binlist[pointer] = binlist[i]; 
	  binlist[i] = -1;
	}	  
    } 

  // now add i to cell array at celli
  s[i].cell = celli;
  
  //first check to see if entry at celli
  if (cells.get(celli) == -1) //if yes, add i to cells gridfield
    cells.get(celli) = i;
  else  // if no, add i to right place in binlist
    {
      int iterater = cells.get(celli);  // now iterate through to end and add i
      int pointer = iterater;
      while (iterater != -1)  // find the end of the list
	{
	  pointer = iterater;
	  iterater = binlist[iterater];
	}
      binlist[pointer] = i;
      binlist[i] = -1; // redundant
    }
}


//==============================================================
// Output event heap...purely used for debugging
//==============================================================
void box::OutputEvents()
{
  h.print();
}


//==============================================================
// Output positions of spheres and their cells...purely used for debugging
//==============================================================
void box::OutputCells()
{
  for (int i=0; i<N; i++)
    std::cout << i << " " << s[i].x << " " << s[i].v << " " << s[i].cell << std::endl;
}


//==============================================================
// Update positions...purely for graphical display
//==============================================================
void box::TrackPositions()
{
  for (int i=0; i<N; i++)
    x[i] = s[i].x + s[i].v*(gtime-s[i].lutime);
}


//==============================================================
// Computes the total energy
//==============================================================
double box::Energy()
{
  double E=0;
  for (int i=0; i<N; i++)
    E += 0.5*M*s[i].v.norm_squared();

  return E/N;
}


//==============================================================
// Calculates the packing fraction
//==============================================================
double box::PackingFraction()
{
  double r_now = r + gtime*growthrate;
  double v = (pow(sqrt(PI)*r_now, DIM))/(exp(lgamma(1.+((double)(DIM))/2.)));
  return N*v/(pow(SIZE, DIM));
}


//==============================================================
// Calculates the optimal ngrids
//==============================================================
int box::Optimalngrids(double maxpf)
{
  double maxr;

  maxr = pow(exp(lgamma(1.+((double)(DIM))/2.))*maxpf/N, 1./DIM)/sqrt(PI);

  return (int)(1./(2.*maxr));
}


//==============================================================
// Processes n events
//==============================================================
void box::Process(int n)
{
  double deltat = gtime;
  for (int i=0; i<n; i++)
    {
      ProcessEvent();
    }
  pf = PackingFraction();   // packing fraction
  deltat = gtime - deltat;
  double oldenergy = energy;
  energy = Energy();        // kinetic energy

  energychange = ((oldenergy - energy)/oldenergy)*100; // percent change in energy

  if (deltat != 0.) 
    pressure = 1+xmomentum/(2.*energy*N*deltat);
 
  // reset to 0
  ncollisions = 0;
  ntransfers = 0;
  nchecks = 0;
  xmomentum = 0.;
  ncycles++;
}


//==============================================================
// Prints statistics for n events
//==============================================================
void box::PrintStatistics()
{
  std::cout << "packing fraction = " << pf << std::endl;
  std::cout << "gtime = " << gtime << std::endl; 
  std::cout << "total time = " << rtime+gtime << std::endl;
  std::cout << "kinetic energy = " << energy << std::endl;
  std::cout << "total # events = " << neventstot << std::endl;
  std::cout << "# events = " << ncollisions+ntransfers+nchecks << ", # collisions = " << ncollisions << ", # transfers = " << ntransfers << ", # checks =" << nchecks << std::endl;
  std::cout << "growthrate = " << growthrate << std::endl;
  std::cout << "reduced pressure = " << pressure << std::endl;
  std::cout << "-----------------" << std::endl;
}


//==============================================================
// Updates spheres to gtime, synchronizes, and can change growth rate
//==============================================================
void box::Synchronize(bool rescale)
{
  double vavg = sqrt(2.*M*energy);

  for (int i=0; i<N; i++)
    {
      s[i].x = s[i].x + s[i].v*(gtime-s[i].lutime);
      s[i].nextevent.time -= gtime;

      if (s[i].nextevent.time < 0.)
	std::cout << "error, event times negative after synchronization" << std::endl;
      if (rescale == true)   // give everyone checks
	{
	  s[i].nextevent = event(0., i, INF); 
	  s[i].v /= vavg;
	}
	  
      s[i].lutime = 0.;
    }
  r += gtime*growthrate;       // r defined at gtime = 0
  rtime += gtime;
  gtime = 0.;
  
  if (rescale == true)
    Process(N);
}


//==============================================================
// Run time
//==============================================================
void box::RunTime()
{
  time(&end);
  std::cout << "run time = " << difftime(end, start) << std::endl;
}


//==============================================================
// Write configuration
//==============================================================
void box::WriteConfiguration(const char* wconfigfile)
{
  if (gtime != 0.)   // synchronize spheres if not currently synchronized
    Synchronize(false);
      
  std::ofstream output(wconfigfile);

  // make header
  output << DIM << "\n";
  output << N << " " << 1 << "\n";
  output << N << "\n";
  output << std::setprecision(16) << 2.*r << "\n";

// FIX-ME  
#if (DIM == 2)
  output << 1 << " " << 0 << " " << 0 << " " << 1 << "\n";
#elif (DIM == 3)
  output << 1 << " " << 0 << " " << 0 << " " << 0 << " " << 1 << " " 
	 << 0 << " " << 0 << " " << 0 << " " << 1 << "\n";
#elif (DIM == 4)
  output << 1 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " "
	 << 1 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " 
	 << 1 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " "
	 << 1 << "\n";
#endif

  output << "T T T" << "\n";
  
  // remember to get 16 digits of accuracy
  for (int i=0; i<N; i++)  
    {
      for (int k=0; k<DIM; k++)
	output << std::setprecision(16) << s[i].x[k] << " ";
      output << "\n";
    }
      
  output.close();
}
