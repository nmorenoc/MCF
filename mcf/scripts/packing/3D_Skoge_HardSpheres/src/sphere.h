#ifndef  SPHERE_H
#define  SPHERE_H


#include "vector.h"


class sphere {

 public:

  // constructor and destructor

  sphere();
  sphere(const sphere& s);
  sphere(int i_i, vector<DIM> x, vector<DIM, int> cell_i, double lutime_i);
  ~sphere();

 //variables
  
  int i;                          // sphere ID

  // impending event
  event nextevent;                // next event...can be collision or transfer
  event nextcollision;            // next collision if next event is transfer
  // maybe nextnext event
  
  // past information
  double lutime;                  // last update time
  vector<DIM, int> cell;          // cell that it belongs to
  vector<DIM, double> x;          // position
  vector<DIM, double> v;          // velocity
  // make sure efficent in memory

 

};

#endif 
