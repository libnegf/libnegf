#include <iostream>
#include <cmath>

int main()
{

  double E, kT;
  int i;
  
  kT = 0.026;
  std::cout.precision(16);
  std::cout.setf( std::ios::scientific, std::ios::floatfield ); 

  for(i=-100;i<100;i++)
  {
     E = i*0.1;
     std::cout << E<< "  "<<1.0/(exp(E/kT)+1.0)<< std::endl;
  }

}
