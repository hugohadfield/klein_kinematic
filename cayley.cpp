#include <iostream>
#include <klein/klein.hpp>
#include "cayley.h"


void test_so3(){    
    float a = 0.1f;
    float b = 0.2f;
    float c = 0.3f;
    float d = 0.4f;
    float e = 0.5f;
    float f = 0.6f;

    kln::branch phi = kln::branch(a, b, c);
    kln::branch omega = kln::branch(d, e, f); 

    kln::rotor R = cayley(phi);
    std::cout << R.scalar() << " " << R.e23() << " " << R.e31() << " " << R.e12() << std::endl;    

    kln::branch kin_output = cayley_kinematic(phi, omega);
    std::cout << kin_output.e23() << " " << kin_output.e31() << " " << kin_output.e12() << std::endl;    
}


void test_se3(){    
    float a = 0.1f;
    float b = 0.2f;
    float c = 0.3f;
    float d = 0.4f;
    float e = 0.5f;
    float f = 0.6f;

    float g = 0.1f;
    float h = 0.2f;
    float i = 0.3f;
    float j = -0.4f;
    float k = -0.5f;
    float l = -0.6f;

    kln::line phi = kln::line(g, h, i, a, b, c);
    kln::line omega = kln::line(j, k, l, d, e, f); 

    kln::motor R = cayley(phi);
    std::cout << R.scalar() << " " << R.e23() << " " << R.e31() << " " << R.e12() 
        << " " << R.e01() << " " << R.e02() << " " << R.e03() << " " << R.e0123() << std::endl;    

    kln::line kin_output = cayley_kinematic(phi, omega);
    std::cout << kin_output.e23() << " " << kin_output.e31() << " " << kin_output.e12() 
            << " " << kin_output.e01() << " " << kin_output.e02() << " " << kin_output.e03() << std::endl;  
  
}


int main(){
    test_se3();
    return 0;
}