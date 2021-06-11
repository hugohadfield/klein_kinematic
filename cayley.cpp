#include <iostream>
#include <klein/klein.hpp>
#include "klein_ops.h"


kln::rotor cayley(kln::branch phi){
    /*
    Implements the simplified so3 cayley map from bivectors to rotors
    */
    kln::rotor phi2 = (phi*phi);
    return (1.0f + phi)*(1.0f + phi)/(1.0f - phi2.scalar());
}


kln::motor cayley(kln::line phi){
    /*
    Implements the simplified se3 cayley map from bivectors to rotors
    */
    kln::motor phi2 = (phi*phi);
    auto denominator = 1.0f - phi2.scalar();
    auto phi2_4 = kln::motor(0,0,0,0,0,0,0,phi2.e0123());
    return (1.0f + phi)*(1.0f + phi)*(denominator + phi2_4)/(denominator*denominator);
}


kln::branch cayley_kinematic(kln::branch phi, kln::branch omega){
    /*
    This is the kinematic equation for the cayley map as found in
    Hadfield H., Lasenby J., Screw Theory in Geometric Algebra for Constrained Rigid Body Dynamics AACA (2021)
    */
    kln::rotor omega_rot = as_rotor(omega);
    return as_branch(0.25f*(1.0f + phi)*omega_rot*(1.0f + (-1.0f*phi)));
}


kln::motor cayley_kinematic(kln::line phi, kln::line omega){
    /*
    This is the kinematic equation for the cayley map as found in
    Hadfield H., Lasenby J., Screw Theory in Geometric Algebra for Constrained Rigid Body Dynamics AACA (2021)
    */
    kln::motor omega_rot = as_motor(omega);
    return (0.25f*(1.0f + phi)*omega_rot*(1.0f + (-1.0f*phi)));
}


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

    auto kin_output = cayley_kinematic(phi, omega);
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

    auto kin_output = cayley_kinematic(phi, omega);
    std::cout << kin_output.e23() << " " << kin_output.e31() << " " << kin_output.e12() 
            << " " << kin_output.e01() << " " << kin_output.e02() << " " << kin_output.e03() << std::endl;    
}


int main(){
    test_se3();
    return 0;
}