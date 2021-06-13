#include <iostream>
#include <klein/klein.hpp>
#include "klein_ops.h"



kln::branch outer_log(kln::rotor  R){
    /*
    For a given rotor this returns the bivector that when outer exponeniated gives the rotor
    */
    return as_branch(R)/R.scalar();
}


kln::line outer_log(kln::motor  R){
    /*
    For a given motor this returns the bivector that when outer exponeniated gives the rotor
    */
    return as_line(R)/R.scalar();
}


kln::rotor outer_exp(kln::branch phi){
    /*
    Implements the so3 outer exponential from bivectors to rotors
    */
    kln::rotor phi2 = (phi*phi);
    return (1.0f + phi)/std::sqrt(1.0f - phi2.scalar());
}


kln::motor outer_exp(kln::line phi){
    /*
    Implements the se3 outer exponential from bivectors to rotors
    */
    kln::motor phi2 = (phi*phi);
    kln::motor phi2_4 = kln::motor(0,0,0,0,0,0,0,phi2.e0123());
    return (1.0f + phi + 0.5f*phi2_4)/std::sqrt(1.0f - phi2.scalar());
}


kln::branch outer_exp_kinematic(kln::branch phi, kln::branch omega){
    /*
    This is the kinematic equation for the outer exponential map as found in
    Hadfield H., Lasenby J., Screw Theory in Geometric Algebra for Constrained Rigid Body Dynamics AACA (2021)
    */
    kln::rotor R = outer_exp(phi);
    kln::rotor omegaR = omega*R;
    return -0.5f*sqrtf(1.0f - (phi*phi).scalar())*(-as_branch(omegaR) + omegaR.scalar()*as_branch(R)/R.scalar());
}


kln::line outer_exp_kinematic(kln::line phi, kln::line omega){
    /*
    This is the kinematic equation for the outer exponential map as found in
    Hadfield H., Lasenby J., Screw Theory in Geometric Algebra for Constrained Rigid Body Dynamics AACA (2021)
    */
    kln::motor R = outer_exp(phi);
    kln::motor omegaR = omega*R;
    return -0.5f*sqrtf(1.0f - (phi*phi).scalar())*(-as_line(omegaR) + omegaR.scalar()*as_line(R)/R.scalar());
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

    kln::rotor R = outer_exp(phi);
    std::cout << R.scalar() << " " << R.e23() << " " << R.e31() << " " << R.e12() << std::endl;    

    kln::branch kin_output = outer_exp_kinematic(phi, omega);
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

    kln::motor R = outer_exp(phi);
    std::cout << R.scalar() << " " << R.e23() << " " << R.e31() << " " << R.e12() 
        << " " << R.e01() << " " << R.e02() << " " << R.e03() << " " << R.e0123() << std::endl;    

    kln::line kin_output = outer_exp_kinematic(phi, omega);
    std::cout << kin_output.e23() << " " << kin_output.e31() << " " << kin_output.e12() 
            << " " << kin_output.e01() << " " << kin_output.e02() << " " << kin_output.e03() << std::endl;  
      
}


int main(){
    test_se3();
    return 0;
}
