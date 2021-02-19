#include <iostream>
#include <klein/klein.hpp>
#include "klein_ops.h"


kln::rotor outer_exp(kln::branch phi){
    kln::rotor phi2 = (phi*phi);
    return (1.0f + phi)/std::sqrt(1.0f - phi2.scalar());
}

kln::branch outer_exp_kinematic(kln::branch phi, kln::branch omega){
    auto R = outer_exp(phi);
    auto omegaR = omega*R;
    return -0.5f*sqrtf(1.0f - (phi*phi).scalar())*(-as_branch(omegaR) + omegaR.scalar()*as_branch(R)/R.scalar());
}

int main(){    
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

    auto kin_output = outer_exp_kinematic(phi, omega);
    std::cout << kin_output.e23() << " " << kin_output.e31() << " " << kin_output.e12() << std::endl;    
    return 0;
}
