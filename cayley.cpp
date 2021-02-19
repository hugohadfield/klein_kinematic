#include <iostream>
#include <klein/klein.hpp>
#include "klein_ops.h"


kln::rotor cayley(kln::branch phi){
    kln::rotor phi2 = (phi*phi);
    return (1.0f + phi)*(1.0f + phi)/(1.0f - phi2.scalar());
}

kln::branch cayley_kinematic(kln::branch phi, kln::branch omega){
    kln::rotor omega_rot = as_rotor(omega);
    return as_branch(0.25f*(1.0f + phi)*omega_rot*(1.0f + (-1.0f*phi)));
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

    kln::rotor R = cayley(phi);
    std::cout << R.scalar() << " " << R.e23() << " " << R.e31() << " " << R.e12() << std::endl;    

    auto kin_output = cayley_kinematic(phi, omega);
    std::cout << kin_output.e23() << " " << kin_output.e31() << " " << kin_output.e12() << std::endl;    
    return 0;
}
