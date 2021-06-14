#pragma once

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
    float denominator = 1.0f - phi2.scalar();
    kln::motor phi2_4 = kln::motor(0,0,0,0,0,0,0,phi2.e0123());
    return (1.0f + phi)*(1.0f + phi)*(denominator + phi2_4)/(denominator*denominator);
}


kln::branch cayley_kinematic(kln::branch phi, kln::branch omega){
    /*
    This is the kinematic equation for the cayley map as found in
    Hadfield H., Lasenby J., Screw Theory in Geometric Algebra for Constrained Rigid Body Dynamics AACA (2021)
    */
    kln::rotor omega_rot = as_rotor(omega);
    return 0.25f*as_branch((1.0f + phi)*omega_rot*(1.0f + -phi));
}


kln::line cayley_kinematic(kln::line phi, kln::line omega){
    /*
    This is the kinematic equation for the cayley map as found in
    Hadfield H., Lasenby J., Screw Theory in Geometric Algebra for Constrained Rigid Body Dynamics AACA (2021)
    */
    kln::motor omega_rot = as_motor(omega);
    return 0.25f*as_line((1.0f + phi)*omega_rot*(1.0f + -phi));
}
