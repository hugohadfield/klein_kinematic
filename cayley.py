
from clifford.g3 import *


def cayley(phi):
    """
    Implements the simplified so3 cayley map from bivectors to rotors
    """
    phi2 = (phi*phi);
    return (1.0 + phi)*(1.0 + phi)/(1.0 - phi2.value[0]);


def cayley_unsimplified(phi):
    """
    Implements the bidirectional cayley map in its unsimplified form
    """
    return (1 + phi)/(1 - phi)


def cayley_kinematic(phi, omega):
    """
    This is the kinematic equation for the cayley map as found in
    Hadfield H., Lasenby J., Screw Theory in Geometric Algebra for Constrained Rigid Body Dynamics AACA (2021)
    """
    return 0.25*(1.0 + phi)*omega*(1.0 - phi)


if __name__ == "__main__":
    a = 0.1
    b = 0.2
    c = 0.3

    d = 0.4
    e = 0.5
    f = 0.6

    phi = a*(e23) + b*(e3*e1) + c*(e12)
    omega = d*(e23) + e*(e3*e1) + f*(e12)

    R = cayley(phi);
    print(R.value[0], R[2,3], R[3,1], R[1,2])

    Ralt = cayley_unsimplified(phi)
    print(Ralt.value[0], Ralt[2,3], Ralt[3,1], Ralt[1,2])

    kin_output = cayley_kinematic(phi, omega);
    print(kin_output[2,3], kin_output[3,1], kin_output[1,2])   
