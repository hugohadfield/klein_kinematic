
import numpy as np
from clifford.g3 import *


def outer_exp(phi):
    phi2 = (phi*phi);
    return (1.0 + phi)/np.sqrt(1.0 - phi2.value[0]);


def outer_exp_kinematic(phi, omega):
    R = outer_exp(phi);
    omegaR = omega*R;
    return -0.5*np.sqrt(1 - (phi*phi).value[0])*(-omegaR(2) + omegaR.value[0]*R(2)/R.value[0]);


if __name__ == "__main__":
    a = 0.1
    b = 0.2
    c = 0.3

    d = 0.4
    e = 0.5
    f = 0.6

    phi = a*(e23) + b*(e3*e1) + c*(e12)
    omega = d*(e23) + e*(e3*e1) + f*(e12)

    R = outer_exp(phi);
    print(R.value[0], R[2,3], R[3,1], R[1,2])

    kin_output = outer_exp_kinematic(phi, omega);
    print(kin_output[2,3], kin_output[3,1], kin_output[1,2])   
