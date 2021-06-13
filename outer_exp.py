
import numpy as np
from clifford.pga import *


def outer_log(R):
    """
    For a given rotor this returns the bivector that when outer exponeniated gives the rotor
    """
    return R(2)/R.value[0]


def outer_exp_so3(phi):
    """
    Implements the so3 outer exponential from bivectors to rotors
    """
    phi2 = (phi*phi);
    return (1.0 + phi)/np.sqrt(1.0 - phi2.value[0]);


def outer_exp_se3(phi):
    """
    Implements the se3 outer exponential from bivectors to rotors
    """
    phi2 = (phi*phi);
    return (1.0 + phi + 0.5*phi2(4))/np.sqrt(1.0 - phi2.value[0]);


def outer_exp_kinematic_so3(phi, omega):
    """
    This is the kinematic equation for the outer exponential map as found in
    Hadfield H., Lasenby J., Screw Theory in Geometric Algebra for Constrained Rigid Body Dynamics AACA (2021)
    """
    R = outer_exp_so3(phi);
    omegaR = omega*R;
    return -0.5*np.sqrt(1 - (phi*phi).value[0])*(-omegaR(2) + omegaR.value[0]*R(2)/R.value[0]);


def outer_exp_kinematic_se3(phi, omega):
    """
    This is the kinematic equation for the outer exponential map as found in
    Hadfield H., Lasenby J., Screw Theory in Geometric Algebra for Constrained Rigid Body Dynamics AACA (2021)
    """
    R = outer_exp_se3(phi);
    omegaR = omega*R;
    return -0.5*np.sqrt(1 - (phi*phi).value[0])*(-omegaR(2) + omegaR.value[0]*R(2)/R.value[0]);


def test_so3():
    a = 0.1
    b = 0.2
    c = 0.3

    d = 0.4
    e = 0.5
    f = 0.6

    phi = a*(e23) + b*(e3*e1) + c*(e12)
    omega = d*(e23) + e*(e3*e1) + f*(e12)

    R = outer_exp_so3(phi);
    print(R.value[0], R[2,3], R[3,1], R[1,2])

    kin_output = outer_exp_kinematic_so3(phi, omega);
    print(kin_output[2,3], kin_output[3,1], kin_output[1,2])   



def test_se3():
    a = 0.1
    b = 0.2
    c = 0.3
    d = 0.4
    e = 0.5
    f = 0.6

    g = 0.1
    h = 0.2
    i = 0.3
    j = -0.4
    k = -0.5
    l = -0.6

    phi = a*(e23) + b*(e3*e1) + c*(e12) + g*e01 + h*e02 + i*e03
    omega = d*(e23) + e*(e3*e1) + f*(e12) + j*e01 + k*e02 + l*e03

    R = outer_exp_se3(phi)
    print(R.value[0], R[2,3], R[3,1], R[1,2], R[0,1], R[0,2], R[0,3], R[0,1,2,3])

    kin_output = outer_exp_kinematic_se3(phi, omega);
    print(kin_output[2,3], kin_output[3,1], kin_output[1,2], kin_output[0,1], kin_output[0,2], kin_output[0,3])   


def test_outer_log():
    for i in range(1000):
        phi = layout.randomMV()(2)
        R = outer_exp_se3(phi)
        phi_log = outer_log(R)
        np.testing.assert_allclose(phi.value, phi_log.value)


if __name__ == "__main__":
    test_se3()
