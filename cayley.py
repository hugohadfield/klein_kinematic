
from clifford.pga import *

def explicit_inverse(X):
    return (~X)*(((X*~X).value[0] - (X*~X)(4))/((X*~X).value[0]**2))


def cayley_so3(phi):
    """
    Implements the simplified so3 cayley map from bivectors to rotors
    """
    phi2 = (phi*phi);
    return (1.0 + phi)*(1.0 + phi)/(1.0 - phi2.value[0]);



def cayley_se3(phi):
    """
    Implements the simplified se3 cayley map from bivectors to rotors
    """
    phi2 = (phi*phi)
    denominator = 1.0 - phi2.value[0]
    return (1.0 + phi)*(1.0 + phi)*(( (1.0 - phi2.value[0]) + phi2(4))/(denominator*denominator));



def cayley_unsimplified(phi):
    """
    Implements the bidirectional cayley map in its unsimplified form
    """
    return (1 + phi)*(1 - phi).inv()


def cayley_kinematic(phi, omega):
    """
    This is the kinematic equation for the cayley map as found in
    Hadfield H., Lasenby J., Screw Theory in Geometric Algebra for Constrained Rigid Body Dynamics AACA (2021)
    """
    return 0.25*(1.0 + phi)*omega*(1.0 - phi)


def test_so3():
    a = 0.1
    b = 0.2
    c = 0.3

    d = 0.4
    e = 0.5
    f = 0.6

    phi = a*(e23) + b*(e3*e1) + c*(e12)
    omega = d*(e23) + e*(e3*e1) + f*(e12)

    R = cayley_so3(phi);
    print(R.value[0], R[2,3], R[3,1], R[1,2])

    Ralt = cayley_unsimplified(phi)
    print(Ralt.value[0], Ralt[2,3], Ralt[3,1], Ralt[1,2])

    kin_output = cayley_kinematic(phi, omega);
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

    R = cayley_se3(phi)
    print(R.value[0], R[2,3], R[3,1], R[1,2], R[0,1], R[0,2], R[0,3], R[0,1,2,3])

    Ralt = cayley_unsimplified(phi)
    print(Ralt.value[0], Ralt[2,3], Ralt[3,1], Ralt[1,2], Ralt[0,1], Ralt[0,2], Ralt[0,3], Ralt[0,1,2,3])

    kin_output = cayley_kinematic(phi, omega);
    print(kin_output[2,3], kin_output[3,1], kin_output[1,2], kin_output[0,1], kin_output[0,2], kin_output[0,3])   



if __name__ == "__main__":
    # print('so3')
    # test_so3()
    print('se3')
    test_se3()
