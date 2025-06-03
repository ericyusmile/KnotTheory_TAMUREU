import csv
import snappy
from snappy import *

def sig_zero(L):
    return (L.signature() == 0)

def linking_matrix_zero(L):
    output = True
    M = L.linking_matrix()
    for i in M:
        for j in i:
            if (j != 0):
                output = False
    return output

#placeholder
def fox_milnor(L):
    return ""

#helper function for eisermann
def signed_determinant(K):
    V = K.seifert_matrix()
    return (-I * (V + V.transpose())).determinant()

def eisermann(L):
    n = len(L.link_components)
    jones_L = L.jones_polynomial()
    q = jones_L.parent().gen()
    jones_unlink = (q + q^(-1))^(n-1)
    if not jones_unlink.divides(jones_L):
        return false
    quotient = jones_L / jones_unlink
    determinant_product = 1
    for component in L.link_components:
        K = L.sublink(component)
        if (len(K.link_components) != 0):
            determinant_product *= real(signed_determinant(K))
    quotientI = quotient(I)
    if (quotientI.imag_part() == 0 and quotientI.real() % 32 == determinant_product % 32):
        return True
    return False


def test(L):
    a = sig_zero(L)
    b = linking_matrix_zero(L)
    c = fox_milnor(L)
    d = eisermann(L)
    if (a and b):
        print("hit!")
        with open("output.csv", mode="a", newline="") as file:
            writer = csv.writer(file)
            writer.writerow([
                L.PD_code(), 
                len(L.link_components),
                len(L.crossings), 
                a, 
                b, 
                c, 
                d
            ])
