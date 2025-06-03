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

def eisermann(L):
    n = len(L.link_components)
    jones_L = L.jones_polynomial()
    jones_unlink = (-(q + q^(-1))^(n-1))
    if not jones_unlink.divides(jones_L):
        return false
    quotient = jones_L / jones_unlink
    components = L.link_components
    determinant_product = 1
    for component in components:
        K = L.sublink(component)
        if (len(K.link_components) != 0):
            determinant_product *= K.determinant()
    quotientI = quotient(I)
    if (quotientI.imag_part() == 0 and (quotientI.real() % 32 == determinant_product % 32 or int(quotientI.real()) % 32 == -determinant_product % 32)):
        return True
    return False
