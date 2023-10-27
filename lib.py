from Compiler.types import *
from Compiler.util import if_else
from Compiler.library import *
from Compiler.instructions import *

def multvec(spA, ijA, x, n):
    # fast multiplication of sparse A with vector x
    b = sfix.Array(n)

    @for_range_opt(n)
    def _(i):
        b[i] = spA[i]*x[i]

    @for_range(n)
    def _(i):
        @for_range(ijA[i]-1, ijA[i+1] - 1)  # index arithmetic
        def _(k):
            b[i] += spA[k] * x[ijA[k] - 1]

    return b

def evalF(x, spG, ijG, spB, ijB, pK, qK, n, slack_voltage):
    uKreal = sfix.Array(n)
    uKreal[0] = slack_voltage
    uKreal[1:] = x[0:n-1]
    uKimag = sfix.Array(n)
    uKimag[0] = 0
    uKimag[1:] = x[n-1:]

    # Compute complex-conjugate node currents
    b = sfix.Array(n)
    b.assign_all(0)
    c = sfix.Array(n)
    c.assign_all(0)
    d = sfix.Array(n)
    d.assign_all(0)
    e = sfix.Array(n)
    e.assign_all(0)

    b = multvec(spG, ijG, uKreal, n)
    c = multvec(spB, ijB, uKimag, n)
    d = multvec(spB, ijB, uKreal, n)
    e = multvec(spG, ijG, uKimag, n)

    iKreal_conj = sfix.Array(n)
    iKreal_conj.assign_all(0)
    iKimag_conj = sfix.Array(n)
    iKimag_conj.assign_all(0)
    iKreal_conj[:] += b[:] - c[:]
    iKimag_conj[:] += (-1 * d[:]) - e[:]

    # Compute node powers while ingnoring slack node 1
    pN = sfix.Array(n-1)
    pN.assign_all(0)
    qN = sfix.Array(n-1)
    qN.assign_all(0)

    @for_range_opt(n-1)
    def _(i):
        pN[i] = (uKreal[i+1] * iKreal_conj[i+1]) - \
            (uKimag[i+1] * iKimag_conj[i+1])
        qN[i] = (uKreal[i+1] * iKimag_conj[i+1]) + \
            (uKimag[i+1] * iKreal_conj[i+1])

    fval = sfix.Array((n-1) * 2)
    fval.assign_part_vector(pN[:] - pK[:])
    fval.assign_part_vector(qN[:] - qK[:], base=n-1)


    return (fval, iKreal_conj, iKimag_conj)

def compute_alpha(x,dx,Fx,spG,ijG,spB,ijB,pK,qK,num_nodes,slack_voltage):
    # Compute an (approximately) optimal step length for Newton's method...
    # using a secant method.

    n = (num_nodes - 1) * 2
    alpha_old = MemValue(sfix(1)) # MemValues so variables can be used in while-loop
    alpha = MemValue(sfix(0.95))
    dalpha = MemValue(sfix(-0.05))
    xnew = sfix.Array(n)
    xnew.assign_all(0)
    Fj = sfix.Array(n)
    Fj.assign_all(0)
    Fjm1 = sfix.Array(n)
    Fjm1.assign_all(0)
    P = sfix(0)
    denum = sfix(0)
    y = sint(0)
    # Below steps are necessary to replace if-condition in while loop
    xnew[:] = x[:] + dx[:]
    Fjm1 = evalF(xnew,spG,ijG,spB,ijB,pK,qK,num_nodes,slack_voltage)[0]

    # To avoid overflow
    Fjm1[:] *= 10**-5

    dx_alpha = sfix.Array(n)
    dx_alpha.assign(dx)
    dx_alpha[:] *= alpha
    xnew[:] = x[:] + dx_alpha[:]
    Fj = evalF(xnew,spG,ijG,spB,ijB,pK,qK,num_nodes,slack_voltage)[0]
    
    # To avoid overflow
    Fj[:] *= 10**-5
    Fx[:] *= 10**-5
    
    Fj_m = sfix.Matrix(n,1)
    Fj_m.assign(Fj)
    P = Fj_m.transpose().dot(Fj[:] - Fx[:])[0]
    Fjm1_m = sfix.Matrix(n,1)
    Fjm1_m.assign(Fjm1)
    @for_range_opt(n)
    def _(i):
        Fjm1_m[i][0] *= alpha
    denum = Fjm1_m.transpose().dot(Fjm1[:] - Fx[:])[0] - P
    dalpha.write(dalpha * P/denum)
    alpha_old.write(alpha)
    alpha.write(alpha + dalpha)
    dalpha_check = dalpha * 10
    y = if_else(dalpha_check**2 > 0.1, 1, 0) 
    @while_do(lambda x: x != 0, y.reveal())
    def _(i):
        # Fjm1 = Fj
        Fj = evalF(xnew, spG, ijG, spB, ijB, pK,
                   qK, num_nodes, slack_voltage)[0] # we have to compute this again for some reason
        
        # To avoid overflow
        Fj[:] *= 10**-5

        Fjm1.assign(Fj)  # reuse value from before
        dx_alpha.assign_all(0)
        dx_alpha.assign(dx)
        dx_alpha[:] *= alpha
        xnew[:] = x[:] + dx_alpha[:]
        Fj=evalF(xnew,spG,ijG,spB,ijB,pK,qK,num_nodes,slack_voltage)[0]
        
        # To avoid overflow
        Fj[:] *= 10**-5

        Fj_m.assign_all(0)
        Fj_m.assign(Fj)
        @for_range_opt(n)
        def _(j):
            Fj_m[j][0] *= alpha_old
        P.update(Fj_m.transpose().dot(Fj[:] - Fx[:])[0])
        @for_range_opt(n)
        def _(j):
            Fjm1_m[j][0] = alpha * Fjm1[j]
        denum.update(Fjm1_m.transpose().dot(Fjm1[:] - Fx[:])[0] - P)
        dalpha.write(dalpha * P/denum)

        #update
        alpha_old.write(alpha)
        alpha.write(alpha + dalpha)
        dalpha_check = dalpha * 10
        i = if_else(dalpha_check**2 > 0.1, 1, 0) # public
        return i.reveal()

    return alpha
