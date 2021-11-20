function f = trans_eqn(E,e,M)
        f = E-e*sin(E)-mod(M,2*pi);