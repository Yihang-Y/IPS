import numpy as np

def A_step( qpxi , h ):
    q,p,xi = qpxi
    
    q = q + h*p
    
    return [q,p,xi]

def B_step(qpxi, h, force: callable):
    q,p,xi = qpxi

    F = force(q)
    p = p + h*F
    return [q,p,xi]

def C_step(qpxi, h):
    q,p,xi = qpxi

    p = p*np.exp(-xi*h)
    return [q,p,xi]

def D_step(qpxi, eta, T, h):
    q,p,xi = qpxi

    xi = xi + eta*(p**2 - T)*h
    return [q,p,xi]

def NHD_ABCDBA(state,h, params:dict,force):
    # this is probably expensive to run at every step but oh well
    if not 'eta' in params or not 'T' in params:
        raise ValueError("Params dictionary must contain 'eta' and 'T'")
    
    qpxi = np.copy(state)
    #print(qpxi)
    eta = params['eta']
    T = params['T']
    
    qpxi = A_step(qpxi , h/2)
    qpxi = B_step(qpxi, h/2, force)
    qpxi = C_step(qpxi, h)
    qpxi = D_step(qpxi, eta, T, h)
    qpxi = B_step(qpxi, h/2, force)
    qpxi = A_step(qpxi , h/2)
    
    return qpxi

def O_step(qpxi, gamma, T, h):
    q,p,xi = qpxi
    xi = xi*np.exp(-gamma*h)
    xi = xi + np.sqrt(T*(1-np.exp(-2*gamma*h)))*np.random.normal(0,1)
    return [q,p,xi]

def NHLD_onestep(state, h, params:dict, force):
    if not 'eta' in params or not 'T' in params or not 'gamma' in params:
        raise ValueError("Params dictionary must contain 'eta', 'gamma', and 'T'")
    
    qpxi = np.copy(state)
    eta = params['eta']
    T = params['T']
    gamma = params['gamma']

    qpxi = A_step(qpxi , h/2)
    qpxi = B_step(qpxi, h/2, force)
    qpxi = D_step(qpxi, eta, T, h/2)
    qpxi = C_step(qpxi, h/2)
    qpxi = O_step(qpxi, gamma, T, h)
    qpxi = C_step(qpxi, h/2)
    qpxi = D_step(qpxi, eta, T, h/2)
    qpxi = B_step(qpxi, h/2, force)
    qpxi = A_step(qpxi , h/2)

    return qpxi
    