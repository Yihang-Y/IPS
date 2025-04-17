import numpy as np

def A_step( qp , h ):
    q,p = qp
    
    q = q + h*p
    
    return [q,p]

def B_step( qp , h, force ):
    q,p = qp
    
    F = force(q)
    
    p = p + h*F
    
    return [q,p]

def O_step( qp , h,gamma, beta ):
    q,p = qp
    
    alpha = np.exp(-h*gamma)
    
    R = np.random.randn( q.size ).reshape( q.shape)
    p = alpha*p + np.sqrt(1/beta) * np.sqrt(1 - alpha**2) * R
    
    return [q,p]

def ld_ABO(state,h,params:dict,force):
    # The algorithm "ABO" does A then B then O 
    if not 'gamma' in params or not 'beta' in params:
        raise ValueError("Params dictionary must contain 'gamma' and 'beta'")
    
    q,p = state
    qp = np.copy([q,p])  #this just translates the separate q and p vectors 
                #into a single vector composed from the pair.
    
    qp = A_step(qp , h )
    qp = B_step(qp, h, force)
    qp = O_step(qp, h, params['gamma'], params['beta'])
    
    return qp