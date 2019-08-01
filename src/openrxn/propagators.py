import numpy as np

def Gillespie(processes,time_range,y0):
    """A propagator function that moves the state vector (y)
    forward in time.

    Inputs:

    processes :  a list of "reactions" with elements of the 
                 form: (rate, q_list, delta_list).  Where 
                 the rate constant is in 1/s, the q_list is a 
                 list of y indices to multiply together with the 
                 rate constant to determine the reactive flux, 
                 and delta_list is a list of tuples (idx, delta)
                 that describe how to update the state vector if
                 that reaction is chosen.
    
    time_range : a tuple (t_init, t_final) describing the range 
                 over which reactions should be processed

    y0 :         the initial value of the state vector

    Returns:

    y_final :    the final value of the state vector
    t_final :    the actual final time
    """

    t = time_range[0]

    y = y0.copy()
    n = len(processes)

    while t < time_range[1]:
        # build vector of reaction rates (r)
        r = np.zeros(n)
        for i,p in enumerate(processes):
            r[i] = p[0]
            for idx in p[1]:
                r[i] *= y[idx]

        # save 1/rsum since multiplication is faster
        oorsum = 1/r.sum()
        
        # choose a reaction to execute
        i = np.random.choice(range(n),p=r*oorsum)

        # update y
        for idx, d in processes[i][2]:
            y[idx] += d

        # update t
        a = np.random.random()
        t += -np.log(a)*oorsum

    return y, t

        


        
        
