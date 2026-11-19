import numpy as np

def Gillespie(processes,update_list,time_range,y0):
    """A propagator function that moves the state vector (y)
    forward in time.

    Inputs:

    processes :  a list of "reactions" with elements of the 
                 form: (rate, q_list, delta_list).  Where 
                 the rate constant is in 1/s, the q_list is a 
                 list of tuples (idx, number) that describe how
                 many of each species are involved in the reaction
                 and delta_list is a list of tuples (idx, delta)
                 that describe how to update the state vector if
                 that reaction is chosen.

    update_list: a list (same size as y0) with the elements 
                 being lists of processes that a given quantity 
                 is involved in.  It is used to make this propagator 
                 more efficient, only updating the processes that 
                 have changed after each reaction.
    
    time_range : a tuple (t_init, t_final) describing the range 
                 over which reactions should be processed

    y0 :         the initial value of the state vector

    Returns:

    y_final :    the final value of the state vector
    t_final :    the final time (time_range[1]; the next event
                 would occur after it)
    """

    t = time_range[0]

    y = y0.copy()
    n = len(processes)

    # build vector of reaction rates (r)
    r = np.zeros(n)
    for i,p in enumerate(processes):
        r[i] = p[0]
        for idx, num in p[1]:
            # for a two-body reaction:  r = k*y*(y-1)
            for j in range(num):
                r[i] *= max(0,y[idx]-j)

    while t < time_range[1]:
        # keep track of processes that will need to be updated
        to_update = []

        # nothing can fire: the state is frozen until the end
        rsum = r.sum()
        if rsum <= 0:
            t = time_range[1]
            break

        # save 1/rsum since multiplication is faster
        oorsum = 1/rsum

        # draw the waiting time first; stop at the end if it passes it
        a = np.random.random()
        t += -np.log(a)*oorsum
        if t >= time_range[1]:
            t = time_range[1]
            break

        # choose a reaction to execute
        i = np.random.choice(range(n),p=r*oorsum)

        # update y
        for idx, d in processes[i][2]:
            y[idx] += d
            to_update += update_list[idx]

        # update only the necessary r values
        for i in set(to_update):
            r[i] = processes[i][0]
            for idx, num in processes[i][1]:
                # for a two-body reaction:  r = k*y*(y-1)
                for j in range(num):
                    r[i] *= max(0,y[idx]-j)

    return y, t

        


        
        
