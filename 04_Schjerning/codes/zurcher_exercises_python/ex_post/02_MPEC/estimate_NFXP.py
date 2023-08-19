#  NFXP class for structural estimation of discrete choice models.

import numpy as np
import scipy.optimize as optimize


ev = np.zeros(1) # Global variable

def estimate(model,solver,data,theta0=[0,0],twostep=0):
    """" Estimate model using NFXP"""
    global ev
    ev = np.zeros((model.n)) 
    
    samplesize = data.shape[0]
    
    # STEP 1: Find p non-parametrically
    tabulate = data.dx1.value_counts() # Count number of observations for each dx1
    p = [tabulate[i]/sum(tabulate) if i < len(tabulate) else 0 for i in range(len(model.p))]

    # STEP 2: Estimate structual parameters
    model.p[:] = p # Use first step estimates as starting values for p
    
    # Estimate RC and C
    pnames = ['RC','c']
    
    # Call BHHH optimizer
    res = optimize.minimize(ll,theta0,args = (model, solver, data, pnames), method = 'trust-ncg',jac = grad, hess = hes, tol=1e-8)
    # Update parameters
    model = updatepar(model,pnames,res.x)
    
    # Estimate RC, c and p
    if twostep == 0:
        pnames = ['RC','c','p']
        theta0 = [model.RC, model.c] + model.p.tolist() # Starting values
        # Call BHHH optimizer
        res = optimize.minimize(ll,theta0, args = (model,solver,data, pnames), method = 'trust-ncg',jac = grad, hess = hes, tol = 1e-8)

        # Update parameters
        model=updatepar(model,pnames,res.x)

    # Converged: "trust-ncg tends to be very conservative about convergence, and will often return status 2 even when the solution is good."
    converged   =   (res.status == 2 or res.status ==0)

    # Compute Variance-Covaiance matrix
    h = hes(res.x, model, solver,data, pnames) # Hessian
    Avar = np.linalg.inv(h*samplesize) # Variance-Covariance matrix from information matrix equality

    theta_hat = res.x # unpack estimates
    
    return model, res, pnames, theta_hat, Avar, converged

def ll(theta, model, solver,data, pnames, out=1, no_guess = True): # out=1 solve optimization
    """ Compute log-likelihood function """
    global ev # Use global variable to store value function to use as starting value for next iteration
    
    #Unpack and convert to numpy array
    x = np.array(data.x - 1) # x is the index of the observed state: We subtract 1 because python starts counting at 0
    d = np.array(data.d) # d is the observed decision
    dx1 = np.array(data.dx1) # dx1 is observed change in x 

    if no_guess == False:
        ev = np.zeros((model.n)) # Use previous value function as starting value
    
    # Update values
    model=updatepar(model,pnames,theta)
    model.p = np.abs(model.p)    # helps BHHH which is run as unconstrained optimization
    model.create_grid() # Update grid
    ev0 = ev # Use previous value function as starting value

    # Solve the model
    ev, pk, dev = solver.poly(model.bellman, V0=ev0 ,beta=model.beta, output=3)

    # Evaluate likelihood function
    lik_pr = pk[x] # Get probability of keeping given observed state    
    choice_prob = lik_pr * (1 - d) + (1-lik_pr) * d # get probability of making observed choice
    log_lik = np.log(choice_prob)  # Compute log-likelihood-contributions
    
    # add on log like for mileage process
    if theta.size>2: # theta > 2 if there are parameters for p
        p = np.append(model.p,1-np.sum(model.p)) # Add residual probability to p
        if any(p<=0):
            log_lik -= 100000*p[dx1] # Penalize if p is negative
        else:
            log_lik += np.log(p[dx1]) # Add log-likelihood contribution
        
    else:
        p = np.nan


    if out == 1:
        # Objective function (negative mean log likleihood)
        return np.mean(-log_lik)

    return model,lik_pr, pk, ev, dev, d,x,dx1

def score(theta, model, solver, data, pnames):
    """ Compute score function """
    global ev # Use global variable to store value function to use as starting value for next iteration
    model,lik_pr, pk, ev, dev, d,x,dx1 = ll(theta, model, solver, data, pnames,9) # Evaluate likelihood function
    F = np.eye(model.n)-dev # Get frechet derivative     
    N = data.x.size # Number of observations
    dc = 0.001*model.grid # Get derivative of cost function in utility wrt c
    pk = pk.reshape((model.n,1))  # Reshape to get correct shape for matrix multiplication



    ##### COMPUTE THE SCORE #######
    # Check if there are parameters for p
    if theta.size>2:
        n_p = len(model.p)
    else:
        n_p = 0


    # STEP 1: compute derivative of contraction operator wrt. parameters

    ## Derivative of utility function wrt. parameters
    dutil_dtheta=np.zeros((model.n, 2 + n_p, 2)) # shape is (gridsize, number of parameters, number of choices in utility function)
    dutil_dtheta[:,0, 0] = 0 # derivative of keeping wrt RC
    dutil_dtheta[:,0, 1] = -1 # derivative of replacing wrt RC
    dutil_dtheta[:,1, 0] = -dc # derivative of keeping wrt c
    dutil_dtheta[:,1, 1] = -dc[0] # derivative of replacing wrt c

    # Derivative of contraction operator wrt. utility parameters
    dbellman_dtheta=np.zeros((model.n, 2 + n_p)) # shape is (gridsize, number of parameters)
    dbellman_dtheta[:,:] =  (pk * dutil_dtheta[:, :, 0] + (1 - pk) * dutil_dtheta[:, :, 1])


    # Derivative of contraction operator wrt. p
    if theta.size>2:        
        vk = -model.cost + model.beta * model.P1 @ ev # Value of keeping
        vr = -model.RC-model.cost[0]+model.beta *model.P2 @ ev # Value of replacing
        vmax = np.maximum(vk,vr) # Get maximum value
        dbellman_dpi = vmax+np.log(np.exp(vk-vmax)+np.exp(vr-vmax)) #Re-centered log-sum: Value functin 
        for i_p in range(n_p): # loop over p
            part1 = dbellman_dpi[i_p:-1]
            part2 = np.hstack((dbellman_dpi[n_p:model.n], np.tile(dbellman_dpi[-1],(n_p-i_p-1))))
            dbellman_dtheta[0:model.n-i_p-1,2+i_p] =part1-part2
        
        invp=np.exp(-np.log(model.p))
        invp = np.vstack((np.diag(invp[0:n_p]),-np.ones((1,n_p))*invp[n_p-1]))
      
    # STEP 2: Compute derivative of Fixed point wrt. parameter
    dev_dtheta = np.linalg.solve(F,dbellman_dtheta)

    # STEP 3: Compute derivative of log-likelihood wrt. parameters


    score = np.zeros((N, 2 + n_p)) # Initialize score function
    for d_loop in range(2): # Loop over decisions (keep=0, replace=1)
        # Get transition matrix
        if d_loop == 0:
            P = model.P1
        else:
            P = model.P2
        dv = dutil_dtheta[:, :, d_loop] + model.beta * P @ dev_dtheta  # derivative of choice-specific value function wrt. parameters
        choice_prob = lik_pr * (1 - d_loop) + (1-lik_pr) * d_loop # get probability of choice in loop
        score += ((d == d_loop) -  choice_prob ).reshape(-1,1) * dv[x, :] # Add derivative of log-likelihood wrt. parameters

    # Add derivative of log-likelihood from mileage process wrt. p
    if theta.size>2:
        for i_p in range(n_p): 
            score[:,2+i_p] = score[:,2+i_p]+invp[dx1,i_p]

    return score

def grad(theta, model, solver,data, pnames):
    """ Compute gradient of log-likelihood function """
    s = score(theta, model, solver, data,pnames)
    return -np.mean(s,0)

def hes(theta, model, solver,data, pnames):
    """ Compute Hessian of log-likelihood function as outer product of scores"""
    s= score(theta, model, solver, data, pnames)

    return s.T@s/data.shape[0]

def updatepar(par,parnames, parvals):
    """ Update parameters """
    for i,parname in enumerate(parnames):
        # First two parameters are scalars
        if i<2:
            parval = parvals[i]
            setattr(par,parname,parval)
        else: # Remaining parameters are lists
            list_val = [None]*(parvals.size-2) 
            for j,parval in enumerate(parvals[2:]):
                list_val[j]=parval
            setattr(par,parname,list_val)
    return par
