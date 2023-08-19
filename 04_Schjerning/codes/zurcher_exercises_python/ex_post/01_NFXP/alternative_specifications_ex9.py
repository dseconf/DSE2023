#  NFXP class for structural estimation of discrete choice models.

import numpy as np
import scipy.sparse as sparse
import scipy.linalg as linalg
from scipy.sparse.linalg import spsolve
import scipy.optimize as optimize

import Solve_NFXP as Solve
import model_zucher as zucher


ev = np.zeros(1) # Global variable


def estimate(model,solver,data,ev_nul,theta0=[0,0],twostep=0):
    global ev
    ev = np.zeros((model.n))
    
    samplesize = data.shape[0]
    # STEP 1: Find p 
    tabulate = data.dx1.value_counts()
    p = [tabulate[i]/sum(tabulate) if i < len(tabulate) else 0 for i in range(len(model.p))]

    # STEP 2: Estimate structual parameters
    model.p[:] = p # Use first step estimates as starting values for p
    
    # Estimate RC and C
    pnames = ['RC','c']
    
    res = optimize.minimize(ll,theta0,args = (model, solver, data, pnames,ev_nul), method = 'trust-ncg',jac = grad, hess = hes, tol=1e-8)
    model=updatepar(model,pnames,res.x)
    
    # Estimate RC, c and p
    if twostep == 0:
        pnames = ['RC','c','p']
        theta0 = [model.RC, model.c] + model.p.tolist()
        res = optimize.minimize(ll,theta0, args = (model,solver,data, pnames,ev_nul), method = 'trust-ncg',jac = grad, hess = hes, tol = 1e-8)

        model=updatepar(model,pnames,res.x)


    theta_hat = res.x
    
    return res,theta_hat

def ll(theta, model, solver,data, pnames, ev_nul, out=1): # out=1 solve optimization, out =2 find hessian or the gradient
    global ev
    
    #Update model parameters
    x = np.array(data.x - 1) # x is the index of the observed state: We subtract 1 because python starts counting at 0
    d = np.array(data.d) # d is the observed decision
    dx1 = np.array(data.dx1) # dx1 is observed change in x 

    model=updatepar(model,pnames,theta)
    model.p = np.abs(model.p)    # helps BHHH which is run as unconstrained optimization

    # Update values
    model.create_grid()
    if ev_nul == 1:
        ev0 = np.zeros((model.n))     # CHANGED HERE
    else:
        ev0 = ev

    # Solve the model
    ev, pk, dev = solver.poly(model.bellman, V0=ev0 ,beta=model.beta, output=3)

    # Evaluate likelihood function
    lik_pr = pk[x]    
    log_lik = np.log(lik_pr+(1-2*lik_pr)*d+1e-15)       # + 1e-15 add a small number, which makes it more robust to log_lik(0)

    # add on log like for mileage process
    if theta.size>2:
        p = np.append(model.p,1-np.sum(model.p))
        if any(p<=0):
            log_lik -= 100000*p[dx1]
        else:
           log_lik += np.log(p[dx1])
        
    else:
        p = np.nan


    if out == 1:
        # Objective function (negative mean log likleihood)
        return np.mean(-log_lik)

    return model,lik_pr, pk, ev, dev, d,x,dx1


def score(theta, model, solver, data, pnames, ev_nul):
    global ev
    model,lik_pr, pk, ev, dev, d,x,dx1 = ll(theta, model, solver, data, pnames,ev_nul,9)
    F = np.eye(model.n)-dev    
    N = data.x.size
    dc = 0.001*model.grid
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

def grad(theta, model, solver,data, pnames,ev_nul):
    s = score(theta, model, solver, data,pnames,ev_nul)
    return -np.mean(s,0)


def hes(theta, model, solver,data, pnames,ev_nul):
    s= score(theta, model, solver, data, pnames,ev_nul)

    return s.T@s/data.shape[0]


def updatepar(par,parnames, parvals):

    for i,parname in enumerate(parnames):
        if i<2:
            parval = parvals[i]
            setattr(par,parname,parval)
        else:
            list_val = [None]*(parvals.size-2) 
            for j,parval in enumerate(parvals[2:]):
                list_val[j]=parval
            setattr(par,parname,list_val)
    return par

