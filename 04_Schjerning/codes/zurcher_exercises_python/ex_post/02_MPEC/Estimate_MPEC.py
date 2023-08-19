#  MPEC class for structural estimation of discrete choice models.

# load general packages
import numpy as np
import Solve_NFXP as nfxp
import model_zucher as zucher
import scipy.optimize as optimize

def sparsity_pattern(Nc,N, M):
    BandDiag = 0
    for i in range(M):
        a = np.eye(N,k=i)
        BandDiag=BandDiag + a
    BandDiag[:,0] = np.ones((N))
    
    JacobSpaPattern = np.append(np.ones((N,Nc)),BandDiag, axis=1)
    
    HessSpaPattern = np.ones((Nc+1,Nc+N))
    HessSpaPattern = np.append(HessSpaPattern, np.append(np.ones((N-1,Nc+1)),np.eye(N-1),axis=1), axis=0)

    return JacobSpaPattern, HessSpaPattern

def estimate(model, data, theta0 = [0,0],twostep=0):
    assert(twostep == 1),'MPEC only implemented for twostep=1'

    # Setup
    pnames = ['RC','c']

    class data_class: pass

    data_class.x = data.x - 1
    data_class.dk = (data.d == 0)
    data_class.dr = (data.d == 1)
    data_class.xd = np.nan+np.zeros((model.n,data.d.size,))
    for i in range(model.n):
        data_class.xd[i,:] = data.x == i + 1
    data_class.xd = data_class.xd.astype(int)
    
    
    data = data_class

    # Step 2: Estimate structual parameters 
    Nc = 2
    J_pattern, _ = sparsity_pattern(Nc,model.n, len(model.p)+1)

    # bounds
    lb = np.zeros((2+model.n))
    ub = np.zeros((2+model.n))

    #bound on c and RC
    lb[0] = 0
    ub[0] = np.inf
    lb[1] = 0
    ub[1] = np.inf
  

    # bounds on EV
    lb[-(model.n):] = -5000
    ub[-(model.n):] = 0

    # Define the objective functions and constraints
    con_bell = lambda theta: con_bellman(theta,model,data, pnames) # Define constratint
    con_Jac = lambda theta: constraint_jac(theta,model,data, pnames)
    con_p_bellman = optimize.NonlinearConstraint(con_bell,0,0, jac = con_Jac, finite_diff_jac_sparsity = J_pattern) 

    theta0 = np.append(theta0,-np.ones((model.n)))  
    res = optimize.minimize(ll,theta0, args=(model,data, pnames), method='trust-constr', jac=True, hess = '2-point', constraints =con_p_bellman, bounds = optimize.Bounds(lb, ub),options={'initial_constr_penalty': 1, 'xtol': 1e-10,'gtol': 1e-10, 'sparse_jacobian': True}) 
 
    theta_hat = res.x[0:2]

    return res, pnames, theta_hat

def ll(theta,model,data,pnames,out=1):
    
    # Unpack
    x = data.x
    xd = data.xd
    dk = data.dk
    dr = data.dr
    # Update values
    model.RC = theta[0]
    model.c = theta[1] 
    ev = theta[-model.n:] 
    model.create_grid()

    # Value of options:
    value_keep = -model.cost + model.beta*ev
    value_replace = -model.RC - model.cost[0] + model.beta*ev[0] 
    pk = 1/(1+np.exp(value_replace-value_keep))  

    # Evaluate the likelihood function 
    lik_pr = pk[x] 

    if out == 2:
        return model, lik_pr
        
    log_lik = np.log(dk*lik_pr+(1-lik_pr)*dr)
    f = -np.mean(log_lik)


    # GRADIENT    
    res = np.array(lik_pr-dk)
    g = np.zeros((2+model.n))
    g[0] = - np.mean(res)    # RC
    g[1] =  np.mean(res*(model.dc[x]-model.dc[0]))  # c
    g[2] = - (model.beta * np.mean(res*(xd[0,:]-1)) ) # ev(0) xd[:,0]-1
    NT = res.size
    g[3:] = -model.beta*np.sum(np.multiply(xd[1:,:],res),1)/NT

    return f, -g

def con_bellman(theta, model, data, pnames, out=1):
    
    # Update parameters
    ev0 = theta[-model.n:]

    ev1, pk, dev = bellman(model, ev0=ev0,output=3)
    if out ==2:
        return pk, dev

    return ev1-ev0

def constraint_jac(theta, model, data, pnames):
    
    pk,dev = con_bellman(theta, model, data, pnames, out=2)
    DCeq = np.zeros((model.n,2+model.n))
    DCeq[:,0] = - model.P1 @(1-pk) 
    DCeq[:,1] = -model.P1@(pk*(model.dc-model.dc[0]))
    DCeq[:,-model.n:] = dev-np.identity(model.n)

    return DCeq





def bellman(model,ev0=np.zeros(1),output=1):
    """" 
    Bellman operator in expected value - operator in model_zucher is in integrated value"""

    # Value of options:
    value_keep = -model.cost + model.beta*ev0
    value_replace = -model.RC - model.cost[0] + model.beta*ev0[0]  

    
    # recenter Bellman by subtracting max(VK, VR)
    maxV = np.maximum(value_keep,value_replace) 
    logsum = (maxV + np.log(np.exp(value_keep-maxV)  +  np.exp(value_replace-maxV)))  # This is the Logsum 
    ev1 = model.P1@logsum

    if output == 1:
        return ev1

    # Compute choice probability
    pk = 1/(1+np.exp(value_replace-value_keep))       
    
    if output == 2:
        return ev1, pk

    # Compute Frechet derivative
    dev1 = dbellman(model, pk)

    return ev1, pk, dev1

def dbellman(model,pk): 
    dev1 = model.beta * model.P1 * pk.transpose()
    dev1[:,0] += model.beta * model.P1 @ (1-pk)
    
    return dev1



