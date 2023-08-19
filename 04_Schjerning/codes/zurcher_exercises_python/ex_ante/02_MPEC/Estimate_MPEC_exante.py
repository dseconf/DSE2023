#  MPEC class for structural estimation of discrete choice models.

# load general packages
import numpy as np
import Solve_NFXP as nfxp
import model_zucher as zucher
import scipy.optimize as optimize

def sparsity_pattern(Nc,N, M):
    """ This function creates the sparsity pattern for the Jacobian and Hessian of the Bellman equation.

        Arguments: 
            Nc: Number of structural parameters
            N: Number of states
            M: Number of states in the band-diagonal matrix
    
    """
    
    BandDiag = 0 # initalize band-diagonal matrix
    for i in range(M): 
        a = np.eye(N,k=i) # create matrix  with ones on the k-th diagonal (k=0 is the main diagonal and k>0 is diagonals above)
        BandDiag=BandDiag + a # add matrix to band-diagonal matrix
    BandDiag[:,0] = np.ones((N)) # set first column to ones
    
    JacobSpaPattern = np.append(np.ones((N,Nc)),BandDiag, axis=1) # append band-diagonal matrix to matrix with NC coulmns of ones
    
    Upper_block = np.ones((Nc+1,Nc+N)) # Construct upper block matrix
    Diagonal_block = np.eye(N-1) # Construct identity matrix to create diagonal of hessian 
    Left_block = np.ones((N-1,Nc+1))
    # Construct hessian of upper_block, identity_block and left_block
    HessSpaPattern = np.append(Upper_block, np.append(Left_block, Diagonal_block, axis=1), axis=0)
    return JacobSpaPattern, HessSpaPattern

def estimate(model, data, theta0 = [0,0],twostep=0):
    assert(twostep == 1),'MPEC only implemented for twostep=1'


    # Setup
    pnames = ['RC','c'] # Define parameter names

    # Step 1: Preprocess data
    class data_class: pass
    data_class.x = data.x - 1 # states in data
    data_class.dk = (data.d == 0) # Dummy for keep 
    data_class.dr = (data.d == 1) # Dummy for replace

    ## We create a dummy for being in a given state in the data
    data_class.xd = np.nan+np.zeros((model.n,data.x.size,)) # Initalize
    for i in range(model.n): # Loop over states
        data_class.xd[i,:] = (data.x == i+1) # Create dummy for being in state i+1 (data starts at 1)
    data_class.xd = data_class.xd.astype(int) # Convert to integer
    


    data = data_class # save data

    # Step 2: Estimate structual parameters 
    Nc = 2 # Number of estimated parameters
    J_pattern, _ = sparsity_pattern(Nc,model.n, len(model.p)+1) # Get sparsity pattern for Jacobian

    ## We are estimating both the parameters RC, c and the expected value EV for each grid point
    # set bounds
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

    # Define constraints and jacobian of constraint
    con_bell = lambda theta: con_bellman(theta,model,data, pnames) # Define constratint
    con_Jac = lambda theta: constraint_jac(theta,model,data, pnames)
    # Give sparsity pattern to constraint and set lower and upper bound to zero
    con_p_bellman = optimize.NonlinearConstraint(con_bell,0,0, jac = con_Jac, finite_diff_jac_sparsity = J_pattern) 

    theta0 = np.append(theta0,-np.ones((model.n)))  # Initial guess for parameters and expected value
    # Estimation
    res = optimize.minimize(ll,theta0, args=(model,data, pnames), method='trust-constr', jac=True, hess = '2-point', constraints =con_p_bellman, bounds = optimize.Bounds(lb, ub),options={'initial_constr_penalty': 1, 'xtol': 1e-10,'gtol': 1e-10, 'sparse_jacobian': True}) 
 
    theta_hat = res.x[0:2] # Estimated parameters

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
    ev = theta[2:]
    # Create grid
    model.create_grid()

    ## FILL IN BELOW

    # Value of options:

    #Hint: We are now solving in expected value function space. This means that unlike in exercise_1, 
    # you should not multiply ev with the transition matrix.
    #v_keep = # Value of keep
    #v_replace = # Value of replace
    
    # Choice probabiltiy of keep for all states
    #pk = 

    # Choice probability of keep given observed state in data
    lik_pr =  pk[data.x] 
    
    if out == 2:
        return model, lik_pr

    # Evaluate log-likelihood
    # log_lik = 

    # Objective function is negative mean log-likelihood
    f = -np.mean(log_lik)


    # GRADIENT: Is written as - (dummy for choice - choice-probability in data) * (value_replace - value_keep)    
    res =  - np.array(dk - lik_pr) # dummy for choice - choice-probability in data
    g = np.zeros((2+model.n)) # Initialize gradient
    g[0] = - np.mean(res)    # Gradient of RC
    g[1] =  np.mean(res*(model.dc[data.x]-model.dc[0])) # Gradient of c
    g[2] = - (model.beta * np.mean(res*(data.xd[0,:]-1)) ) # Gradient of first point in EV
    NT = res.size # Number of observations
    g[3:] = -model.beta*np.sum(np.multiply(data.xd[1:,:],res),1)/NT # Gradient of remaining points in EV

    return f, -g
    
def con_bellman(theta, model, data, pnames, out=1):
    
    # Update parameters
    ev0 = theta[-model.n:]

    # FILL IN BELOW
    ev1, pk, dev = bellman(model,ev0, output=3)

    if out ==2:
        return pk, dev

    return ev1-ev0

def constraint_jac(theta, model, data, pnames):
    
    pk,dev = con_bellman(theta, model, data, pnames, out=2) # Get pk and dev
    DCeq = np.zeros((model.n,2+model.n)) # Initialize Jacobian
    DCeq[:,0] = - model.P1 @(1-pk)  # Jacobian of RC
    DCeq[:,1] = -model.P1@(pk*(model.dc-model.dc[0])) # Jacobian of c
    DCeq[:,-model.n:] = dev-np.identity(model.n) # Jacobian of EV

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
