# load general packages
import numpy as np
import scipy.sparse as sparse
import scipy.optimize as optimize

# Import 
pk = 0.0 
def setup_data(data): 
   
    # setup
    class data_class: pass
    data_class.x = data.x -1
    data_class.dx1= data.dx1
    data_class.dk = (data.d == 0)
    data_class.dr = (data.d == 1)
    return data_class

def estimate(model, data, Kmax = 100):
    '''Estimate model using Nested Psuedo Likelihood (NPL)'''
    #Load previous choice probability as global variable
    global pk    

    # Find transition probabilities, p (non-parametrical 1-step estimator)
    tabulate = data.dx1.value_counts() # Count number of observations for each dx1
    p = [tabulate[i]/sum(tabulate) if i < len(tabulate) else 0 for i in range(len(model.p))]
    model.p[:] = p # Use first step estimates as starting values for p
    
    
    # Find state transition matrix
    model.state_transition()  
    
    #Set starting valiues
    pk0 = np.ones((model.n))*0.99  # starting value for CCP's
    theta0 = [0,0] # starting value for parameters 

    #Outer loop: Will update CCPs until convergence or Kmax iterations are reached
    for _ in range(Kmax): #

        # Step 0)  Pre-compute unconditional transition matrix Fu and the inverse of I-beta*Fu to used in step 1 
        model.unc_state_transition(pk0) # Calculate Fu and Finv

        #Inner loop:
        # Step 1)  Maximize the pseudo-likelihood function given step K-1 CCPs
        res = optimize.minimize(ll,theta0,args =(model, data, pk0), method='Newton-CG', jac = grad, hess= hes, tol = 1e-6)
        theta_hat = res.x # save parameters
        NPL_metric = np.abs(theta0-theta_hat) #save distance between parameters

        #Outer loop step
        # Step 2)  Update CCPs using theta_npl from step 1)
        pk0 = pk
        theta0 = theta_hat
        if NPL_metric.all() < 1e-6: # check convergence
            return res, theta_hat,pk
    
    print(f'The function did not converge after {K} iterations')
    
    return res, theta_hat, pk

def ll(theta, model, data, pk0,out=1):
    '''Log-likelihood function for NPL'''
    #Load previous choice probability as global variable
    global pk

    # update parameters
    model.RC = theta[0]
    model.c = theta[1]
    model.create_grid()

    # Update CCPs
    #Fill in
    #pk =
    
    #Map choice probabilities to data
    pKdata = pk[data.x] 

    #Return CCPs if out = 2
    if out == 2:
        return pk, pKdata

    #Calculate log-likelihood
    #Fill in
    #log_lik = 
    
    #Return log-likelihood if out = 1
    f = -np.mean(log_lik)
    return f

def score(theta, model, data, pk0):
    global pk
    pk, pkdata = ll(theta, model, data, pk0,out=2)
    res = data.dk-pkdata

    NT = ((pkdata.size))
    score = np.zeros((NT,theta.size))
    dP = model.P1[0,:]-model.P1  

    dvdRC = np.ravel(-1+model.beta*dP@model.Finv@(1-pk[:,np.newaxis])*(-1))
    dvdc = np.ravel(model.dc + model.beta*dP@model.Finv@(pk*(-model.dc)))
    score[:,0] = res*dvdRC[data.x]
    score[:,1] = res*dvdc[data.x]

    return score 

def grad(theta, model, data, pk0):
    s = score(theta, model, data, pk0)
    return np.mean(s,0)

def hes(theta, model, data, pk0):
    s = score(theta, model, data, pk0)
    return s.T@s/data.x.shape[0]

def solve(model):
    '''Solve model by successive approximation in choice probabilitiy space'''
    #Set starting valiues
    pk0 = np.ones((model.n))*0.99 
    
    #Allocate
    pk = np.nan+np.zeros((100,model.n))
    pk[0,:] = pk0
    
    #Update grids
    model.create_grid()

    #Solve using successive approximation
    for i in range(1,100): 
        pk[i,:]  = model.Psi(pk[i-1,:])    
    return pk
