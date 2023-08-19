# Zurcher class: Contains model parts for Rust's engine repplacement model Rust(Ecta, 1987)

#import packages
import numpy as np
import time
import pandas as pd

class zurcher():
    def __init__(self,**kwargs):
        self.setup(**kwargs)

    def setup(self,**kwargs):     
  
        # a) parameters
        # Spaces
        self.n = 175                      # Number of grid points
        self.max = 450                    # Max of mileage

        # structual parameters
        self.p = np.array([0.0937, 0.4475, 0.4459, 0.0127])   # Transition probability
        self.RC = 11.7257                                     # Replacement cost
        self.c = 2.45569                                      # Cost parameter
        self.beta = 0.9999                                    # Discount factor

        # b. update baseline parameters using keywords
        for key,val in kwargs.items():
            setattr(self,key,val) 

        # c. Create grid
        self.create_grid()

    def create_grid(self):
        self.grid = np.arange(0,self.n) # milage grid
        self.cost = 0.001*self.c*self.grid  # cost function
        self.dc = 0.001*self.grid
        self.state_transition() 

    def state_transition(self):
        '''Compute transition probability matrixes conditional on choice'''
        p = np.append(self.p,1-np.sum(self.p)) # Get transition probabilities
        P1 = np.zeros((self.n,self.n)) # Initialize transition matrix
        # Loop over rows
        for i in range(self.n):
            # Check if p vector fits entirely
            if i <= self.n-len(p):
                P1[i][i:i+len(p)]=p
            else:
                P1[i][i:] = p[:self.n-len(p)-i]
                P1[i][-1] = 1.0-P1[i][:-1].sum()

        # conditional on d=1, replacement
        P2 = np.zeros((self.n,self.n))
        # Loop over rows
        for i in range(self.n):
            P2[i][:len(p)]=p
        self.P1 = P1
        self.P2 = P2
        
    def unc_state_transition(self,pk):
        #Calculate unconditional transition matrix
        self.Fu = self.P1*pk[:,np.newaxis] + self.P1[0,:]*(1-pk[:,np.newaxis]) 
        
        #Calculate [I-Fu]^{-1}
        self.Finv = np.linalg.inv(np.identity(self.n)-self.beta*self.Fu)
        
    def Psi(self, pk0,Finv=None):
        '''Evaluate Psi function (mapping from initial CCP to updated CCP)'''
        
        # Evaluate psi function (find Vsigma implied by initial CCP)
        self.psi(pk0, Finv)
        
        # Evaluate lambda function (find updated CCP implied by Vsigma)
        pk = self.lambdaa()
        
        return pk

    def psi(self,pk,Finv=None):
        '''Evaluate psi function (mapping from initial CCP to Vsigma)'''
        
        #Calculate [I-Fu]^{-1} if not given
        if Finv is None:
            self.unc_state_transition(pk)
        
        # Euler constant
        eulerc = np.euler_gamma
        
        # Compute Vsigma
        #Fill in
        #Hint:  Use the small psi function from slides
        #       Make sure to get E[\epsilon(a)|a,x] correct

        #self.Vsigma = 

    def lambdaa(self):
        '''Evaluate lambda function (mapping from Vsigma to updated CCP)'''
        
        #Compute updated CCP (pk)
        #Fill in
        #Hint:  Use the big lambda function from slides
        #       This is just the choice probability expression as used in NFXP
        
        #pk=
        
        return  pk

    def bellman(self,ev0,output=1):
        '''Evaluate Bellman operator, choice probability and Frechet derivative - written in integrated value form'''

        # Value of options:
        value_keep = -self.cost + self.beta * self.P1 @ ev0 # nx1 matrix
        value_replace = -self.RC - self.cost[0] + self.beta * self.P2 @ ev0   # 1x1

        # recenter Bellman by subtracting max(VK, VR)
        maxV = np.maximum(value_keep, value_replace) 
        logsum = (maxV + np.log(np.exp(value_keep-maxV)  +  np.exp(value_replace-maxV)))  # Compute logsum to handle expectation over unobserved states
        ev1 = logsum # Bellman operator as integrated value

        if output == 1:
            return ev1

        # Compute choice probability of keep
        pk = 1/(1+np.exp(value_replace-value_keep))       
        
        if output == 2:
            return ev1, pk

        # Compute derivative of Bellman operator
        dev1 = self.dbellman(pk)

        return ev1, pk, dev1

    def dbellman(self,pk): 
        '''Compute derivative of Bellman operator'''
        dev1 = np.zeros((self.n,self.n))
        for d in range(2): # Loop over choices 
            if d == 0:
                P = self.P1
                choice_prob =  pk
            else:
                P = self.P2
                choice_prob = 1-pk

            dev1 += self.beta * choice_prob.reshape(-1, 1) * P 
        
        return dev1

    def read_busdata(self, bustypes = [1,2,3,4]): 
        data = np.loadtxt(open("busdata1234.csv"), delimiter=",")
        idx = data[:,0]             # bus id
        bustype = data[:,1]         # bus type
        dl = data[:,4]              # laggend replacement dummy
        d = np.append(dl[1:], 0)    # replacement dummy
        x = data[:,6]               # Odometer

        # Discretize odometer data into 1,2,...,n
        x = np.ceil(x*self.n/(self.max*1000))

        # Montly mileage
        dx1 = x-np.append(0,x[0:-1])
        dx1 = dx1*(1-dl)+x*dl
        dx1 = np.where(dx1>len(self.p),len(self.p),dx1) # We limit the number of steps in mileage

        # change type to integrer
        x = x.astype(int)
        dx1 = dx1.astype(int)

        # Collect in a dataframe
        remove_first_row_index=idx-np.append(0,idx[:-1])
        data = {'id': idx,'bustype':bustype, 'd': d, 'x': x, 'dx1': dx1, 'boolean': remove_first_row_index}
        df= pd.DataFrame(data) 

        # Remove observations with missing lagged mileage
        df = df.drop(df[df.boolean!=0].index)

        # Select bustypes 
        for j in [1,2,3,4]:
            if j not in bustypes:
                df = df.drop(df[df.bustype==j].index) 

        # save data
        dta = df.drop(['id','bustype','boolean'],axis=1)
        
        return dta

    def sim_data(self,N,T,pk): 

        np.random.seed(2020)
        
        # Index 
        idx = np.tile(np.arange(1,N+1),(T,1))  
        t = np.tile(np.arange(1,T+1),(N,1)).T
        
        # Draw random numbers
        u_init = np.random.randint(self.n,size=(1,N)) # initial condition
        u_dx = np.random.rand(T,N) # mileage
        u_d = np.random.rand(T,N) # choice
        
        # Find states and choices
        csum_p = np.cumsum(self.p)
        dx1 = 0
        for val in csum_p:
            dx1 += u_dx>val
        
        x = np.zeros((T,N),dtype=int)
        x1 =  np.zeros((T,N),dtype=int)
        d = np.nan + np.zeros((T,N))
        x[0,:] = u_init # initial condition
        for it in range(T):
            d[it,:] = u_d[it,:]<1-pk[x[it,:]]  # replace = 1 , keep = 0   
            x1[it,:] = np.minimum(x[it,:]*(1-d[it,:]) + dx1[it,:] , self.n-1) # State transition, minimum to avoid exceeding the maximum mileage
            if it < T-1:
                x[it+1,:] = x1[it,:]
                
        
        # reshape 
        idx =  np.reshape(idx,T*N,order='F')
        t = np.reshape(t,T*N,order='F')
        d = np.reshape(d,T*N,order='F')
        x = np.reshape(x,T*N,order='F') + 1 # add 1 to make index start at 1 as in data - 1,2,...,n
        x1 = np.reshape(x1,T*N,order='F') + 1 # add 1 to make index start at 1 as in data - 1,2,...,n
        dx1 = np.reshape(dx1,T*N,order='F')


        data = {'id': idx,'t': t, 'd': d, 'x': x, 'x1': x1, 'dx1': dx1}
        df= pd.DataFrame(data) 

        return df

    def eqb(self, pk):
        # Inputs
        # pk: choice probability

        # Outputs    
        # pp: Pr{x} (Equilibrium distribution of mileage)
        # pp_K: Pr{x,i=Keep}
        # pp_R: Pr{x,i=Replace}
        tmp =self.P1[:,1:self.n] * pk[1:self.n]
        pl = np.hstack(((1-np.sum(tmp,1,keepdims=True)), tmp)) 

        pp = self.ergodic(pl)

        pp_K = pp.copy()    
        pp_K[0] = self.p[0]*pp[0]*pk[0]
        pp_R = (1-pk)*pp

        return pp, pp_K, pp_R

    def ergodic(self,p):
        #ergodic.m: finds the invariant distribution for an NxN Markov transition probability: q = qH , you can also use Succesive approximation
        n = p.shape[0]
        if n != p.shape[1]:
            print('Error: p must be a square matrix')
            ed = np.nan
        else:
            ap = np.identity(n)-p.T
            ap = np.concatenate((ap, np.ones((1,n))))
            ap = np.concatenate((ap, np.ones((n+1,1))),axis=1)

            # find the number of linearly independent columns
            temp, _ = np.linalg.eig(ap)
            temp = ap[temp==0,:]
            rank = temp.shape[1]
            if rank < n+1:
                print('Error: transition matrix p is not ergodic')
                ed = np.nan
            else:
                ed = np.ones((n+1,1))
                ed[n] *=2
                ed = np.linalg.inv(ap)@ed
                ed = ed[:-1]
                ed = np.ravel(ed)

        return ed
