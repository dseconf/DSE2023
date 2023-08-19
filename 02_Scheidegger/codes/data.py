import numpy as np
import pandas as pd

class Data:
    def __init__(self):
        self.X_tr = None
        self.X_te = None
        self.X_va = None

        self.y_tr = None
        self.y_te = None
        self.y_va = None

    def load_lending_club_data(self):
        X = pd.read_csv('data/LCx.csv', index_col=0)
        y = pd.read_csv('data/LCy.csv', index_col=0)

        N = X.shape[0]
        ##################
        # split the data
        ##################
        # define training and testing size in % of total sample
        tr_p = 0.9
        te_p = 0.05
        # then transform these numbers into number of sample
        tr = int(np.ceil(tr_p*N))
        te = int(np.ceil(te_p*N))

        # transform the pandas data frame into numpy
        X = X.values
        y = y.values

        # shuffle the indexes
        ind = np.arange(0,N)
        np.random.shuffle(ind)

        self.X_tr = X[:tr,:]
        self.X_te = X[tr:(tr+te),:]
        self.X_va = X[(tr+te):,:]

        self.y_tr = y[:tr,:]
        self.y_te = y[tr:(tr+te),:]
        self.y_va = y[(tr+te):,:]


    def gen_random_data(self, N = 100*1000, D = 10, LL = [10,10], output_dim = 1, sigma=1):
        assert  type(LL) == list, 'the complexity parameter LL should be a python list []'
        assert output_dim>0, 'ouput_dim must be larger than 0'

        # never use random function without a seed, science is replicable!
        np.random.seed(1234)

        # we create some random data by computing a pseudo inverse neural network
        # its pseudo because the activation are not inverted. You can figure out why with pen and paper.
        # we start by defining the target y
        y = np.random.normal(size=(N,output_dim))
        if output_dim > 1:
            # if we have more than one dim, we pass y into a softmax function to get probabilities
            e = np.exp(y)
            y = e/np.tile(np.sum(e,1).reshape(-1, 1), output_dim)

        # we add to the inverse archtecture D, the input dimension X to get this at the end
        LL.append(D)
        # we initialize our X with the value of y
        X = y
        for i, L in enumerate(LL):
            # for each layer L in LL, we generate random weights
            W = np.random.uniform(-1, 1, size=(X.shape[1],L))  # randomly draw the weight matrix
            X = X @ W
            # X = 1/(1+np.exp(-X))
            X=np.maximum(X,0)
            # if i ==0:
            #     # if we are not at the first layer, we compute a relu
            #     X = np.maximum(X@W,0)
            # else:
            #     X = X@W

        y = y + np.random.normal(scale=sigma, size=y.shape)
        if output_dim>1:
            # if we have more than one dim, set highest to 1 (correct class) and other to 0
            e = np.max(y,1)
            y = (np.tile(e.reshape(-1,1),output_dim)==y)*1

        ##################
        # split the data
        ##################
        # define training and testing size in % of total sample
        tr_p = 0.9
        te_p = 0.05
        # then transform these numbers into number of sample
        tr = int(np.ceil(tr_p*N))
        te = int(np.ceil(te_p*N))

        # shuffle the indexes
        ind = np.arange(0,N)
        np.random.shuffle(ind)

        self.X_tr = X[:tr,:]
        self.X_te = X[tr:(tr+te),:]
        self.X_va = X[(tr+te):,:]

        self.y_tr = y[:tr,:]
        self.y_te = y[tr:(tr+te),:]
        self.y_va = y[(tr+te):,:]
