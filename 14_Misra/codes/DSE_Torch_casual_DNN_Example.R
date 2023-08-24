# Simple Example estimating Heterogeneity using DNNs
# Sanjog Misra
# DSE Summer School 2023
# The model will be y = alpha(x) + beta(x)t + eps
# We will ignore cross fitting in this simple example
# The propensity score will be fixed (but estimated)

## ----------------------------------------------------------------------------------------------
# Load torch 
library(torch)

### generate training data -----------------------------------------------------
# input dimensionality (number of input features)
d_in <- 5
# output dimensionality (dim of [alpha,beta] : for now 2)
d_out <- 2
# number of observations in training set
n <- 10000

# Treatments
dim_z = 1L

# Set a seed
set.seed(1212)
torch_manual_seed(1212)

## ----------------------------------------------------------------------------------------------
# create random data
x = torch_tensor(matrix(rnorm(n*d_in),n,d_in))

#  True CATE
trub = 2-x[, 2, NULL]+.25*(x[, 1, NULL])^3
plot(as.numeric(trub)~as.numeric((x[, 1])),pch=19,col="#00000030")

# True Baseline
trua= x[, 1, NULL] * 0.2 - x[, 2, NULL] * 1.3 - x[, 3, NULL] * 0.5

# Binary Treatment z
z <- 0+1*(torch_randn(n, 1L)>0)
table(as.numeric(z))

# Outcomes are linear in treatments
y <- trua+ trub*z+ torch_randn(n, 1)

# Logit? Try changing to 
# trupr = 1/(1+exp(-trua-trub*z))
# y=rbinom(n,1,prob=as.numeric(trupr))
# Then convert torch tensor
# y = torch_tensor(matrix(y,n,1),dtype = torch_float64())


# The Estimator
## ----------------------------------------------------------------------------------------------
# The architecture (all layers are relu)
# Two layers with 20 nodes
arch=c(20,20)
# Can chan ge to what you need e.g.
# arch = c(10,10,10) 
# would eb three layers with 10 nodes each
# Ideally start with simple architectures then grow


# Model constructor
model <- nn_sequential(nn_linear(d_in, arch[1]),nn_relu())

j=1 # in case of single layer
# Loop through architecture
if(length(arch)>1){
  for(j in 2:length(arch)){
    model$add_module(name = paste("layer_",j-1),module=nn_linear(arch[j-1],arch[j]))
    model$add_module(name = paste0("relu_",j-1),nn_relu())
  }
}

# Output layer
model$add_module("ab",nn_linear(arch[j],d_out))
# Can print this to see what's in it
model

## ----------------------------------------------------------------------------------------------
# Optimization/Learning framework
# for ADAM, need to choose a reasonable learning rate (subjective)
learning_rate <- 0.01
optimizer <- optim_adam(model$parameters, lr = learning_rate)
# We use adam but you can change as needed


## ----------------------------------------------------------------------------------------------
# Training Loop
NEPOCH = 1000
intv = NEPOCH/100
cat("Begining training...\n")
pb <- txtProgressBar(min=1,max=100,style=3)
pct=0
st = proc.time()

for (t in 1:NEPOCH) {

  ### -------- Forward pass --------
  ### Causal Model
  ab <- model(x)
  alpha=ab[,1]$reshape(c(n,1L))
  beta=ab[,2]$reshape(c(n,1L))
  ab$retain_grad()
  alpha$retain_grad()
  beta$retain_grad()
  y_pred <- alpha+beta*z
  loss <- nnf_mse_loss(y_pred, y, reduction = "mean")

  ### -------- Backpropagation --------
  # Need to zero out the gradients before the backward pass
  # otherwise it accumulates
  optimizer$zero_grad()

  # gradients are computed on the loss tensor
  loss$backward()

  ### -------- Update weights --------
  # use the optimizer to update model parameters
  optimizer$step()
  # progress bar update
  if(t%%intv==0) {pct=pct+1; setTxtProgressBar(pb, pct)}
}
et = proc.time()
elapsed = (et-st)[3]
elapsed


# Now let's do inference
## ----------------------------------------------------------------------------------------------
# predict theta(x) we call it ab here
ab <- model(x)
# Split
alpha=ab[,1]$reshape(c(n,1L))
beta=ab[,2]$reshape(c(n,1L))
# We'll do some processing to tell torch we'll be taking gradients
ab$retain_grad()
alpha$retain_grad()
beta$retain_grad()

# Do we recover function?
plot(as.numeric(trub)~as.numeric((x[, 1])),
     pch=19,col='#00000030',xlab="X",ylab="beta")
points(as.numeric(beta)~as.numeric((x[, 1])),pch=19,col='#ff000030')

# Inference
## ----------------------------------------------------------------------------------------------
# Predict y
y_pred <- alpha+beta*z

# Compute loss
loss <- .5*nnf_mse_loss(y_pred, y, reduction = "none")

# Gradient wrt ab
# Trick! we sum loss and take derivatives which gives us the correct answer 
abg = autograd_grad(loss$sum(),ab,create_graph = TRUE,retain_graph = TRUE)

# Extract dim
K=dim(abg[[1]])[2]

# Trick again! sum gradient and take derivatives for hessian
abg1 = abg[[1]]$sum(1)

# Hessian 
# Here we usea single unconditional hessian since treatment is randomized
hess=matrix(0,K,K)
for(j in 1:K){
  tmp=autograd_grad(abg1[j],ab,retain_graph = TRUE)[[1]]
  hess[,j] =  as.numeric(tmp$mean(1))
}

# Lambda inverse
lami=-solve(hess)

# Save gradients
gg = as.array(abg[[1]])


# Automatically differentiate H construct
ab <- model(x)
ab$retain_grad()
Hi = ab[,2]

H = Hi$sum()
Hab = autograd_grad(H,ab)
Hab2 = (as.array(Hab[[1]]))

# Compute IF adjusted estomator
psi.auto = as.numeric(Hi)+sapply(1:n,function(j) Hab2[j,]%*%lami%*%(t(gg)[,j]))
# Estimator and se
au.est= mean(psi.auto)
au.se= sqrt(var(psi.auto)/n)

# Doubly Robust SE
Z = as.numeric(z)
Y = as.numeric(y)
Yhat = as.numeric(y_pred)
mu0 = as.numeric(alpha)
mu1 = mu0+as.numeric(beta)
e = as.numeric(mean(Z))

# Doubly Robust Analytical Formula
psi.dr = (mu1-mu0) +Z*((Y-mu1)/e) - (1-Z)*((Y-mu0)/(1-e))
DR.est = mean(psi.dr)
DR.se = sqrt(var(psi.dr)/n)
DR.se

dr = c(Est=DR.est,se = DR.se,
       CI.L=DR.est-1.96*DR.se,CI.U=DR.est+1.96*DR.se)
au = c(Est=au.est,se = au.se,
       CI.L=au.est-1.96*au.se,au.U=au.est+1.96*au.se)

rbind(DRobust=dr,Automatic =au)


# For the sake of comparison lets run a causal forest
library(grf) # install if needed
cf = causal_forest(X =as.matrix(x),Y = Y,W=Z,num.trees = 100 )
tau = (cf$predictions)

# Plot densities
mm = data.frame(Estimate = c(as.numeric(trub),as.numeric(tau),as.numeric(beta)))
mm$model = c(rep('Truth',1000),rep('Causal Forest',1000),rep('Deep Net',1000))
library(ggplot2)
library(ggridges)

ggplot(data=mm,aes(x=Estimate,y=model))+
  geom_density_ridges(fill = "lightblue") +
    scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  theme_ridges()+ theme(legend.position = "none")
