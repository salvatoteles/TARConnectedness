setwd("~/Doctorado/Artículos/Connectedness/TFG/EconomicLetters/LatentPanelGravity")
library(pracma)
data <- read.csv('my_spill_data.csv', header = TRUE,sep=";")
data[is.na(data)] <- 0

N <- round(sqrt(dim(data)[1]))  # number of countries
T <- N
K <- dim(data)[2] - 4  # number of regressors

# define NxN matrix of outcomes Y, and KxNxN multi-matrix of regressors X:
Y0 <- matrix(data[, 3], nrow = N, ncol = N)
Y1 <- exp(Y0)*(Y0 > 0)
Y <- Y1 / sd(Y1)

X <- array(0, dim = c(K, N, N))
for (k in 1:K) {
  X[k, , ] <- matrix(data[, 4+k], nrow = N, ncol = N)
}

# Eliminar el primer regresor
# X <- X[-1,,]
# X <- X[-10,,]
# X <- X[-8,,]
# K <- K - 3


# redefine island and landlock regressor:
  
  #X(K-2,:,:)=( squeeze(X(K-2,:,:))==2 );   %island is dummy for island == 2 (i.e. both islands)
#X(K,:,:)=( squeeze(X(K,:,:))==2 );       %landlock is dummy for landlock == 2 (i.e. both landlocked)

# check which regressors are binary (just for information at this point)
# binary_regressor_list <- c()
# 
# for (k in 1:K) {
#   if (all(isTRUE(X[k,,]== as.integer(as.logical(X[k,,]== 1))))) {
#     binary_regressor_list <- c(binary_regressor_list, k)
#   }
# }
# 
# binary_regressor_list
binary_regressor_list <- c(2,3)

# add constant
# X2 <- array(0, dim = c(K+1, N, N))
# X2[1:K,,] <- X
# 
# X2[K+1,,] <- matrix(1, nrow = N, ncol = N)
# K <- K + 1

lambda_known <- matrix(1, nrow = N, ncol = 1)  # known factor loading, corresponding to time dummies
f_known <- matrix(1, nrow = N, ncol = 1)  # known factors, corresponding to individual dummies
Rex1 <- dim(lambda_known)[2]
Rex2 <- dim(f_known)[2]
weight <- matrix(1, nrow = N, ncol = N) - diag(N)  # all diagonal observations are missing, so get zero weight


precision_beta <- 10^-2; repMIN<- 100; repMAX<- 101; MAX_STEPS<- 450;
R<-0           #number of additional factors
dist='Poisson';            #logit or probit or poisson 

xt <- 1:T
One <- rep(1, T)
X0 <- matrix(xt, ncol = 1)

#based on orthogonal polinomials of 1:T
if (R > 1) {
  for (i in 2:R) {
    X0 <- cbind(X0, xt^i)
  }
}

library(MASS)
Xtt <- One
for (i in 1:ncol(X0)) {
  Xtt <- cbind(Xtt, X0[,i] - Xtt %*%ginv(Xtt)%*% X0[,i])
}
fhat <- Xtt[, 2:(R+1)]

# (ii) Based on Principal Components 
# Obtain initial values by linear principal components 
# YtY <- t(Y) %*% Y
# eigen_result <- eigen(YtY)
# fhat2 <- eigen_result$vectors[, 1:R]  # eigenvectors corresponding to the R largest eigenvalues
# D <- eigen_result$values[1:R]  # eigenvalues corresponding to the R largest eigenvalues
# fhatN <- fhat2 %*% sqrt(T)


# Estimation: without bias correction
#[beta,alpha,gamma,lambda,f,exitflag,obj,APE,Var_beta,Var_APE] = FactorMLE(dist, Y, X, weight, lambda_known, f_known, R, precision_beta, repMIN, repMAX, MAX_STEPS);
lolo <- FactorMLE(dist, Y, X, weight, lambda_known, f_known, R, precision_beta, repMIN, repMAX, MAX_STEPS, jacknife=FALSE);

beta <- lolo$beta

obj <- lolo$obj*N*N / sum(sum(weight));

std_beta<-sqrt(diag(lolo$Var_beta))
lolo$APE
std_APE<-sqrt(diag(lolo$Var_APE))



FactorMLE <- function(dist, Y, X, weight, lambda_known, f_known, R, precision_beta, repMIN, repMAX, MAX_STEPS, jacknife=FALSE,beta_start, alpha_start, gamma_start, lambda_start, f_start) {

  
 #  %%%% Maximum Likelihood Estimation of Panel Factor Models
 #  %%%% This version October 6th
 #  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 #    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 #    % INPUT PARAMETERS: 
 #    %     dist = 'probit' for probit model
 #    %          = 'logit' for logit model
 #    %          = 'LS' for (weighted) least squares estimation (i.e. Gaussian MLE)
 #    %          = 'Poisson' for Poisson
 #    %     Y = NxT matrix of outcomes
 #    %     X = KxNxT multi-matrix of regressors
 #    %     weight = NxT matrix that specifies fixed weights for each observation,
 #    %              for standard MLE we set weight(i,t)=1 if (i,t) observed
 #    %              and weight(i,t)=0 if (i,t) unobserved,
 #    %              but more general weights are also possible.
 #    %     lambda_known = N x Rex1 matrix of KNOWN factor loadings,
 #    %                    e.g. lambda_known=ones(N,1) to control standard time dummies
 #    %         ... if NO KNOWN factor loadings, then set lambda_known=zeros(N,0) !!!
 #      %     f_known = T x Rex2 matrix of KNOWN factors,
 #      %                    e.g. f_known=ones(T,1) to control standard individual specific fixed effects
 #      %         ... if NO KNOWN factors, then set f_known=zeros(T,0) !!!
 #        %     R = positive integer 
 #        %         ... number of interactive fixed effects in the estimation
 #        %     precision_beta = defines stopping criteria for numerical optimization,
 #        %                      namely optimization is stopped when difference in beta
 #        %                      relative to previous opimtization step is smaller than
 #        %                      "precision_beta" (uniformly over all K components of beta)
 #        %     repMIN = positive integer
 #        %              ... minimal number of runs of optimization with different starting points
 #        %     repMAX = positive integer
 #        %              ... maximal number of runs of optimization (in case numerical optimization
 #                                                                   %              doesn't terminate properly, we do multiple runs even for repMIN=1)
 # %     MAX_STEPS = maximum number of iteration steps for each optimization
 # %                 run (e.g. MAX_STEPS=300)
 # %
 # % OPTIONAL INPUT PARAMTERS:
 # %     beta_start, alpha_start, gamma_start, lambda_start, f_start
 # %           ... starting values for first numerical optimization
 # %               (afterwards random starting values)
 # %           ... same dimensions as output parameters below
 # %
 # % OUTPUT PARAMETERS:
 # %     beta   = Kx1 vector of parameter estimate 
 # %     alpha  = N x Rex2 matrix of factor loadings corresponding to f_known
 # %     gamma  = T x Rex1 matrix of factors corresponding to lambda_known 
 # %     lambda = N x R matrix of estimates for factor loading
 # %     f      = T x R matrix of estimates for factors
 # %     exitflag = 1 if iteration algorithm properly converged at optimal beta
 # %              = -1 if iteration algorithm did not properly converge at optimal beta
 # %     obj = objective function at optimum
 # %     Var_beta = estimated variance-covariance matrix of beta
 # %     APE = Kx1 vector of average partial effects
 # %     Var_APE = estimated variance-covariance matrix of APE's
 #   %     beta_corrected = bias corrected estimator for beta
 #  %     APE_corrected = bias corrected estimator for APE's
  
  
  K <- dim(X)[1]  # number of regressors
  N <- dim(X)[2]  # cross-sectional dimension
  T <- dim(X)[3]  # time-serial dimension
  Rex1 <- dim(lambda_known)[2]
  Rex2 <- dim(f_known)[2]
  
  # Verificar que se especifique 'logit', 'probit', 'LS' o 'Poisson'
  if (!(dist == "logit" || dist == "probit" || dist == "LS" || dist == "Poisson")) {
    cat("ERROR in FactorMLE: el primer argumento de entrada debe ser 'logit', 'probit', 'Poisson' o 'LS'.\n")
    exitflag <- -1
    return()
  }
  
  # if (missing(beta_start)) {
  #   beta <- beta_start
  #   alpha <- alpha_start
  #   gamma <- gamma_start
  #   lambda <- lambda_start
  #   f <- f_start
  #   
  #   if (dim(beta)[1] != K || dim(beta)[2] != 1 ||
  #       dim(alpha)[1] != N || dim(alpha)[2] != dim(f_known)[2] ||
  #       dim(gamma)[1] != T || dim(gamma)[2] != dim(lambda_known)[2] ||
  #       dim(lambda)[1] != N || dim(lambda)[2] != R ||
  #       dim(f)[1] != T || dim(f)[2] != R) {
  #     exitflag <- -1
  #     cat("ERROR in FactorMLE: PARAMETER STARTING VALUES HAVE WRONG DIMENSIONS\n")
  #     return()
  #   }
  # } else
  beta <- matrix(0, nrow = K, ncol = 1)
  alpha <- matrix(0, nrow = N, ncol = dim(f_known)[2])
  gamma <- matrix(0, nrow = T, ncol = dim(lambda_known)[2])
  if(R==0){
    lambda <- matrix(vector(), nrow = N, ncol = R)
    f <- matrix(vector(), nrow = T, ncol = R)
  } else {
  lambda <- matrix(0, nrow = N, ncol = R)
  f <- matrix(0, nrow = T, ncol = R)
  }
  #}
  
  # Asegurarse de que todos los valores de los regresores para datos faltantes sean iguales a cero
  for (i in 1:N) {
    for (t in 1:T) {
      if (weight[i, t] == 0) {
        X[, i, t] <- 0
      }
    }
  }
  
  
  if(missing(repMIN)){
    variance_requested=0;
  } else{
    variance_requested=1;
  }
  # if nargout>8
  # variance_requested=1;   %user wants variance
  # else
   # variance_requested=0;
  # end  
  
    
    # OPTIMIZACIÓN NUMÉRICA:
    
    obj_best <- -Inf
    
    count_successful_runs <- 0
    count_total_runs <- 0
    
    while ((count_successful_runs < repMIN) && (count_total_runs < repMAX)) {
      count_total_runs <- count_total_runs + 1
      # ITERACIÓN:
      
      step_size <- Inf
      count <- 0
      count_final <- 0
      v_reference <- NULL
      para_reference <- NULL
      b_reference <- 0
      
      while (count_final <= 5 && count < MAX_STEPS) {
        count <- count + 1
        betaOLD <- beta
        
        if (count <= 5) {
          cat("hola")
          result <- step(dist, Y, X, weight, lambda_known, f_known, beta, alpha, gamma, lambda, f)
          beta <- result$beta
          alpha <- result$alpha
          gamma <- result$gamma
          lambda <- result$lambda
          f <- result$f
          obj <- result$obj
        } else {
          cat("adios")
          result <- stepNR(dist, Y, X, weight, lambda_known, f_known, beta, alpha, gamma, lambda, f, v_reference, para_reference, b_reference)
          beta <- result$beta
          alpha <- result$alpha
          gamma <- result$gamma
          lambda <- result$lambda
          f <- result$f
          obj <- result$obj
          v_reference <- result$v_reference
          para_reference <- result$para_reference
          b_reference <- result$b_reference
        }
        
        step_size <- max(abs(beta - betaOLD))
        
        if (step_size <= precision_beta) {
          count_final <- count_final + 1
        } else {
          count_final <- 0
        }
        
        cat("Iteration:", count, "\n")
        cat("Objective:", obj, "\n")
      }
      
      # CONTAR ITERACIONES EXITOSAS:
      status <- -1
      if (count < MAX_STEPS) {
        count_successful_runs <- count_successful_runs + 1
        status <- 1
      }
      
      # VERIFICAR SI EL OBJETIVO ES MEJOR QUE EL MEJOR OBJETIVO ANTERIOR,
      # SI ES ASÍ, GUARDAR ESOS PARÁMETROS:
      if (obj > obj_best) {
        obj_best <- obj
        beta_best <- beta
        alpha_best <- alpha
        gamma_best <- gamma
        lambda_best <- lambda
        f_best <- f
        exitflag <- status  # exitflag indica si esos valores óptimos de los parámetros corresponden a una ejecución de iteración "exitosa", con la definición de "exitosa" anterior
      }
      
      # GENERAR NUEVOS VALORES INICIALES ALEATORIOS:
      #beta <- 2 * runif(K) - 1
      beta <- matrix(2 * runif(K) - 1, ncol = 1)
      #alpha <- 2 * runif(N, min = -1, max = 1)
      alpha <- matrix(2*runif(N)-1,ncol=dim(f_known)[2])
      gamma <- matrix(2 * runif(T) -1, ncol = dim(lambda_known)[2])
      if(R==0){
        lambda <- matrix(vector(), nrow = dim(X)[2], ncol = 0)
        f <- matrix(vector(), nrow = dim(X)[2], ncol = 0)
      }else{
      lambda <- matrix(2 * runif(R*N) - 1,ncol=R)
      f <- matrix(2 * runif(R*T) - 1,ncol=R)
      }
    }
  
    ##########################hasta aqui estamos bien
    
    #REPORT PARAMETERS WITH MAXIMUM OBJECTIVE THAT WAS FOUND:
    beta<-beta_best;
    alpha<-alpha_best;
    gamma<-gamma_best;
    lambda<-lambda_best;
    f<-f_best;
    
    #CALCULATE VARAINCE OF ESTIMATOR VIA HESSIAN:
    if (variance_requested==1){
    obj_score_hessian <- SampleLogL(dist, Y, X, weight, lambda_known, f_known, beta, alpha, gamma, lambda, f, 0)
    obj <- obj_score_hessian$obj
    score <- obj_score_hessian$s
    Hessian <- obj_score_hessian$H
    
    V <- pracma::pinv(-Hessian * sum(weight))
    Var_beta <- V[1:K,1:K]
    
    }else{
      obj_score_hessian <- SampleLogL(dist, Y, X, weight, lambda_known, f_known, beta, alpha, gamma, lambda, f, 0)
      obj <- obj_score_hessian$obj
    }

    Z <- mult(X,beta) + alpha %*% t(f_known) + lambda_known %*% t(gamma) + lambda %*% t(f);
    
    if (dist == "logit" | dist == "probit" | dist == "Poisson") {
      # Cálculo del APE para REGRESORES CONTINUOS (DE HECHO, LO USAMOS PARA TODOS LOS REGRESORES NO BINARIOS):
      if (dist == "logit") {
        Z <- pmax(pmin(Z, 100), -100)  # truncar el índice en -100 y 100
        bc_pdf <- exp(Z) / (1 + exp(Z))^2
        bc_der_pdf <- exp(Z) * (1 - exp(Z)) / (1 + exp(Z))^3
      } else if (dist == "probit") {
        Z <- pmax(pmin(Z, 6), -6)  # truncar el índice en -6 y 6
        bc_pdf <- dnorm(Z)
        bc_der_pdf <- dnorm(Z) * dnorm(-Z)
      } else {
        bc_pdf <- exp(Z)
        bc_der_pdf <- exp(Z)
      }
    
    
    bc_pdf <- weight * bc_pdf
    bc_der_pdf <- weight * bc_der_pdf
    
    F1 <- sum(sum(bc_pdf)) / sum(sum(weight))
    APE_CONT <- beta * F1
    
    #calculate variance of APE by delta method (using Jacobian J)
    if (variance_requested == 1) {
      s <- rep(0, K + N * (R + Rex2) + T * (R + Rex1))  # score of F1 = average(F'(Z))
      ff <- cbind(f, f_known)
      ll <- cbind(lambda, lambda_known)
      
      for (k in 1:K) {
        XX <- X[k,,]
        s[k] <- sum(sum(bc_der_pdf * XX)) / sum(sum(weight))
      }
      
      for (r in 1:(R + Rex2)) {
        s[K + (r-1) * N + (1:N)] <- bc_der_pdf %*% ff[, r] / sum(sum(weight))
      }
      
      for (r in 1:(R + Rex1)) {
        s[K + N * (R + Rex2) + (r-1) * T + (1:T)] <- t(bc_der_pdf) %*% ll[, r] / sum(sum(weight))
      }
      
      #aqui hay que crear una matriz 322x8 con las 8 primeras de la diagonal ocupadas
      matriz_suple <- matrix(0, nrow = K + N * (R + Rex2) + T * (R + Rex1), ncol = K)
      matriz_suple[1:K,1:K] <- diag(rep(1,K))
      #J <- diag(K + N * (R + Rex2) + T * (R + Rex1), K) * F1 + matrix(s, nrow = K + N * (R + Rex2) + T * (R + Rex1), ncol = K)
      J <- matriz_suple*F1+ matrix(s, nrow = K + N * (R + Rex2) + T * (R + Rex1), ncol = K)
      Var_APE_CONT <- t(J) %*% V %*% J
    }
    
    APE <- APE_CONT
    if (variance_requested == 1) {
      Var_APE <- diag(Var_APE_CONT)
    }
    } else {
      APE<-beta;
      if (variance_requested==1){
      diag(Var_APE) = Var_beta;
      }
    }
    
    #vamos por aqui
    #CALCULATE BIAS CORRECTED ESTIMATORS VIA JACKKNIFE (assuming N=T for trade application)
    
    #if (nargsout < 10) {
    #if (missing(repMAX) && missing(MAX_STEPS)) {   #esto no es correcto, es si piden mas de 10 outputs, que de momento no lo hacemos
    if(jacknife==TRUE){
     if (N == T) {
        repMIN2 <- 1
        repMAX2 <- 1
        MAX_STEPS2 <- 150
        precision_beta2 <- precision_beta
        betaJK <- matrix(0, nrow = K, ncol = N)
        APEJK <- matrix(0, nrow = K, ncol = N)
        
        for (i in 1:N) {
          cat(paste0("Jackknife run ", i, " / ", N, "\n"))
          
          Y2 <- Y
          X2 <- X
          weight2 <- weight
          
          Y2 <- Y2[-i, ]
          Y2 <- Y2[,-i ]
          
          X2 <- X2[, -i, ] 
          X2 <- X2[, , -i]
          
          weight2 <- weight2[-i, ] 
          weight2 <- weight2[, -i]
          
          beta2 <- beta
          alpha2 <- alpha
          gamma2 <- gamma
          lambda2 <- lambda
          f2 <- f
          lambda_known2 <- lambda_known
          f_known2 <- f_known
          
          alpha2 <- alpha2[-i, ]
          gamma2 <- gamma2[-i, ] 
          lambda2 <- lambda2[-i, ]
          f2 <- f2[-i, ] 
          lambda_known2 <- lambda_known2[-i, ] 
          f_known2 <- f_known2[-i, ] 
          
          f_known2 <- matrix(f_known2,nrow=length(f_known2),ncol=1)
          lambda_known2 <- matrix(lambda_known2,nrow=length(lambda_known2),ncol=1)
          
          #lo siguiente es raro, es la funcion dentro de la funcion. Quizás el exotglag2 tenga algo que ver
          result <- FactorMLE(dist, Y2, X2, weight2, lambda_known2, f_known2, R, precision_beta2, repMIN2, repMAX2, MAX_STEPS2, beta2, alpha2, gamma2, lambda2, f2)
          beta2 <- result$beta
          APE2 <- result$APE
          gamma2 <- result$gamma
          lambda2 <- result$lambda
          f2 <-result$f
          obj2 <- result$obj
        
          
          betaJK[, i] <- beta2
          APEJK[, i] <- APE2
        }
        
        beta_corrected <- N * beta - (N - 1) * colMeans(betaJK)
        APE_corrected <- N * APE - (N - 1) * colMeans(APEJK)
      } else {
        beta_corrected <- rep(0, K)
        APE_corrected <- rep(0, K)
      }
    } else {
      beta_corrected <- rep(0, K)
      APE_corrected <- rep(0, K)
    }
     #no sacamos beta_corrected ni APE_corrected de momento, porque depende de no tenerlos en el input
    return(list(beta = beta, alpha = alpha, gamma = gamma,lambda=lambda,f=f,exitflag=exitflag,
                obj=obj,APE=APE,Var_beta=Var_beta,Var_APE=Var_APE,beta_corrected=beta_corrected,APE_corrected=APE_corrected,v_reference=v_reference))
    
       
}



##########Function step, que calcula un paso de iteracion similar al EM
mult <- function(X, beta) {
  mat <- array(0, dim = c(dim(X)[2], dim(X)[3]))
  for (k in 1:dim(X)[1]) {
    #mat <- mat + beta[k] * X[k, , , drop = FALSE]
    #mat <- mat + beta[k]*matrix(aperm(X, c(dim(mat)[3], 1, dim(X)[1]))[,,k], nrow = dim(mat)[2], ncol = dim(mat)[3])
    mat <- mat + beta[k]*X[k,,]
    }
  return(mat)
}
 
# Xs <- array(0,dim=c(3,4,2))
# 
# Xs[,,1] <- matrix(c(0.8147,0.9058,0.1270,0.9134,0.6324,0.0975,0.2785,0.5469,0.9575,0.9649,0.1576,0.9706),nrow=3,ncol=4)
# Xs[,,2] <- matrix(c(0.9572,0.4854,0.8003,0.1419,0.4218,0.9157,0.7922,0.9595,0.6557,0.0357,0.8491,0.9340),nrow=3,ncol=4)
# 
# betas <- c(1.5,2,0.5)
# 
# matrix(aperm(Xs, c(dim(Xs)[3], 1, dim(Xs)[1])), nrow = dim(Xs)[2], ncol = dim(Xs)[3])

model <- function(dist, y, z, weight) {
  if (dist == "logit") {
    # LOGIT LOG LIKELIHOOD AND DERIVATIVES:
    z <- pmax(pmin(z, 100), -100)  # truncate index at -100 and 100
    logL <- y * z - log(1 + exp(z))  # log-likelihood for each observation
    
    # if user also requests derivatives:
    if (length(names(match.call())) > 4) {
      dlogL <- y - exp(z) / (1 + exp(z))  # score wrt single index for each observation
      ddlogL <- -exp(z) / ((1 + exp(z))^2)  # Hessian wrt single index for each observation
      
    }
  } else if (dist == "probit") {
    # PROBIT LOG LIKELIHOOD AND DERIVATIVES:
    z <- pmax(pmin(z, 6), -6)  # truncate index at -6 and 6
    F <- pnorm(z)
    logL <- y * log(F) + (1 - y) * log(1 - F)  # log-likelihood for each observation
    
    if (length(names(match.call())) > 4) {
      F1 <- dnorm(z)
      dlogL <- y * F1 / F - (1 - y) * F1 / (1 - F)
      F2 <- -z * F1  # derivative of normal PDF
      ddlogL <- y * (F2 / F - F1^2 / F^2) - (1 - y) * (F2 / (1 - F) + F1^2 / (1 - F)^2)
      
    }
  } else if (dist == "Poisson") {
    # POISSON LOG LIKELIHOOD AND DERIVATIVES:
    logL <- y * z - exp(z)
    
    #if (length(names(match.call())) > 4) {
      dlogL <- y - exp(z)
      ddlogL <- -exp(z)
      
    #}
  } else {
    # MINUS LEAST SQUARES OBJECTIVE AND DERIVATIVES (i.e. STANDARD NORMAL LOG LIKELIHOOD):
    logL <- -(y - z)^2 / 2
    
    #if (length(names(match.call())) > 4) {
      dlogL <- y - z
      ddlogL <- -1 + 0 * y  # 0*y only to guarantee ddlogL has the same dimensions as y
      
    #}
  }
  
  # Weigh all observations:
  logL <- weight * logL
  
  #if (length(names(match.call())) > 4) {
    dlogL <- weight * dlogL
    ddlogL <- weight * ddlogL
    
    return(list(logL = logL, dlogL = dlogL, ddlogL = ddlogL))
  #}
  
  #return(logL)
}




step <- function(dist, Y, X, weight, lambda_known, f_known, beta0, alpha0, gamma0, lambda0, f0) {
  
  K <- dim(X)[1]
  N <- dim(X)[2]
  T <- dim(X)[3]
  R <- dim(lambda0)[2]
  
  fct <- 1  # if no progress, try to make a smaller step, normal step corresponds to fct = 1
  stop <- FALSE
  
  while (!stop) {
    # Definition of single index
    Z <- mult(X, beta0) + alpha0 %*% t(f_known) + lambda_known %*% t(gamma0) + lambda0 %*% t(f0)
    
    # Get corresponding log-likelihood and derivatives for each observation
    result <- model(dist, Y, Z, weight)
    logL <- result$logL
    dlogL <- result$dlogL
    ddlogL <- result$ddlogL
    obj0 <- mean(logL)
    
    # Add dlogL./wt to single index
    Z <- Z + fct * dlogL / mean(mean(-ddlogL))
    if(max(Z)>100){Z<- Z/100000}
    
    # Principal components step to get new alpha, gamma, lambda, f
    M_lambda_known <- diag(N) - lambda_known %*% solve(t(lambda_known) %*% lambda_known) %*% t(lambda_known)
    M_f_known <- diag(T) - f_known %*% solve(t(f_known) %*% f_known) %*% t(f_known)
    res <- M_lambda_known %*% (Z - mult(X, beta0)) %*% M_f_known
    
    result <- eigen(t(res) %*% res)
    if(R==0){
      f <- matrix(vector(), nrow = dim(X)[2], ncol = 0)
      
      #f <- matrix(0, nrow = dim(X)[2], ncol = 1)
    }else{
    f <- result$vectors[, 1:R]
    for (r in 1:R) {
      if(r==1){
        f <- as.matrix(result$vectors[, 1:R])
      }
      f[, r] <- f[, r] / norm(f[, r],type="2")
      if (mean(f[, r]) < 0) {
        f[, r] <- -f[, r]
      }
    }
    }
    lambda <- res %*% f
    f <- f * sqrt(T)
    lambda <- lambda / sqrt(T)
    
    res <- Z - mult(X, beta0) - lambda %*% t(f)
    #gamma <- solve(t(lambda_known) %*% lambda_known) %*% t(res)
    gamma <- t(res) %*% lambda_known / sum(lambda_known^2)
    res <- res - lambda_known %*% t(gamma)
    alpha <- res %*% f_known / sum(f_known^2)
    #alpha <- solve(t(f_known) %*% f_known) %*% t(res)
    
    # Update Z
    Z <- mult(X, beta0) + alpha %*% t(f_known) + lambda_known %*% t(gamma) + lambda %*% t(f)
    if(max(Z)>100){Z<- Z/100000}
    result <- model(dist, Y, Z, weight)                                                         #NOS QUEDAMOS POR AQUI, QUE HAY ERROR EN MODEL AL CREAR dlogL
    logL <- result$logL
    dlogL <- result$dlogL
    ddlogL <- result$ddlogL
    obj1 <- mean(logL)
    
    # WLS step to get new beta (possibility 1)
    Z <- Z + fct * dlogL / mean(mean(-ddlogL))
    if(all(f == 0)){
      lambda_all <- lambda_known
    } else {
    lambda_all <- cbind(lambda, lambda_known)
    }
    if(all(f == 0)){
      f_all <- f_known
    }else {
    f_all <- cbind(f, f_known)
    }
    Mlambda_all <- diag(N) - lambda_all %*% solve(t(lambda_all) %*% lambda_all) %*% t(lambda_all)
    Mf_all <- diag(T) - f_all %*% solve(t(f_all) %*% f_all) %*% t(f_all)
    W <- matrix(0, nrow = K, ncol = K)
    V <- matrix(0, nrow = K, ncol = 1)
    for (k1 in 1:K) {
      Xk1 <- Mlambda_all %*% X[k1, , ] %*% Mf_all    #no coinciden
      V[k1] <- 1 / (N * T) * sum(diag(Xk1 %*% t(Z)))
      for (k2 in 1:K) {
        Xk2 <- Mlambda_all %*% X[k2, , ] %*% Mf_all
        W[k1, k2] <- 1 / (N * T) * sum(diag(Xk1 %*% t(Xk2)))
      }
    }
    beta <- solve(W, V)
    
    # WLS step to get new beta (possibility 2)
    # Z <- mult(X, beta0) + fct * dlogL / pmax(-ddlogL, 1e-7)
    # W <- matrix(0, nrow = K, ncol = K)
    # V <- matrix(0, nrow = K, ncol = 1)
    # for (k1 in 1:K) {
    #   Xk1 <- X[k1, , ] * (-ddlogL)
    #   V[k1] <- 1 / (N * T) * sum(diag(Xk1 %*% t(Z)))
    #   for (k2 in 1:K) {
    #     Xk2 <- X[k2, , ]
    #     W[k1, k2] <- 1 / (N * T) * sum(diag(Xk1 %*% t(Xk2)))
    #   }
    # }
    # beta <- solve(W, V)

    obj2 <- mean(model(dist, Y, mult(X, beta) + alpha %*% t(f_known) + lambda_known %*% t(gamma) + lambda %*% t(f), weight)$logL)
    
    if (obj2 > obj0) {
      stop <- TRUE  # Stop if we indeed increased the likelihood function in this step
    } else {
      fct <- fct / 2  # Reduce the "steplength" if no increase in likelihood
      if (fct < 1e-4) {
        stop <- TRUE  # Final stopping criterion
      }
    }
  }
  
  if (obj2 > obj0) {
    # If progress was made in likelihood maximization
    obj <- obj2  # Report the new objective and updated parameters
  } else {
    # Otherwise, report the old parameters and objective
    beta <- beta0
    alpha <- alpha0
    gamma <- gamma0
    f <- f0
    lambda <- lambda0
    obj <- obj0
  }
  
  return(list(beta = beta, alpha = alpha, gamma = gamma, lambda = lambda, f = f, obj = obj))
}


#Function stepNR, que calcular un paso de iteracion con Newton-Raphson


stepNR <- function(dist, Y, X, weight, lambda_known, f_known, beta0, alpha0, gamma0, lambda0, f0, v_reference, para_reference, b_reference) {
  K <- dim(X)[1]
  N <- dim(X)[2]
  T <- dim(X)[3]
  R <- dim(lambda0)[2]
  Rex1 <- dim(lambda_known)[2]
  Rex2 <- dim(f_known)[2]
  
  # get log-likelihood, score, Hessian for each observation
  obj_score_hessian <- SampleLogL(dist, Y, X, weight, lambda_known, f_known, beta0, alpha0, gamma0, lambda0, f0, 1)
  obj0 <- obj_score_hessian$obj
  score0 <- obj_score_hessian$s
  Hessian0 <- obj_score_hessian$H
  
  # rescale score and Hessian to account for parameter rescaling
  rescale <- c(rep(1, K), rep(sqrt(N), N * (R + Rex2)), rep(sqrt(T), T * (R + Rex1)))
  score0 <- rescale * score0
  Hessian0 <- outer(rescale, rescale) * Hessian0
  
  # check if v_reference is still appropriate
  if (b_reference > 0) {
    v <- flat_directions(alpha0, gamma0, lambda0, f0, lambda_known, f_known, K)
    ev_min <- min(abs(eigen(t(v) %*% v_reference %*%t(v_reference)%*% v)$values))   #cuidado aqui
  } else {
    ev_min <- 0
  }
  
  if (ev_min < 0.1) {   # criterion for updating reference parameter
    cat("updating reference parameter for penalized objective function\n")
    alpha0 <- normalize(alpha0, gamma0, lambda0, f0, lambda_known, f_known)[["alpha"]]
    gamma0 <- normalize(alpha0, gamma0, lambda0, f0, lambda_known, f_known)[["gamma"]]
    lambda0 <- normalize(alpha0, gamma0, lambda0, f0, lambda_known, f_known)[["lambda"]]
    f0 <- normalize(alpha0, gamma0, lambda0, f0, lambda_known, f_known)[["f"]]
    para_reference <- para_transform2(beta0, alpha0, gamma0, lambda0, f0)
    v_reference <- flat_directions(alpha0, gamma0, lambda0, f0, lambda_known, f_known, K)
    b_reference <- abs(pracma::Trace(-Hessian0) / nrow(Hessian0))
  }
  
  # penalized objective
  para0 <- para_transform2(beta0, alpha0, gamma0, lambda0, f0)
  objP <- obj0 - b_reference / 2 * crossprod((para0 - para_reference), v_reference) %*% t(v_reference) %*% (para0 - para_reference)
  score <- score0 - b_reference * v_reference %*% (t(v_reference) %*% (para0 - para_reference))
  Hessian <- Hessian0 - b_reference * v_reference %*% t(v_reference)
  
  epsilon_max <- 1
  epsilon_min <- 10^(-4)
  epsilon <- 10^(-5)
  
  stop2 <- 0
  while (stop2==0) {
    #direction <- solve(-Hessian + b_reference * epsilon * diag(nrow(Hessian)), score)
    direction <- pracma::mldivide(-Hessian + b_reference * epsilon * diag(nrow(Hessian)),score)  ############ESTO DIFIERE DE MATLAB
    direction <- direction/sqrt(sum(direction^2))
    #step_para <- direction %*% (t(direction) %*% score) / (t(direction) %*% (-Hessian) %*% direction)
    step_para <- mrdivide(direction %*% (t(direction) %*% score),(t(direction) %*% (-Hessian) %*% direction))
    
    # translate step_para into original parameters
    step_beta <- NULL
    step_alpha <- NULL
    step_gamma <- NULL
    step_lambda <- NULL
    step_f <- NULL
    step_transform <- para_transform1(step_para, K, N, T, R, Rex1, Rex2)
    step_beta <- step_transform$beta
    step_alpha <- step_transform$alpha
    step_gamma <- step_transform$gamma
    step_lambda <- step_transform$lambda
    step_f <- step_transform$f
    
    # try to update parameters, but making sure that penalized objective increases
    fct <- 1   # if we cannot make progress, then we try to make "a smaller step",
    # but normal step corresponds to fct=1
    stop <- 0
    while (stop==0) {
      para <- para0 + fct * step_para
      beta <- beta0 + fct * step_beta
      lambda <- lambda0 + fct * step_lambda
      alpha <- alpha0 + fct * step_alpha
      f <- f0 + fct * step_f
      gamma <- gamma0 + fct * step_gamma
      
      obj <- SampleLogL(dist, Y, X, weight, lambda_known, f_known, beta, alpha, gamma, lambda, f)$obj
      objPnew <- obj - b_reference / 2 * crossprod((para - para_reference), v_reference) %*% t(v_reference) %*% (para - para_reference)
      
      if (objPnew > objP) {
        stop <- 1
        stop2 <- 1
      } else {
        if (fct < 1/2 && epsilon < epsilon_max * 0.9999) {
          fct <- 1
          epsilon <- max(epsilon_min, epsilon * 10)
          stop <- 1
        } else {
          fct <- fct / 2   # if we cannot make progress, then we reduce the "steplength"
          if (fct < 10^(-4))   # need some final stopping criterion
            stop <- 1
          stop2 <- 1
        }
      }
    }
  }
  
  if (objPnew < objP) {   # if no progress was made in the likelihood maximization
    beta <- beta0   # then just report the old parameters and objective
    alpha <- alpha0
    gamma <- gamma0
    f <- f0
    lambda <- lambda0
    obj <- obj0
  }
  
  cat(fct, "\n")
  cat(epsilon, "\n")
  
  return(list(beta = beta, alpha = alpha, gamma = gamma, lambda = lambda, f = f, obj = obj, v_reference = v_reference, para_reference = para_reference, b_reference = b_reference))
}



SampleLogL <- function(dist, Y, X, weight, lambda_known, f_known, beta, alpha, gamma, lambda, f, Hlin=0, logL=NULL, dlogL=NULL, ddlogL=NULL) {
#SampleLogL <- function(dist, Y, X, weight, lambda_known, f_known, beta, alpha, gamma, lambda, f, Hlin = TRUE, logL = NULL, dlogL = NULL, ddlogL = NULL) {
  # score and Hessian use the following order of parameters:
  # beta, lambda, alpha, f, gamma
  
  K <- dim(X)[1]
  N <- dim(X)[2]
  T <- dim(X)[3]
  R <- dim(lambda)[2]
  Rex1 <- dim(lambda_known)[2]
  Rex2 <- dim(f_known)[2]
  
  ff <- cbind(f, f_known)
  ll <- cbind(lambda, lambda_known)
  
  #if (missing(logL) | missing(dlogL) | missing(ddlogL)) {
  if(is.null(logL)|is.null(dlogL)|is.null(ddlogL)){
    # if logL, dlogL, ddlogL not already provided as inputs, then compute them
    Z <- mult(X,beta) + alpha %*% t(f_known) + lambda_known %*% t(gamma) + lambda %*% t(f)
    if(max(Z)>100){Z<- Z/100000}
    logL <- model(dist, Y, Z, weight)$logL
    dlogL <- model(dist, Y, Z, weight)$dlogL
    ddlogL <- model(dist, Y, Z, weight)$ddlogL
  }
  
  # calculate log-likelihood:
  obj <- sum(logL) / (N * T)
  
  
  # calculate score:
  s <- matrix(0, nrow = K + N * (R + Rex2) + T * (R + Rex1), ncol = 1)
  for (k in 1:K) {
    XX <- X[k, ,]
    s[k] <- sum(dlogL * XX) / (N * T)
  }
  for (r in 1:(R + Rex2)) {
    s[K + (r - 1) * N + (1:N)] <- dlogL %*% ff[, r] / (N * T)
  }
  for (r in 1:(R + Rex1)) {
    s[K + N * (R + Rex2) + (r - 1) * T + (1:T)] <- t(dlogL) %*% ll[, r] / (N * T)
  }

  
  # calculate Hessian:
  H <- matrix(0, nrow = K + N * (R + Rex2) + T * (R + Rex1), ncol = K + N * (R + Rex2) + T * (R + Rex1))
  for (k1 in 1:K) {
    X1 <- X[k1, , ]
    for (k2 in 1:K) {
      X2 <- X[k2, , ]
      H[k1, k2] <- sum(ddlogL * X1 * X2) / (N * T)
    }
    for (r in 1:(R + Rex2)) {
      ind <- K + (r - 1) * N + (1:N)
      H[k1, ind] <- (X1 * ddlogL) %*% ff[, r] / (N * T)
      H[ind, k1] <- t(H[k1, ind])
    }
    for (r in 1:(R + Rex1)) {
      ind <- K + N * (R + Rex2) + (r - 1) * T + (1:T)
      H[k1, ind] <- (X1 * ddlogL) %*% ll[, r] / (N * T)
      H[ind, k1] <- t(H[k1, ind])
    }
  }
  
  for (r1 in 1:(R + Rex2)) {
    ind1 <- K + (r1 - 1) * N + (1:N)
    for (r2 in 1:(R + Rex2)) {
      ind2 <- K + (r2 - 1) * N + (1:N)
      diag(H[ind1, ind2]) <-ddlogL %*% (ff[, r1] * ff[, r2]) / (N * T)
      #H[ind1, ind2] <- diag(ddlogL %*% (ff[, r1] * ff[, r2]) / (N * T)) 
    }
  }
  
  for (r1 in 1:(R + Rex1)) {
    ind1 <- K + N * (R + Rex2) + (r1 - 1) * T + (1:T)
    for (r2 in 1:(R + Rex1)) {
      ind2 <- K + N * (R + Rex2) + (r2 - 1) * T + (1:T)
      diag(H[ind1, ind2]) <- t(ddlogL) %*% (ll[, r1] * ll[, r2]) / (N * T)
      #H[ind1, ind2] <- diag(ddlogL %*% (ll[, r1] * ll[, r2]) / (N * T))
    }
  }
  
  #aqui
  
  for (r1 in 1:(R + Rex2)) {
    ind1 <- K + (r1 - 1) * N + (1:N)
    for (r2 in 1:(R + Rex1)) {
      ind2 <- K + N * (R + Rex2) + (r2 - 1) * T + (1:T)
      H[ind1, ind2] <- ddlogL * (ll[, r2] %*% t(ff[, r1])) / (N * T)
      H[ind2, ind1] <- t(H[ind1, ind2])
    }
  }
  
  if (Hlin==1 && R!=0) {
    # add terms in the Hessian, which are only there because index itself is non-linear.
    # those terms depend on the first derivative of the log-likelihood (dlogL)
    for (r in 1:R) {
      ind1 <- K + (r - 1) * N + (1:N)
      ind2 <- K + N * (R + Rex2) + (r - 1) * T + (1:T)
      H[ind1, ind2] <- H[ind1, ind2] + dlogL / (N * T)
      H[ind2, ind1] <- H[ind2, ind1] + t(dlogL) / (N * T)
    }
  }
  
  return(list(obj = obj, s = s, H = H))
}



########################33
#translate parameter vector of length K + N*(R+Rex2) + T*(R+Rex1) back into
#the original parameters:
para_transform1 <- function(para, K, N, T, R, Rex1, Rex2) {
  beta <- para[1:K]
  lambda <- matrix(0, nrow = N, ncol = R)
  f <- matrix(0, nrow = T, ncol = R)
  alpha <- matrix(0, nrow = N, ncol = Rex2)
  gamma <- matrix(0, nrow = T, ncol = Rex1)
  
  if(R!=0){
  for (r in 1:R) {
    lambda[, r] <- para[K + (r - 1) * N + (1:N)] * sqrt(N)
  }
  }
  for (r in 1:Rex2) {
    alpha[, r] <- para[K + (R + r - 1) * N + (1:N)] * sqrt(N)
  }
  if(R!=0){
  for (r in 1:R) {
    f[, r] <- para[K + N * (R + Rex2) + (r - 1) * T + (1:T)] * sqrt(T)
  }
  }
  for (r in 1:Rex1) {
    gamma[, r] <- para[K + N * (R + Rex2) + (R + r - 1) * T + (1:T)] * sqrt(T)
  }
  
  return(list(beta = beta, alpha = alpha, gamma = gamma, lambda = lambda, f = f))
}



##############################
#translate original parameters into single parameter vector of length 
#K + N*(R+Rex2) + T*(R+Rex1)

para_transform2 <- function(beta, alpha, gamma, lambda, f){
  K <- length(beta)
  N <- nrow(lambda)
  T <- nrow(f)
  R <- ncol(lambda)
  Rex1 <- ncol(gamma)
  Rex2 <- ncol(alpha)
  
  para <- rep(0, K + N * (R + Rex2) + T * (R + Rex1))
  
  para[1:K] <- beta
  
  if(R!=0){
  for (r in 1:R) {
    para[K + (r - 1) * N + (1:N)] <- lambda[, r] / sqrt(N)
  }
  }
  for (r in 1:Rex2) {
    para[K + (R + r - 1) * N + (1:N)] <- alpha[, r] / sqrt(N)
  }
  if(R!=0){
  for (r in 1:R) {
    para[K + N * (R + Rex2) + (r - 1) * T + (1:T)] <- f[, r] / sqrt(T)
  }
  }
  for (r in 1:Rex1) {
    para[K + N * (R + Rex2) + (R + r - 1) * T + (1:T)] <- gamma[, r] / sqrt(T)
  }
  
  return(para)
}




#######################################
#define the matrix of flat directions of the likelihood:

flat_directions <- function(alpha, gamma, lambda, f, lambda_known, f_known, K) {
  N <- nrow(lambda)
  T <- nrow(f)
  R <- ncol(lambda)
  Rex1 <- ncol(lambda_known)
  Rex2 <- ncol(f_known)
  
  ff <- cbind(f, f_known)
  ll <- cbind(lambda, lambda_known)
  
  v <- matrix(0, nrow = K + N * (R + Rex2) + T * (R + Rex1), ncol = (R + Rex2) * (R + Rex1))
  cnt <- 0
  
  for (r1 in 1:(R + Rex2)) {
    for (r2 in 1:(R + Rex1)) {
      ind1 <- K + (r1 - 1) * N + (1:N)
      ind2 <- K + N * (R + Rex2) + (r2 - 1) * T + (1:T)
      vv <- matrix(0, K + N * (R + Rex2) + T * (R + Rex1),1)
      vv[ind1] <- ll[, r2]
      vv[ind2] <- -ff[, r1]
      vv <- vv / sqrt(sum(vv^2))   #en lugar de norm
      cnt <- cnt + 1
      v[, cnt] <- vv
    }
  }
  
  return(v)
}



#####################################################
#normalize parameters
normalize <- function(alpha, gamma, lambda, f, lambda_known, f_known) {
  N <- nrow(lambda)
  T <- nrow(f)
  
  M_lambda_known <- diag(N) - lambda_known %*% solve(t(lambda_known) %*% lambda_known) %*% t(lambda_known)
  M_f_known <- diag(T) - f_known %*% solve(t(f_known) %*% f_known) %*% t(f_known)
  
  # Find normalized alpha and gamma:
  Z <- alpha %*% t(f_known) + lambda_known %*% t(gamma) + lambda %*% t(f)
  Z <- Z - M_lambda_known %*% Z %*% M_f_known  # Part of index that can be explained by lambda_known and f_known
  #gamma <- solve(t(lambda_known) %*% lambda_known) %*% t(lambda_known) %*% Z
  gamma <- t(Z)%*%lambda_known%*% solve(t(lambda_known)%*%lambda_known)
  Z <- Z - lambda_known %*% t(gamma)
  alpha <- Z %*% f_known / sum(f_known^2)
  
  # Find normalized lambda and f:
  f <- M_f_known %*% f
  lambda <- M_lambda_known %*% lambda  # We already took care of everything that can be explained by f_known and lambda_known
  
  A <- pracma::sqrtm(t(f) %*% f / T)$B  # Note that this is symmetric
  if(dim(A)[1]==0 && dim(A)[2]==0){
    f <- matrix(vector(), nrow = T, ncol = dim(A)[2])
  }else{
  f <- f / A[[1]]
  }
  lambda <- lambda %*% A
  
  return(list(alpha = alpha, gamma = gamma, lambda = lambda, f = f))
}

