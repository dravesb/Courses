library(igraph)
library(pscl)

TAU2 <- 100 # default prior variance
build.null.prior <- function (d, tau2=TAU2)
  list(mean = numeric(d), Sigma.inv = rep(1 / tau2, d))

lprior <- function (beta, prior) {
  if (is.vector(prior$Sigma.inv))
    return(-.5 * sum((beta - prior$mean) ^ 2 * prior$Sigma.inv))
  C.inv <- chol(prior$Sigma.inv)
  -.5 * crossprod(C.inv %*% (beta - prior$mean))
}

# supported exponential families:
exp.families <- c("binomial", "poisson", "gaussian")

# `g` is graph, `vars` is a list with each predictor data, `factors` is a list with
# each factor predictor data, `varfuns` and `facfuns` are vectors specifying
# the kernel function: 0 (FALSE) for '*' or 1 (TRUE) for '+', for each
# predictor, `response` is the graph attribute for the response variable, and
# `family` is the exponential family to be used.
gglmParamSet <- setRefClass("gglmParamSet",
  fields = c(".pset", "graph", "family.index", "nlevel", "slevel",
             "varfun", "facfun", "vardata", "facdata", "response", "size"),
  methods = list(
    # TODO: check if family in c('binomial',...)
    # 'size' is valid only for 'binomial' family (usual interpretation)
    initialize = function (g, vars=NULL, factors=NULL,
                           varfuns=NULL, facfuns=NULL,
                           response=NULL, family="binomial", size=0) {
      nvar <- length(vars)
      nfac <- length(factors)
      if (nvar == 0 && nfac == 0) stop("no variables to fit")
      graph <<- g
      fam.id <- match(family, exp.families)
      if (is.na(fam.id))
        stop(paste("family '", family, "' is not supported"))
      if (!is.null(response) && !(response %in% list.edge.attributes(g)))
        stop("'", response, "' is not a valid edge attribute")

      family.index <<- fam.id - 1L # zero-based
      if (length(factors)==0) {  #Liz changes this so factors can be an empty list.  
        nlevel <<- rep(0L, 0)
        slevel <<- 0L
      }
      else {
        nlevel <<- sapply(factors, function (f) length(levels(f)))
        slevel <<- sum(nlevel)
      }
      if (is.null(varfuns)) varfuns <- rep(0L, nvar)
      if (is.null(facfuns)) facfuns <- rep(0L, nfac)
      varfun <<- as.integer(varfuns)
      facfun <<- as.integer(facfuns)
      n <- vcount(graph)
      vardata <<- matrix(0, nrow=n, ncol=nvar)
      facdata <<- matrix(0L, nrow=n, ncol=nfac)
      if (nvar > 0)
        for (i in 1:nvar) vardata[,i] <<- vars[[i]]
      if (nfac > 0) {
        for (i in 1:nfac)
          facdata[,i] <<- as.integer(factors[[i]]) - 1L # zero-based
      }
      response <<- as.character(response)
      size <<- as.numeric(size)

      .pset <<- .Call(G_new_gglmparamset, graph, nvar, nfac, family.index,
                      nlevel, slevel, varfun, facfun, vardata, facdata,
                      response, size, PACKAGE="gglm")
    },

    reset = function () {
      nvar <- ncol(vardata)
      nfac <- ncol(facdata)
      .pset <<- .Call(G_new_gglmparamset, graph, nvar, nfac, family.index,
                      nlevel, slevel, varfun, facfun, vardata, facdata,
                      response, size, PACKAGE="gglm")
    },

    setfactor = function (fac, z) {
      if (fac < 1 || fac > ncol(facdata)) stop("invalid factor index")
      if (length(z) != vcount(graph)) stop("invalid vector length")
      .Call(G_setfactor_gglmparamset, .pset, as.integer(fac - 1),
            as.integer(z), PACKAGE="gglm")
    },

    setvariable = function (var, z) {
      if (var < 1 || var > ncol(vardata)) stop("invalid variable index")
      if (length(z) != vcount(graph)) stop("invalid vector length")
      .Call(G_setvariable_gglmparamset, .pset, as.integer(var - 1),
            as.double(z), PACKAGE="gglm")
    },


    predict = function (beta, offset = NULL, type = "link") {
      n <- vcount(graph)
      r <- match(type, c("link", "response", "variance", "cumulant"))
      if (is.na(r))
        stop("invalid type")
      r <- as.integer(r - 1)
      matrix(.Call(G_predict_gglmparamset, .pset, as.numeric(beta),
                   as.numeric(offset), r, numeric(n ^ 2), PACKAGE="gglm"),
             nrow=n, ncol=n)
    },

    # 'stats' computes sufficient statistics for 'fit':
    # st$m = X' * (y - mu(beta)) and st$V = X' * diag(W(beta)) * X
    # if 'initialize', st$m = X' * (y - mu) + X' * diag(W(mu)) * eta(mu)
    stats = function (beta, offset = NULL, initialize = FALSE) {
      p <- ncol(vardata) + sum(nlevel)
      m <- numeric(p)
      V <- numeric(p ^ 2)
      # write m and V in-place:
      .Call(G_stats_pairs_gglmparamset, .pset, as.numeric(beta),
            as.numeric(offset), as.integer(initialize),
            m, V, PACKAGE="gglm")
      list(m=m, V=matrix(V, nrow=p))
    },

    lhood = function (beta, offset=NULL) {
      .Call(G_lhood_gglmparamset, .pset, as.numeric(beta), as.numeric(offset),
            PACKAGE="gglm")
    },


    fit = function (start = NULL, prior = NULL, offset = NULL, initialize = FALSE,
                    epsilon = 1e-8, maxit = 25, trace = FALSE) {
      d <- ncol(vardata) + sum(nlevel)
      if (!is.null(prior) && length(prior$mean) != d)
        stop("invalid 'prior' dimension")
      if (!is.null(start))
        beta <- start
      else
        beta <- numeric(d)

      lhood.cur <- -Inf
      it <- ifelse(initialize, 0, 1)
      repeat {
        st <- stats(beta, offset, initialize)
        if (!is.null(prior)) {
          if (is.vector(prior$Sigma.inv)) {
            if (initialize)
              st$m <- st$m + prior$Sigma.inv * prior$mean
            else
              st$m <- st$m - prior$Sigma.inv * (beta - prior$mean)
            diag(st$V) <- diag(st$V) + prior$Sigma.inv
          }
          else {
            if (initialize)
              st$m <- st$m + drop(prior$Sigma.inv %*% prior$mean)
            else
              st$m <- st$m - drop(prior$Sigma.inv %*% (beta - prior$mean))
            st$V <- st$V + prior$Sigma.inv
          }
        }

        C <- chol(st$V)
        if (initialize)
          beta.new <- backsolve(C, backsolve(C, st$m, transpose=TRUE))
        else
          beta.new <- beta + backsolve(C, backsolve(C, st$m, transpose=TRUE))

        lhood.new <- lhood(beta.new, offset)
        if (!is.null(prior))
          lhood.new <- lhood.new + lprior(beta.new, prior)
        if (is.infinite(lhood.new)) warning("log likelihood is infinity")
        if (!initialize && !is.infinite(lhood.new) &&
            (lhood.new < lhood.cur ||
             abs((lhood.new - lhood.cur) / lhood.new) < epsilon))
          break
        beta <- beta.new; lhood.cur <- lhood.new
        if (trace) message("[", it, "] lhood = ", lhood.cur)
        if (initialize || it >= maxit) break
        it <- it + 1
      }
      list(coef = beta, C = C, lhood = lhood.cur)
    }
  )
)


gglm.fit <- function (g, vars = NULL, factors = NULL,
                      varfuns = NULL, facfuns = NULL,
                      prior = NULL, # list(mean, Sigma.inv) or function (beta)
                      start = NULL, offset = NULL,
                      epsilon = 1e-8, maxit = 25, trace = FALSE, ...) {
  gp <- gglmParamSet(g, vars, factors, varfuns, facfuns, ...)
  d <- ncol(gp$vardata) + sum(gp$nlevel)

  fprior <- prior
  if (is.function(prior)) fprior <- prior(start)
  if (is.null(start)) { # initialize?
    res <- gp$fit(start, fprior, offset, initialize = TRUE, trace = trace)
    beta <- res$coef; lhood <- res$lhood;
  }
  else {
    beta <- start; lhood <- -Inf
  }
  if (!is.null(fprior))
    lhood <- lhood + lprior(beta, fprior) # log joint

  it <- 1
  repeat {
    if (is.function(prior)) fprior <- prior(beta)
    res <- gp$fit(beta, fprior, offset, maxit = 1)
    beta.new <- res$coef; lhood.new <- res$lhood
    if (!is.null(fprior))
      lhood.new <- lhood.new + lprior(beta.new, fprior) # log joint

    if (lhood.new < lhood || abs((lhood.new - lhood) / lhood.new) < epsilon)
      break
    beta <- beta.new; lhood <- lhood.new; C <- res$C
    if (trace) message("[", it, "] lhood = ", lhood)
    if (it >= maxit) break
    it <- it + 1
  }
  list(coef = beta, C = C, lhood = lhood, gp = gp)
}

#added by Liz.  gglm object and functions node, factor.node, and ed.a
#if listing response, do so like response = "namehere".  It compares with list.edge.attributs.  

gglm <- function(formula, response = NULL, family = c("binomial"),  epsilon = 1e-8, maxit = 25, offset = NULL, prior = NULL, start = NULL, ...) {
  
  call <- match.call() #returns input
  terms <- terms(formula)
  g <-try(eval(terms[[2]], attr(terms,".Environment"))) #we are evaluating terms[[2]] in the environment of the formula declaration.  Thus HERE we set g = t if formula was t ~ ...  
  
  if(inherits(g,"try-error")){
    stop("Invalid formula")
  }
  
  ####set up covariates based on formula###
  termlist <- as.list(attr(terms,"variables"))[-(1:2)]  #this does not include interactions. termlist will be empty if only intercept.
  
  ###checking intercept first.  It will be first node attribute.###
  vars <- list()
  varfuns <- c()
  factors <- list()
  facfuns <- c()
  edgep <- list()
  
  if (attr(terms, "intercept") == 1) {
    vars[[1]] <- 1
    varfuns[1] <- TRUE #*
  }
  
  ###evaluating each term in formula and put into correct list and corresponding vector
  ###we have vars, varfuns, factors, facfuns, MUST ADD edge attributes
  
  m <- c("node","node.factor","ed.a")
  for (i in seq_along(termlist))   { 
    t <- charmatch(termlist[[i]][1], m, nomatch = 4)
    
    if (t == 1){
      vars <- c(vars, list(eval(termlist[[i]], attr(terms,".Environment"))$vec))
      varfuns <- c(varfuns, eval(termlist[[i]], attr(terms,".Environment"))$fun)
    }else{
      if (t == 2) {
        factors <- c(factors, list(eval(termlist[[i]], attr(terms,".Environment"))$vec)) 
        facfuns <- c(facfuns, eval(termlist[[i]], attr(terms,".Environment"))$fun)
      }else{
        if (t == 3) {
          edgep <- c(edgep, list(eval(termlist[[i]],attr(terms,".Environment")))) #TBD
        }else{ 
          print(termlist[[i]])
          stop("predictor not recognized") 
        }
      }
    }
  }
  
  gglmfit <- gglm.fit(g, response = response, vars = vars, factors = factors, 
                      varfuns = varfuns, facfuns = facfuns, 
                      prior = prior, start = start, 
                      offset = offset, epsilon = epsilon, 
                      maxit = maxit, trace = FALSE, family = family, ...)
  
  
  #changing output, this is a work in progress!  
  out <- list(coefficients = gglmfit$coef, 
              residuals = 1, #(y-mu)/mu.eta(eta)
              fitted.values = gglmfit$gp$predict(gglmfit$coef, type = 'response'),
              effects = 0,
              rank = 1,
              qr = 1,
              #family = family(),
              linear.predictors = 1, #eta
              deviance = 1, #family$deviance
              null.deviance = 1,
              aic = 1, #family$aic
              iter = 1, #eventually gglmfit$it
              weights = 1,
              prior.weights = 1,
              df.residual = 1,
              df.null = 1,
              y = g,
              converged = 1,
              boundary = 1, 
              call = call,
              formula = formula, 
              terms = terms,
              data = environment(formula),
              offset = offset,
              control = 1,
              method = 'gglm.fit',
              contrasts = 1,
              xlevels = 1)
  
  class(out) <- c("glm","lm")
  out

  #TO DO:
  #incorporate prior
  #incorporate edge predictor
  #return all of the appropriate gglm pieces in out (and then rename accordingly)
  #incorporate summary(out)
}

node.factor <- function(x, fun = c('+','*', 'i', 'j'))  { 
  
  x <- as.factor(x)
  
  x.na <- is.na(x)
  if (any(x.na)) 
    stop("The vertex attribute argument 'x' contains missing values.")
  
  if (missing(fun))
    fun <- c('+')
  
  preds <- list()
  preds$vec <- x
  
  preds$fun <- switch(fun, 
                      '+' = FALSE, 
                      '*' = TRUE, 
                      'i' = 3, #to fix
                      'j' = 4) #to fix
  return(preds)
  
  #TO DO:
  #add prior
  #add constraint
  #add different functionality
}

node <- function(x, fun = c('+','*', 'i', 'j'))  { 
  
  if (!is.numeric(x)) 
    stop("The vertex attribute argument 'x' must be a numeric vertex attribute.")
  
  x.na <- is.na(x)
  if (any(x.na)) 
    stop("The vertex attribute argument 'x' contains missing values.")
  
  if (missing(fun))
    fun <- c('+')
  
  preds <- list()
  preds$vec <- x
  
  preds$fun <- switch(fun, 
                      '+' = FALSE, 
                      '*' = TRUE, 
                      'i' = 3, #to fix
                      'j' = 4) #to fix
  return(preds)
}

ed.a <- function(x)  {  
    
    #x is an edge attribute.  
    if (!is.numeric(x)) 
      stop("The edge attribute 'x' must be a numeric edge attribute.")
    
    x.na <- is.na(x)
    if (any(x.na)) 
      stop("The edge attribute 'x' contains missing values.")
    
    preds <- x
    return(preds)
  }
  
  #function used to organize the model parts.  It creates
  #lists for the zero (Y/N) counts and the count model (poisson).
  
  getpreds <- function(terms = terms, formula) {
  
  terms2 <- terms(formula)
  termlist <- as.list(attr(terms2,"variables"))[-(1:2)]  #this does not include interactions. termlist will be empty if only intercept.
  
  ###checking intercept first.  It will be first node attribute.###
  vars <- list()
  varfuns <- c()
  factors <- list()
  facfuns <- c()
  edgep <- list()
  
  if (attr(terms2, "intercept") == 1) {
    vars[[1]] <- 1
    varfuns[1] <- TRUE 
  }
  
  ###evaluating each term in formula and put into correct list and corresponding vector
  ###we have vars, varfuns, factors, facfuns, MUST ADD edge attributes
  
  m <- c("node","node.factor","ed.a")
  for (i in seq_along(termlist))   { 
    t <- charmatch(termlist[[i]][1], m, nomatch = 4)
    
    if (t == 1){
      vars <- c(vars, list(eval(termlist[[i]], attr(terms,".Environment"))$vec))
      varfuns <- c(varfuns, eval(termlist[[i]], attr(terms,".Environment"))$fun)
    }else{
      if (t == 2) {
        factors <- c(factors, list(eval(termlist[[i]], attr(terms,".Environment"))$vec)) 
        facfuns <- c(facfuns, eval(termlist[[i]], attr(terms,".Environment"))$fun)
      }else{
        if (t == 3) {
          edgep <- c(edgep, list(eval(termlist[[i]],attr(terms,".Environment")))) #TBD
        }else{ 
          print(termlist[[i]])
          stop("predictor not recognized") 
        }
      }
    }
  }
  
  list(vars = vars, varfuns = varfuns, factors = factors, facfuns = facfuns)
}

#####the ZIP version###
gglm2 <- function(formula, response = NULL, family = c("binomial"),  epsilon = 1e-8, maxit = 25, offset = NULL, prior = NULL, start = NULL, ...) {
  #do we need two separate offsets?  
  
  cl <- match.call()  
  terms <- terms(formula)
  
  g <-try(eval(terms[[2]], attr(terms,".Environment"))) #we are evaluating terms[[2]] in the environment of the formula declaration.  Thus HERE we set g = t if formula was t ~ ...  
  
  if(inherits(g,"try-error")){
    stop("Response must be a network.")
  }
  
  #if we are doing a zero inflated model.  Floating integers just pick correct piece of formula.   
  if (length(formula[[3]]) > 1 && identical(formula[[3]][[1]],  #if we are doing a zero inflated model!  
                                            as.name("|"))) {
    ff <- formula #full formula
    
    ffc <- . ~ . #count formula (will have to update the response to only include non-zeros at some cut off
    ffz <- . ~ . #zero inflated formula
    
    ffc[[2]] <- ff[[2]]
    ffc[[3]] <- ff[[3]][[2]]
    
    ffz[[2]] <- ff[[2]]
    ffz[[3]] <- ff[[3]][[3]]
    
    predsz <- getpreds(terms = terms, formula = ffz)
    predsc <- getpreds(terms = terms, formula = ffc)

  }else{
  
    preds <- getpreds(terms = terms, formula = formula)
    }
  
  if(exists("preds")) { 
    
    gglmfit <- gglm.fit(g, response = response, vars = preds$vars, factors = preds$factors, 
                        varfuns = preds$varfuns, facfuns = preds$facfuns, 
                        prior = prior, start = start, 
                        offset = offset, epsilon = epsilon, 
                        maxit = maxit, trace = FALSE, family = family, ...)
  }
  
  to.upper <- function(X) t(X)[lower.tri(X)]
  logit <- function (p) log(p / (1 - p))
  
  if(exists("predsz")) {
  #currently, if doing  a zip, you do not supply a family.  It is taken to be
  #poisson for count and binomial for zeros
  #note, pis = # when no edge.  pis = 0 when edge.  
  
  
    #initial estimates
	cgf <- gglm.fit(g, response = response, vars = predsc$vars, factors = predsc$factors, 
                        varfuns = predsc$varfuns, facfuns = predsc$facfuns, 
                        prior = prior, start = start, 
                        offset = offset, epsilon = epsilon, 
                        maxit = maxit, trace = FALSE, family = c("poisson"), ...)
    betas <- cgf$coef
    betas <- rep(1, length = length(cgf$coef)) #testing start value
    
    cfit <- to.upper(cgf$gp$predict(betas, type = "response")) #need upper triangle of the predict
        
    zgf <- gglm.fit(g, response = NULL, vars = predsz$vars, factors = predsz$factors, 
                        varfuns = predsz$varfuns, facfuns = predsz$facfuns, 
                        prior = prior, start = start, 
                        offset = offset, epsilon = epsilon, 
                        maxit = maxit, trace = FALSE, family = c("binomial"), ...)
    gammas <-zgf$coef
    zeta <- logit(to.upper(zgf$gp$predict(gammas, type = "response")))  #need inverse link function and upper triangle of this 
            
    #here, we turn the upper right triangle of adjacency matrix into a vector.  This allows us 
    #to create an index denoting where edges exist/do not exist 
            
    adj <- to.upper(as_adjacency_matrix(g, type = "upper", attr = NULL, sparse = FALSE))
    edge <- which(adj == 1)
    pis <- vector(,length(adj))
	t <- complementer(g)
    itermax <- 5 #change this
	lhood <- -Inf 
	tol <- 1e-4


	for (iter in 1:itermax){                 
    
    #estimates for pi
		pis <- (1 + exp(-cfit - zeta))^-1   #if this is 1 then the offset, log(1-pis), will break.  
	#when there is an edge
		pis[edge] <- 0
    
    #update beta
		cgf <- gglm.fit(g, response = response, vars = predsc$vars, factors = predsc$factors, 
                        varfuns = predsc$varfuns, facfuns = predsc$facfuns, 
                        prior = prior, start = start, 
                        offset = log(1 - pis), epsilon = epsilon, 
                        maxit = maxit, trace = FALSE, family = c("poisson"), ...)
		betas <- cgf$coef
		cfit <- to.upper(cgf$gp$predict(betas, type = "response"))
		blhood <- cgf$lhood  #deviance in the future
		
    #update gamma
		E(t)$resp <- pis[-edge]

		zgf <- gglm.fit(t, response = "resp", vars = predsz$vars, factors = predsz$factors, 
                        varfuns = predsz$varfuns, facfuns = predsz$facfuns, 
                        prior = prior, start = start, 
                        offset = offset, epsilon = epsilon, 
                        maxit = maxit, trace = FALSE, family = c("binomial"), ...)
		gammas <-zgf$coef
		zfit <- to.upper(zgf$gp$predict(gammas, type = "response")) 
		zeta <- logit(zfit)
		zlhood <- zgf$lhood
	
	#look at likelihood	
		lhood.new <- zlhood + blhood
	
		if (lhood.new < lhood || abs((lhood.new - lhood) / lhood.new) < epsilon) break
			
		lhood <- lhood.new
		if (iter >= itermax) break
		}
   }      
    
  #changing output. Will need to consider ZIP formulation. this is a work in progress!  
  if(exists("preds")) { 
	out <- list(coefficients = gglmfit$coef, 
              residuals = 1, #(y-mu)/mu.eta(eta)
              fitted.values = gglmfit$gp$predict(gglmfit$coef, type = "response"),
              effects = 0,
              rank = 1,
              qr = 1,
              #family = family(),
              linear.predictors = 1, #eta
              deviance = 1, #family$deviance
              null.deviance = 1,
              aic = 1, #family$aic
              iter = 1, #eventually gglmfit$it
              weights = 1,
              prior.weights = 1,
              df.residual = 1,
              df.null = 1,
              y = g,
              converged = 1,
              boundary = 1, 
              call = cl,
              formula = formula, 
              terms = terms,
              data = environment(formula),
              offset = offset,
              control = 1,
              method = 'gglm.fit',
              contrasts = 1,
              xlevels = 1)
  
	class(out) <- c("glm","lm")
	}

	if(exists("predsz")) {
	
	fitted.values.vec = (1-pis)*cfit
	
	n <- length(pis)
	fitted.values.mat = matrix(0, n, n)
	mindex <- matrix(sequence(n^2), n, n, byrow = TRUE)
	fitted.values.mat[mindex[lower.tri(mindex)]] <- fitted.values.vec
	fitted.values.mat[lower.tri(fitted.values.mat)] <- t(fitted.values.mat)[lower.tri(t(fitted.values.mat))]

	
	out <- list(coefficients = list(count = betas, zero = gammas), 
			residuals = 1, 
			fitted.values.mat = fitted.values.mat,    #change to matrix
			fitted.values.vec = fitted.values.vec,
			optim = 1, 
			method = 'EM', 
			control = 1, 
			start = 1, 
			weights = 1,
			#if (identical(as.vector(weights),rep.int(1L, n))) NULL else weights, 
			offset = offset, #two offsets?
			#zero = if (identical(offsetz, rep.int(0, n))) NULL else offsetz), 
			#n = nobs, 
			#df.null = nobs - 2, 
			#df.residual = nobs - (kx + kz + (dist == "negbin")), 
			#terms = list(count = mtX, zero = mtZ, full = mt), 
			#theta = theta, 
			#SE.logtheta = SE.logtheta, 
			loglik = lhood, 
			#vcov = vc, 
			#dist = dist, link = linkstr, linkinv = linkinv, 
			converged = 1, 
			call = cl, formula = formula)
			#levels = .getXlevels(mt, mf), 
			#contrasts = list(count = attr(X, "contrasts"), 
            #zero = attr(Z, "contrasts")))
    
    #class(out) <- "zeroinfl"
    }

  out
  
}
  


