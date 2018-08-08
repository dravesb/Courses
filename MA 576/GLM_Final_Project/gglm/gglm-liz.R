library(Matrix) # sparse matrices
library(igraph)

EPSILON <- 1e-6 # for regularization
TAU2 <- 100 # default prior variance

SWEEP <- function (n, m, A, ind) # A is n-by-m double
  .Call(G_sweep, as.integer(n), as.integer(m), A, as.integer(ind - 1))

build.null.prior <- function (d, tau2=TAU2)
  list(mean = numeric(d), Sigma.inv = rep(1 / tau2, d))

lprior <- function (beta, prior) {
  if (is.vector(prior$Sigma.inv))
    return(-.5 * sum((beta - prior$mean) ^ 2 * prior$Sigma.inv))
  C.inv <- chol(prior$Sigma.inv)
  -.5 * crossprod(C.inv %*% (beta - prior$mean))
}


# [ update based on *active set method* for beta[const.ind] > 0 ]
update.beta <- function (d, const.ind, active, beta, st) {
  m <- st$m; S <- st$V
  lambda <- -m[active] # Lagrange multipliers
  if (any(lambda < 0)) {
    j <- which(lambda < 0)
    active <- active[-j] # remove constraints
  }

  # compute step:
  non.active <- which(is.na(match(1:d, active)))
  C <- chol(S)
  V <- chol2inv(C)
  if (length(active) > 0)
    V <- matrix(SWEEP(d, d, as.double(V), active), nrow = d)
  #if (any(is.na(V))) {
  #  warning("NAs in covariance: setting to zero")
  #  V[is.na(V)] <- 0 # FIXME: workaround
  #}
  p <- V[non.active, non.active] %*% m[non.active] # step
  if (any(is.na(p))) stop("invalid step")

  alpha.min <- 1
  cand <- which(non.active %in% const.ind & p < 0)
  if (length(cand) > 0) {
    offset <- EPSILON
    alpha <- (offset - beta[non.active[cand]]) / p[cand]
    alpha.min <- min(1, alpha)
    if (alpha.min < 1) { # any blocking constraints?
      j <- which(alpha == alpha.min)
      active <- c(active, non.active[cand[j]]) # activate constraint
    }
  }
  if (is.na(alpha.min))
    stop("invalid min alpha")
  beta[non.active] <- beta[non.active] + alpha.min * p

  list(beta=beta, C=C, active=active)
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

    reorder = function (fac, beta) {
      if (length(beta) != nlevel[fac]) stop("invalid beta length")
      ix <- sort(beta, index=TRUE)$ix
      f <- integer(vcount(graph))
      for (i in 1:length(beta))
        f[facdata[,fac] == ix[i] - 1] <- i - 1 # zero-based
      .Call(G_setfactor_gglmparamset, .pset, as.integer(fac - 1),
            as.integer(f), PACKAGE="gglm")
      return(ix)
    },

    remap = function (fac) {
      k <- nlevel[fac]
      .Call(G_remap_gglmparamset, .pset, as.integer(fac - 1), integer(k),
            integer(k), PACKAGE="gglm")
    },

    ppl = function (beta, offset = NULL, type = "edge") {
      t <- match(type, c("edge", "degree"))
      if (is.na(t)) stop("invalid PPL type")
      t <- as.integer(t - 1)
      .Call(G_ppl_gglmparamset, .pset, as.numeric(beta), as.numeric(offset),
            t, PACKAGE="gglm")
    },

    predict = function (beta, offset = NULL, type = "link",
                        ignore.non.edges = FALSE) {
      n <- vcount(graph)
      r <- match(type, c("link", "response", "variance", "cumulant"))
      if (is.na(r))
        stop("invalid type")
      r <- as.integer(r - 1)

      if (ignore.non.edges) {
        eta <- .Call(G_predict_gglmparamset, .pset, as.numeric(beta),
                   as.numeric(offset), r, as.integer(ignore.non.edges),
                   numeric(ecount(graph)), PACKAGE="gglm")
        el <- get.edgelist(graph)
        Matrix::sparseMatrix(dims=c(n, n), i=pmin(el[,1], el[,2]),
                             j=pmax(el[,1], el[,2]), x=eta, symmetric=TRUE)
      }
      else {
        eta <- .Call(G_predict_gglmparamset, .pset, as.numeric(beta),
                   as.numeric(offset), r, as.integer(ignore.non.edges),
                   numeric(n ^ 2), PACKAGE="gglm")
        matrix(eta, nrow=n, ncol=n)
      }
    },

    # 'stats' computes sufficient statistics for 'fit':
    # st$m = X' * (y - mu(beta)) and st$V = X' * diag(W(beta)) * X
    # if 'initialize', st$m = X' * (y - mu) + X' * diag(W(mu)) * eta(mu)
    stats = function (beta, offset = NULL, ignore.non.edges = FALSE,
                      initialize = FALSE) {
      p <- ncol(vardata) + sum(nlevel)
      m <- numeric(p)
      V <- numeric(p ^ 2)
      # write m and V in-place:
      if (ignore.non.edges)
        .Call(G_stats_edges_gglmparamset, .pset, as.numeric(beta),
              as.numeric(offset), as.integer(initialize),
              m, V, PACKAGE="gglm")
      else
        .Call(G_stats_pairs_gglmparamset, .pset, as.numeric(beta),
              as.numeric(offset), as.integer(initialize),
              m, V, PACKAGE="gglm")
      list(m=m, V=matrix(V, nrow=p))
    },

    lhood = function (beta, offset=NULL, ignore.non.edges=FALSE) {
      .Call(G_lhood_gglmparamset, .pset, as.numeric(beta), as.numeric(offset),
            as.integer(ignore.non.edges), PACKAGE="gglm")
    },

    label.stats = function (fac, beta, pp=NULL, offset=NULL,
                            ignore.non.edges=FALSE) {
      if (fac < 0 || fac > length(nlevel)) stop("invalid factor")
      if (is.null(pp)) pp <- rep(1, nlevel[fac])
      n <- vcount(graph)
      l <- .Call(G_labelstats_gglmparamset, .pset, as.integer(fac - 1),
                 as.numeric(beta), as.numeric(offset), log(pp),
                 as.integer(ignore.non.edges),
                 integer(n), numeric(n * nlevel[fac]), PACKAGE = "gglm")
      matrix(l, ncol=n)
    },

    factor.optimize = function (optfac, beta, pp=NULL, offset=NULL,
                            min.card.constraint=1) {
      if (is.null(pp)) pp <- rep(1, nlevel[optfac])
      l <- label.stats(optfac, beta, pp, offset)
      lhood <- 0
      sigma <- .Call(G_labeloptim_gglmparamset, as.integer(nlevel[optfac]),
                     as.integer(vcount(graph)), as.numeric(l),
                     as.integer(min.card.constraint), integer(nlevel[optfac]),
                     facdata[, optfac], lhood, PACKAGE = "gglm")
      setfactor(optfac, sigma) # facdata[, optfac] <<- sigma
      return(l)
    },


    fit = function (start = NULL, prior = NULL, constrained = NULL,
                    offset = NULL, ignore.non.edges = FALSE,
                    initialize = FALSE,
                    epsilon = 1e-8, maxit = 25, trace = FALSE) {
      d <- ncol(vardata) + sum(nlevel)
      const.ind <- integer(0) # constrained indices
      if (!is.null(constrained)) {
        if (length(constrained) != d) stop("invalid 'constrained' length")
        const.ind <- which(constrained)
      }
      if (!is.null(prior) && length(prior$mean) != d)
        stop("invalid 'prior' dimension")
      if (!is.null(start))
        beta <- start
      else
        beta <- numeric(d)
      active <- const.ind[which(beta[const.ind] <= 0)]

      lhood.cur <- -Inf
      it <- ifelse(initialize, 0, 1)
      repeat {
        st <- stats(beta, offset, ignore.non.edges, initialize)
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

        if (initialize) {
          # TODO: constrained
          C <- chol(st$V)
          beta.new <- backsolve(C, backsolve(C, st$m, transpose=TRUE))
        }
        else {
          if (is.null(constrained)) {
            C <- chol(st$V)
            beta.new <- beta + backsolve(C, backsolve(C, st$m, transpose=TRUE))
          }
          else {
            u <- update.beta(d, const.ind, active, beta, st)
            C <- u$C; active <- u$active
            beta.new <- u$beta
          }
        }

        lhood.new <- lhood(beta.new, offset, ignore.non.edges)
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
                      constrained = NULL,
                      optim.factors = NULL, # list with probs
                      fit.factor.probs = NULL, # T/F for each optim factor
                      start = NULL, offset = NULL,
                      epsilon = 1e-8, maxit = 25, trace = FALSE, ...) {
  gp <- gglmParamSet(g, vars, factors, varfuns, facfuns, ...)
  d <- ncol(gp$vardata) + sum(gp$nlevel)
  if (!is.null(optim.factors)) {
    if (length(optim.factors) != length(factors))
      stop("'factors' and 'optim.factors' have inconsistent lengths")
    if (!is.null(fit.factor.probs) &&
        length(fit.factor.probs) != length(factors))
      stop("'optim.factors' and 'fit.factor.probs' have inconsistent lengths")
    pz <- optim.factors # normalized weights
    for (f in seq_along(optim.factors)) {
      if (!is.null(optim.factors[[f]]))
        pz[[f]] <- pz[[f]] / sum(pz[[f]])
    }
  }
  if (!is.null(constrained) && length(constrained) != d)
    stop("inconsistent length of 'constrained'")

  fprior <- prior
  if (is.function(prior)) fprior <- prior(start)
  if (is.null(start)) { # initialize?
    res <- gp$fit(start, fprior, constrained, offset,
                  initialize = TRUE, trace = trace)
    beta <- res$coef; lhood <- res$lhood;
  }
  else {
    beta <- start; lhood <- -Inf
  }
  if (!is.null(fprior))
    lhood <- lhood + lprior(beta, fprior) # log joint

  it <- 1
  repeat {
    # [ fit factors ]
    if (!is.null(optim.factors)) {
      for (f in seq_along(optim.factors)) {
        if (!is.null(optim.factors[[f]]))
          gp$factor.optimize(f, beta, pz[[f]])
      }
    }
    # [ fit beta ]
    if (is.function(prior)) fprior <- prior(beta)
    res <- gp$fit(beta, fprior, constrained, offset, maxit = 1)
    beta.new <- res$coef; lhood.new <- res$lhood
    if (!is.null(fprior))
      lhood.new <- lhood.new + lprior(beta.new, fprior) # log joint

    # [ fit factor probs ]
    if (!is.null(optim.factors) && !is.null(fit.factor.probs)) {
      for (f in seq_along(optim.factors)) {
        if (!is.null(pz[[f]]) && fit.factor.probs[f]) {
          z <- gp$facdata[, f]
          cs <- tapply(z, z, length) + optim.factors[[f]] - 1
          pz[[f]] <- cs / sum(cs)
          lhood.new <- lhood.new + sum((optim.factors[[f]] - 1) *
                                       log(pz[[f]]))
        }
      }
    }
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
gglm <- function(formula, response = NULL, family = "binomial",  epsilon = 1e-8, maxit = 25, offset = NULL, prior = NULL, start = NULL, ...) {
  
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
              family = family(),
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


