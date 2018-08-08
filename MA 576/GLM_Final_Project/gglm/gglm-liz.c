#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#ifdef _OPENMP
#include <omp.h>
#include <R_ext/MathThreads.h>
#endif

#include <math.h>
#include <stdlib.h>
#include <igraph/igraph.h>

/*
 * _GGLMParamSet_: implements `node(x, fun=c('+', '*'))` and
 * `node.factor(x, fun=c('+','*'))`
 */

typedef struct {
  igraph_t g;
  int n; /* vcount(g) */
  int nvar; /* # non-factor predictors */
  int nfac; /* # factor predictors */
  int family; /* exponential family */
  int *nlevel; /* for factors: #levels */
  int slevel; /* sum(nlevel), for computational convenience (caching) */
  int *varfun; /* is function '*' [0] or '+' [1]? (for each variable) */
  int *facfun; /* is function '*' [0] or '+' [1]? (for each factor) */
  double *vardata; /* nobs by nvar */
  int *facdata; /* nobs by nfac; zero-based index */
  int nresponse; /* if nresponse > 0, use 'response' */
  double wresponse; /* in binomial case, 'wresponse' is size */
  igraph_vector_t response;
} GGLMParamSet;

#define NOEDGE (-1)
#define GGLMPARAMSET(sp) ((GGLMParamSet *) R_ExternalPtrAddr(sp))  /*L: the external pointer*/


/* [===[ Exponential Families ]===] */

static double identity (double x) { return x; }

/* binomial */
static double linkfun_binomial (double mu) { return log(mu / (1. - mu)); }
static double linkinv_binomial (double eta) { return 1. / (1 + exp(-eta)); }
static double variance_binomial (double mu) { return mu * (1. - mu); }
# define LOGEPS2 (-37) /* floor(log(eps / 2)) */
static double cumulant_binomial (double eta) { /* log(1 + exp(eta)) */
  double l = 0, x = eta;
  if (eta > 0) { l = eta; x = -eta; }
  if (x < LOGEPS2) return l;
  return l + log(1 + exp(x));
}
static double initialize_binomial (double y) { return (y + .5) / 2.; }

/* poisson */
#define linkfun_poisson log
#define linkinv_poisson exp
#define variance_poisson identity
#define cumulant_poisson exp
static double initialize_poisson (double y) { return y + .1; }

/* gaussian */
#define linkfun_gaussian identity
#define linkinv_gaussian identity
static double variance_gaussian (double mu) { (void) mu; return 1.; }
static double cumulant_gaussian (double eta) { return .5 * eta * eta; }
#define initialize_gaussian identity

/* TODO: gamma! */

typedef double (*FamilyFunc) (double);
static FamilyFunc linkfun[3] = {linkfun_binomial, linkfun_poisson,
  linkfun_gaussian};
static FamilyFunc linkinv[3] = {linkinv_binomial, linkinv_poisson,
  linkinv_gaussian};
static FamilyFunc variance[3] = {variance_binomial, variance_poisson,
  variance_gaussian};
static FamilyFunc cumulant[3] = {cumulant_binomial, cumulant_poisson,
  cumulant_gaussian};
static FamilyFunc initialize[3] = {initialize_binomial, initialize_poisson,
  initialize_gaussian};


static double compute_etaij_factor (GGLMParamSet *P, int i, int j,
    double *beta, double *offset, int loffset, int fac) {
  double xi, xj, sijk, etaij = 0;
  int fi, fj, c, k;

  for (k = 0; k < P->nvar; k++) {
    xi = P->vardata[k * P->n + i];
    xj = P->vardata[k * P->n + j];
    sijk = (P->varfun[k]) ? xi * xj : xi + xj; /*L: varfun 0 or 1, change this when we include other options!*/
    etaij += sijk * beta[k];
  }
  c = P->nvar; /* factor column baseline in beta */
  for (k = 0; k < P->nfac; k++) {
    if (k != fac) {
      fi = P->facdata[k * P->n + i];
      fj = P->facdata[k * P->n + j];
      if (P->facfun[k]) /* '*' ? */
        etaij += (fi == fj) * beta[c + fi]; /*facfun 0 or 1, change this if we include other options.*/
      else /* '+' */
        etaij += beta[c + fi] + beta[c + fj];
    }
    c += P->nlevel[k];
  }
  if (loffset > 0) {
    if (loffset == 1) etaij += *offset;
    else {
      int s = i * P->n - i * (i + 1) / 2; /* index pair shift */
      etaij += offset[s + j - i - 1];
    }
  }
  return etaij;
}

#define compute_etaij(P,i,j,beta,offset,loffset) \
  compute_etaij_factor(P,i,j,beta,offset,loffset,-1)



/* [===[ igraph ]===] (from rinterface.c) */

static int R_SEXP_to_vector_copy (SEXP sv, igraph_vector_t *v) {
  return igraph_vector_init_copy(v, REAL(sv), length(sv));
}

static int R_SEXP_to_vector (SEXP sv, igraph_vector_t *v) {
  v->stor_begin=REAL(sv);
  v->stor_end=v->stor_begin+length(sv);
  v->end=v->stor_end;
  return 0;
}

static int R_SEXP_to_igraph (SEXP graph, igraph_t *res) {
  res->n=REAL(VECTOR_ELT(graph, 0))[0];
  res->directed=LOGICAL(VECTOR_ELT(graph, 1))[0];
  R_SEXP_to_vector(VECTOR_ELT(graph, 2), &res->from);
  R_SEXP_to_vector(VECTOR_ELT(graph, 3), &res->to);
  R_SEXP_to_vector(VECTOR_ELT(graph, 4), &res->oi);
  R_SEXP_to_vector(VECTOR_ELT(graph, 5), &res->ii);
  R_SEXP_to_vector(VECTOR_ELT(graph, 6), &res->os);
  R_SEXP_to_vector(VECTOR_ELT(graph, 7), &res->is);
  
  /* attributes */
  REAL(VECTOR_ELT(VECTOR_ELT(graph, 8), 0))[0] = 1; /* R objects refcount */
  REAL(VECTOR_ELT(VECTOR_ELT(graph, 8), 0))[1] = 0; /* igraph_t objects */
  res->attr=VECTOR_ELT(graph, 8);
  
  return 0;
}

static SEXP R_igraph_getListElement(SEXP list, const char *str) {
  SEXP names = getAttrib(list, R_NamesSymbol);
  int i;
  for (i = 0; i < length(list); i++)
    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) /*L: if attribute name matches, accesses element i of the vector*/
      return VECTOR_ELT(list, i);  /*L: returns element i in vector list*/
  return R_NilValue;
}

/* numeric edge attributes for the whole graph */    /*is this the WHOLE graph?*/
static int R_igraph_attribute_get_edge_attr (const igraph_t *graph,
    const char *name, igraph_vector_t *value) {
  SEXP eal = VECTOR_ELT(graph->attr, 3); /*L: creates vector, eal, which is the 3rd element in graph$attr*/
  SEXP ea = R_igraph_getListElement(eal, name); /*L: ea = element i in vector list*/
  igraph_vector_t newvalue;

  if (ea == R_NilValue)
    IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);
  PROTECT(ea = coerceVector(ea, REALSXP));
  R_SEXP_to_vector_copy(coerceVector(ea, REALSXP), &newvalue); /*L: store in newvalue*/
  igraph_vector_destroy(value); /*L: destroy old value*/
  *value = newvalue; 
  UNPROTECT(1);
  return 0;
}


/* [===[ Auxiliary ]===] */

/* remap 'sigma' (length 'n'); 'ix' and buffer 'revind' have length 'k' */
static void remap (int n, int k, int *sigma, int *ix, int *revind) {
  int i, l, c = 0;

  /* [ init ] */
  for (l = 0; l < k; l++) {
    ix[l] = 0;
    revind[l] = -1;
  }

  /* [ find order ] */
  for (i = 0; i < n; i++) {
    l = sigma[i];
    if (ix[l] == 0) /* first visit? */
      revind[l] = c++; /* reverse order of appearance of `l` */
    ix[l] = 1; /* mark as visited */
  }
  for (l = 0; l < k; l++) {
    if (revind[l] == -1) /* no appearance? */
      revind[l] = c++;
  }

  /* [ permute ] */
  for (i = 0; i < n; i++) sigma[i] = revind[sigma[i]];
  for (l = 0; l < k; l++) ix[revind[l]] = l + 1; /* one-based */
}

/* SWEEP operator */
static int sweepaux (int n, int m, double *A, int k) {
  int i;
  double D = A[k * (n + 1)];
  if (D == 0) return 1; /* error */
  D = 1.0 / D;
  F77_NAME(dscal)(&m, &D, A + k, &n); /* A[k,] = A[k,] * D */
  for (i = 0; i < n; i++) {
    if (i != k) {
      double B = -A[i + k * n]; /* -A[i,k] */
      /* A[i,] = A[i,] + B * A[k,]: */
      F77_NAME(daxpy)(&m, &B, A + k, &n, A + i, &n);
      A[i + k * n] = B * D;
    }
  }
  A[k * (n + 1)] = D; /* A[k,k] = D */
  return 0;
}

SEXP sweep (SEXP sn, SEXP sm, SEXP sa, SEXP sind) {
  int n = INTEGER(sn)[0];
  int m = INTEGER(sm)[0];
  double *A = REAL(sa);
  int *ind = INTEGER(sind);
  int i, k = length(sind);
  for (i = 0; i < k; i++) sweepaux(n, m, A, ind[i]);
  return sa;
}



/* [===[ GGLMParamSet ]===] */

static void finalize_gglmparamset (SEXP sp) {
  GGLMParamSet *P = GGLMPARAMSET(sp);
  if (P->nresponse > 0)
    igraph_vector_destroy(&P->response);   /* L: we destroy the response we created? */
  free(P);
}

/* stored in object slots: nlevel, varfun, facfun, vardata, facdata */
SEXP new_gglmparamset (SEXP sg, SEXP snvar, SEXP snfac, SEXP sfamily,
    SEXP snlevel, SEXP sslevel, SEXP svarfun, SEXP sfacfun,
    SEXP svardata, SEXP sfacdata, SEXP sresponse, SEXP swresponse) {
  SEXP sp;
  int nresp = length(sresponse);
  GGLMParamSet *P = (GGLMParamSet *) malloc(sizeof(GGLMParamSet)  /*L: new struct that we have defined above, GGLMParamSet*/
      + nresp * sizeof(igraph_vector_t)); /*L: malloc allows us to make it an unknown size*/
  R_SEXP_to_igraph(sg, &P->g); /*L: similar to P$n notation in R*/
  P->n = igraph_vcount(&P->g); /*L: Remember, P is a pointer to the new struct. To get pieces of it we need &P*/
  P->nvar = asInteger(snvar);
  P->nfac = asInteger(snfac);
  P->family = asInteger(sfamily);
  P->nlevel = INTEGER(snlevel);
  P->slevel = asInteger(sslevel);
  P->varfun = INTEGER(svarfun);
  P->facfun = INTEGER(sfacfun);
  P->vardata = REAL(svardata);
  P->facdata = INTEGER(sfacdata);
  P->nresponse = nresp;
  P->wresponse = asReal(swresponse);
  if (nresp > 0) { /* L: response was provided? */
    int e, ne = igraph_ecount(&P->g);
    igraph_vector_init(&P->response, ne); /*L: initialize response vector*/
    R_igraph_attribute_get_edge_attr(&P->g, CHAR(STRING_ELT(sresponse, 0)), /*L: updates the response attribute in P->response*/
        &P->response); /*L: uses function define above*/
    /* adjust response for binomial case: y -> y / w */
    if (P->family == 0 && P->wresponse > 0) { /* binomial with size? */
      for (e = 0; e < ne; e++)
        VECTOR(P->response)[e] /= P->wresponse;
    }
  }

  PROTECT(sp = R_MakeExternalPtr(P, install("gglm"), R_NilValue)); /*L: huh?*/
  R_RegisterCFinalizerEx(sp, finalize_gglmparamset, TRUE); /*runs finalize_gglmparamset*/
  UNPROTECT(1);
  return(sp); /*pointer to all of P-> */
}


/* compute X * beta;
 * type = 0 (link), 1 (response), 2 (variance), 3 (cumulant) */
SEXP predict_gglmparamset (SEXP sp, SEXP sbeta, SEXP soffset, SEXP stype,
    SEXP signore, SEXP seta) {
  GGLMParamSet *P = GGLMPARAMSET(sp);
  double *beta = REAL(sbeta);
  double *offset = REAL(soffset);
  int type = asInteger(stype);
  int ignore_non_edges = asInteger(signore);
  double *eta = REAL(seta);
  int loffset = length(soffset);
  int i, j;
  double etaij;

  if (ignore_non_edges) {
    int ne = igraph_ecount(&P->g);
    for (e = 0; e < ne; e++) { /* for all edges */
      igraph_edge(&P->g, e, &i, &j);
      etaij = compute_etaij(P, i, j, beta, offset, loffset); /* link */
      if (type > 0) {
        if (type == 3) /* cumulant? */
          etaij = cumulant[P->family](etaij);
        else {
          etaij = linkinv[P->family](etaij); /* response */
          if (type == 2) /* variance? */
            etaij = variance[P->family](etaij);
        }
      }
      eta[e] = etaij;
    }
  }
  else {
    for (i = 0; i < P->n; i++) { /* for all pairs of vertices (i, j) */
      eta[i * (P->n + 1)] = 0; /* diagonal entries */
      for (j = i + 1; j < P->n; j++) {
        etaij = compute_etaij(P, i, j, beta, offset, loffset);
        if (type > 0) {
          if (type == 3) /* cumulant? */
            etaij = cumulant[P->family](etaij);
          else {
            etaij = linkinv[P->family](etaij); /* response */
            if (type == 2) /* variance? */
              etaij = variance[P->family](etaij);
          }
        }
        eta[i + P->n * j] = eta[j + P->n * i] = etaij;
      }
    }
  }
  return seta;
}

/* Posterior Predictive Loss: type = 0 (edges) or 1 (degrees) */
SEXP ppl_gglmparamset (SEXP sp, SEXP sbeta, SEXP soffset, SEXP stype) {
  GGLMParamSet *P = GGLMPARAMSET(sp);
  double *beta = REAL(sbeta);
  double *offset = REAL(soffset);
  int type = asInteger(stype);
  int loffset = length(soffset);
  int i, j, e;
  double yij, etaij, muij, muetaij;
  double ppl = 0.;
  double wij = (P->family == 0 && P->wresponse > 0) ? P->wresponse : 1.;
  igraph_vector_t adj;
  int *isneighbor = (int *) malloc(P->n * sizeof(int));
  double *dmu, *dvar;
  double pplnorm = P->n;
  if (type == 0) pplnorm *= (P->n - 1) / 2.; /* #pairs */

  if (type == 1) {
    dmu = (double *) malloc(P->n * sizeof(double));
    dvar = (double *) malloc(P->n * sizeof(double));
    for (i = 0; i < P->n; i++) dmu[i] = dvar[i] = 0.;
  }
  igraph_vector_init(&adj, P->n);
  for (i = 0; i < P->n; i++) {
    igraph_neighbors(&P->g, &adj, i, IGRAPH_ALL);
    for (j = 0; j < P->n; j++) isneighbor[j] = 0;
    for (j = 0; j < igraph_vector_size(&adj); j++)
      isneighbor[(int) VECTOR(adj)[j]] = 1;

    for (j = i + 1; j < P->n; j++) {
      yij = 0.;
      if (isneighbor[j]) {
        yij = 1.;
        if (P->nresponse > 0) {
          igraph_get_eid(&P->g, &e, i, j, IGRAPH_UNDIRECTED, 1);
          yij = VECTOR(P->response)[e];
        }
      }

      etaij = compute_etaij(P, i, j, beta, offset, loffset);
      muij = linkinv[P->family](etaij);
      muetaij = wij * variance[P->family](muij);
      yij -= wij * muij;
      if (type == 0)
        ppl += (yij * yij + muetaij) / pplnorm;
      else { /* type == 'degree' */
        dmu[i] += yij; dmu[j] += yij;
        dvar[i] += muetaij; dvar[j] += muetaij;
      }
    }
  }
  free(isneighbor);
  igraph_vector_destroy(&adj);
  if (type == 1) { /* type == 'degree'? */
    ppl = 0.;
    for (i = 0; i < P->n; i++)
      ppl += (dmu[i] * dmu[i] + dvar[i]) / pplnorm;
    free(dmu); free(dvar);
  }
  return ScalarReal(ppl);
}


/* compute stats: m := X' * (y - mu(beta)) and V := X' * W(beta) * X */
/* `ignore_non_edges` is TRUE */
SEXP stats_edges_gglmparamset (SEXP sp, SEXP sbeta, SEXP soffset, SEXP sinit,
    SEXP sm, SEXP sV) {
  GGLMParamSet *P = GGLMPARAMSET(sp);
  double *beta = REAL(sbeta);
  double *offset = REAL(soffset);
  int loffset = length(soffset);
  int initial_stats = asInteger(sinit);
  double wij;
  int e, k, j, l, c;
  double *m = REAL(sm);
  double *V = REAL(sV);
  double *mu, *mueta;
  int *edgei, *edgej;
  int ne = igraph_ecount(&P->g);
  int nt = P->nvar + P->slevel;
#ifdef _OPENMP
  int nthreads = omp_get_num_procs();
#endif

  /* TODO: more general weights */
  wij = (P->family == 0 && P->wresponse > 0) ? P->wresponse : 1.;

  /* [ init ] */
  edgei = (int *) malloc(ne * sizeof(int));
  edgej = (int *) malloc(ne * sizeof(int));
  mu = (double *) malloc(ne * sizeof(double));
  mueta = (double *) malloc(ne * sizeof(double));
  for (e = 0; e < ne; e++)
    igraph_edge(&P->g, e, &edgei[e], &edgej[e]);

  /* FIXME: move to R, at gglm.fit */
  /* [ init ] */
  for (k = 0; k < P->nvar; k++) {
    m[k] = 0;
    for (j = k; j < nt; j++) V[k + nt * j] = 0;
  }
  c = P->nvar; /* factor column baseline */
  for (k = 0; k < P->nfac; k++) {
    for (l = 0; l < P->nlevel[k]; l++) {
      m[c + l] = 0;
      for (j = c + l; j < nt; j++) V[(c + l) + nt * j] = 0;
    }
    c += P->nlevel[k];
  }


#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads) default(none) \
  shared(P, wij, initial_stats, ne, nt, edgei, edgej, mu, mueta, m, V, \
      initialize, variance, linkfun, linkinv, beta, loffset, offset)
#endif
  {
    int k, l, c, d, e, f, r;
    int fi, fj;
    double yij, xi, xj, sijk, sijl, etaij, muij, muetaij;
    double mt, Vt;

    /* [ init: m <- m + X' * y ] */
    /* FIXME: move to R, at gglm.fit */
    /* var block */
    for (k = 0; k < P->nvar; k++) {
      mt = 0;
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (e = 0; e < ne; e++) { /* for all edges */
        yij = (P->nresponse > 0) ? VECTOR(P->response)[e] : 1.;
        yij *= wij; /* restore counts in binomial case */
        xi = P->vardata[k * P->n + edgei[e]];
        xj = P->vardata[k * P->n + edgej[e]];
        sijk = (P->varfun[k]) ? xi * xj : xi + xj;
        mt += sijk * yij;
      }
#ifdef _OPENMP
#pragma omp atomic
#endif
      m[k] += mt;
    }
    /* fac block */
    c = P->nvar;
    for (k = 0; k < P->nfac; k++) {
      for (l = 0; l < P->nlevel[k]; l++) {
        mt = 0;
#ifdef _OPENMP
#pragma omp for nowait
#endif
        for (e = 0; e < ne; e++) { /* for all edges */
          yij = (P->nresponse > 0) ? VECTOR(P->response)[e] : 1.;
          yij *= wij; /* restore counts in binomial case */
          fi = P->facdata[k * P->n + edgei[e]];
          fj = P->facdata[k * P->n + edgej[e]];
          sijk = (P->facfun[k]) ? (fi == l) * (fj == l) :
                                  (fi == l) + (fj == l);
          mt += sijk * yij;
        }
#ifdef _OPENMP
#pragma omp atomic
#endif
        m[c + l] += mt;
      }
      c += P->nlevel[k];
    }


    /* [ compute mu and mueta ] */
    if (initial_stats) {
#ifdef _OPENMP
#pragma omp for nowait 
#endif
      for (e = 0; e < ne; e++) { /* for all edges */
        yij = (P->nresponse > 0) ? VECTOR(P->response)[e] : 1.;
        muij = initialize[P->family](yij);
        muetaij = variance[P->family](muij);
        muij -= linkfun[P->family](muij) * muetaij;
        mu[e] = muij * wij;
        mueta[e] = muetaij * wij;
      }
    }
    else {
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (e = 0; e < ne; e++) { /* for all edges */
        etaij = compute_etaij(P, edgei[e], edgej[e], beta, offset, loffset);
        muij = linkinv[P->family](etaij);
        mu[e] = muij * wij;
        mueta[e] = variance[P->family](muij) * wij;
      }
    }


    /* [ update mean ] */
    /* var block */
    for (k = 0; k < P->nvar; k++) {
      mt = 0;
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (e = 0; e < ne; e++) {
        xi = P->vardata[k * P->n + edgei[e]];
        xj = P->vardata[k * P->n + edgej[e]];
        sijk = (P->varfun[k]) ? xi * xj : xi + xj;
        mt += sijk * mu[e];
      }
#ifdef _OPENMP
#pragma omp atomic
#endif
      m[k] -= mt;
    }

    /* fac block */
    c = P->nvar;
    for (k = 0; k < P->nfac; k++) {
      for (l = 0; l < P->nlevel[k]; l++) {
        mt = 0;
#ifdef _OPENMP
#pragma omp for nowait
#endif
        for (e = 0; e < ne; e++) {
          fi = P->facdata[k * P->n + edgei[e]];
          fj = P->facdata[k * P->n + edgej[e]];
          sijk = (P->facfun[k]) ? (fi == l) * (fj == l) :
                                  (fi == l) + (fj == l);
          mt += sijk * mu[e];
        }
#ifdef _OPENMP
#pragma omp atomic
#endif
        m[c + l] -= mt;
      }
      c += P->nlevel[k];
    }


    /* [ update variance ] */
    /* var x var block */
    for (k = 0; k < P->nvar; k++) {
      for (l = k; l < P->nvar; l++) {
        Vt = 0;
#ifdef _OPENMP
#pragma omp for nowait
#endif
        for (e = 0; e < ne; e++) {
          xi = P->vardata[k * P->n + edgei[e]];
          xj = P->vardata[k * P->n + edgej[e]];
          sijk = (P->varfun[k]) ? xi * xj : xi + xj;
          xi = P->vardata[l * P->n + edgei[e]];
          xj = P->vardata[l * P->n + edgej[e]];
          sijl = (P->varfun[l]) ? xi * xj : xi + xj;
          Vt += mueta[e] * sijk * sijl;
        }
#ifdef _OPENMP
#pragma omp atomic
#endif
        V[k + nt * l] += Vt;
      }
    }

    /* var x fac block */
    for (k = 0; k < P->nvar; k++) {
      c = P->nvar;
      for (f = 0; f < P->nfac; f++) {
        for (l = 0; l < P->nlevel[f]; l++) {
          Vt = 0;
#ifdef _OPENMP
#pragma omp for nowait
#endif
          for (e = 0; e < ne; e++) {
            xi = P->vardata[k * P->n + edgei[e]];
            xj = P->vardata[k * P->n + edgej[e]];
            sijk = (P->varfun[k]) ? xi * xj : xi + xj;
            fi = P->facdata[f * P->n + edgei[e]];
            fj = P->facdata[f * P->n + edgej[e]];
            sijl = (P->facfun[f]) ? (fi == l) * (fj == l) :
                                    (fi == l) + (fj == l);
            Vt += mueta[e] * sijk * sijl;
          }
#ifdef _OPENMP
#pragma omp atomic
#endif
          V[k + nt * (c + l)] += Vt;
        }
        c += P->nlevel[f];
      }
    }

    /* fac x fac block */
    c = P->nvar;
    for (k = 0; k < P->nfac; k++) {
      for (l = 0; l < P->nlevel[k]; l++) {
        d = c;
        for (f = k; f < P->nfac; f++) {
          for (r = 0; r < P->nlevel[f]; r++) {
            Vt = 0;
#ifdef _OPENMP
#pragma omp for nowait
#endif
            for (e = 0; e < ne; e++) {
              fi = P->facdata[k * P->n + edgei[e]];
              fj = P->facdata[k * P->n + edgej[e]];
              sijk = (P->facfun[k]) ? (fi == l) * (fj == l) :
                                      (fi == l) + (fj == l);
              fi = P->facdata[f * P->n + edgei[e]];
              fj = P->facdata[f * P->n + edgej[e]];
              sijl = (P->facfun[f]) ? (fi == r) * (fj == r) :
                                      (fi == r) + (fj == r);
              Vt += mueta[e] * sijk * sijl;
            }
#ifdef _OPENMP
#pragma omp atomic
#endif
            V[(c + l) + nt * (d + r)] += Vt;
          }
          d += P->nlevel[f];
        }
      }
      c += P->nlevel[k];
    }
  }

  free(edgei); free(edgej);
  free(mu); free(mueta);
  return sm;
}


/* compute stats: m := X' * (y - mu(beta)) and V := X' * W(beta) * X */
/* `ignore_non_edges` is FALSE */
#ifdef _OPENMP
SEXP stats_pairs_gglmparamset (SEXP sp, SEXP sbeta, SEXP soffset, SEXP sinit,
    SEXP sm, SEXP sV) {
  GGLMParamSet *P = GGLMPARAMSET(sp);
  double *beta = REAL(sbeta);
  double *offset = REAL(soffset);
  int loffset = length(soffset);
  int initial_stats = asInteger(sinit);
  double wij;
  int e, k, i, j, l, c, s;
  double *m = REAL(sm);
  double *V = REAL(sV);
  double *mu, *mueta;
  int *edgei, *edgej;
  int *edgeij = NULL; /* if edgeij[p] > 0 for pair p, edgeij holds eid */
  int ne = igraph_ecount(&P->g);
  int np = P->n * (P->n - 1) / 2;
  int nt = P->nvar + P->slevel;
#ifdef _OPENMP
  int nthreads = omp_get_num_procs();
#endif

  /* TODO: more general weights */
  wij = (P->family == 0 && P->wresponse > 0) ? P->wresponse : 1.;

  /* [ init ] */
  edgei = (int *) malloc(ne * sizeof(int));
  edgej = (int *) malloc(ne * sizeof(int));
  mu = (double *) malloc(np * sizeof(double));
  mueta = (double *) malloc(np * sizeof(double));
  for (e = 0; e < ne; e++)
    igraph_edge(&P->g, e, &edgei[e], &edgej[e]);
  
  if (initial_stats) {
    edgeij = (int *) malloc(np * sizeof(int));
    for (i = 0; i < np; i++) edgeij[i] = NOEDGE;
    for (e = 0; e < ne; e++) {
      i = edgei[e]; j = edgej[e];
      if (i > j) { s = i; i = j; j = s; } /* (i, j) = sort(i, j) */
      s = i * P->n - i * (i + 1) / 2; /* index pair shift */
      edgeij[s + j - i - 1] = e;
    }
  }

  /* TODO: move to R, at gglm.fit */
  /* [ init ] */
  for (k = 0; k < P->nvar; k++) {
    m[k] = 0;
    for (j = k; j < nt; j++) V[k + nt * j] = 0;
  }
  c = P->nvar; /* factor column baseline */
  for (k = 0; k < P->nfac; k++) {
    for (l = 0; l < P->nlevel[k]; l++) {
      m[c + l] = 0;
      for (j = c + l; j < nt; j++) V[(c + l) + nt * j] = 0;
    }
    c += P->nlevel[k];
  }


#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads) default(none) \
  shared(P, wij, initial_stats, ne, nt, edgei, edgej, edgeij, mu, mueta, \
      m, V, initialize, variance, linkfun, linkinv, beta, loffset, offset)
#endif
  {
    int i, j, k, l, c, d, e, f, r, s;
    int fi, fj, gi, gj;
    double yij, xi, xj, zi, zj, sijk, sijl, etaij, muij, muetaij;
    double mt, Vt;

    /* [ init: m <- m + X' * y ] */
    /* FIXME: move to R, at gglm.fit; remove dependency on edgei, edgej */
    /* var block */
    for (k = 0; k < P->nvar; k++) {
      mt = 0;
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (e = 0; e < ne; e++) { /* for all edges */
        yij = (P->nresponse > 0) ? VECTOR(P->response)[e] : 1.;
        yij *= wij; /* restore counts in binomial case */
        xi = P->vardata[k * P->n + edgei[e]];
        xj = P->vardata[k * P->n + edgej[e]];
        sijk = (P->varfun[k]) ? xi * xj : xi + xj;
        mt += sijk * yij;
      }
#ifdef _OPENMP
#pragma omp atomic
#endif
      m[k] += mt;
    }
    /* fac block */
    c = P->nvar;
    for (k = 0; k < P->nfac; k++) {
      for (l = 0; l < P->nlevel[k]; l++) {
        mt = 0;
#ifdef _OPENMP
#pragma omp for nowait
#endif
        for (e = 0; e < ne; e++) { /* for all edges */
          yij = (P->nresponse > 0) ? VECTOR(P->response)[e] : 1.;
          yij *= wij; /* restore counts in binomial case */
          fi = P->facdata[k * P->n + edgei[e]];
          fj = P->facdata[k * P->n + edgej[e]];
          sijk = (P->facfun[k]) ? (fi == l) * (fj == l) :
                                  (fi == l) + (fj == l);
          mt += sijk * yij;
        }
#ifdef _OPENMP
#pragma omp atomic
#endif
        m[c + l] += mt;
      }
      c += P->nlevel[k];
    }


    /* [ compute mu and mueta ] */
    if (initial_stats) {
#ifdef _OPENMP
#pragma omp for nowait schedule(dynamic)
#endif
      for (i = 0; i < P->n - 1; i++) {
        s = i * P->n - i * (i + 1) / 2; /* index pair shift */
        for (j = i + 1; j < P->n; j++) {
          e = edgeij[s + j - i - 1];
          if (e != NOEDGE)
            yij = (P->nresponse > 0) ? VECTOR(P->response)[e] : 1.;
          else
            yij = 0.;
          muij = initialize[P->family](yij);
          muetaij = variance[P->family](muij);
          muij -= linkfun[P->family](muij) * muetaij;
          mu[s + j - i - 1] = muij * wij; /* note: relative to i */
          mueta[s + j - i - 1] = muetaij * wij;
        }
      }
    }
    else {
#ifdef _OPENMP
#pragma omp for nowait schedule(dynamic)
#endif
      for (i = 0; i < P->n - 1; i++) {
        s = i * P->n - i * (i + 1) / 2; /* index pair shift */
        for (j = i + 1; j < P->n; j++) {
          etaij = compute_etaij(P, i, j, beta, offset, loffset);
          muij = linkinv[P->family](etaij);
          mu[s + j - i - 1] = muij * wij; /* note: relative to i */
          mueta[s + j - i - 1] = variance[P->family](muij) * wij;
        }
      }
    }


    /* [ update mean ] */
    /* var block */
    for (k = 0; k < P->nvar; k++) {
      mt = 0;
#ifdef _OPENMP
#pragma omp for nowait schedule(dynamic)
#endif
      for (i = 0; i < P->n - 1; i++) {
        s = i * P->n - i * (i + 1) / 2; /* index pair shift */
        xi = P->vardata[k * P->n + i];
        for (j = i + 1; j < P->n; j++) {
          xj = P->vardata[k * P->n + j];
          sijk = (P->varfun[k]) ? xi * xj : xi + xj;
          mt += sijk * mu[s + j - i - 1]; /* note: relative to i */
        }
      }
#ifdef _OPENMP
#pragma omp atomic
#endif
      m[k] -= mt;
    }

    /* fac block */
    c = P->nvar;
    for (k = 0; k < P->nfac; k++) {
      for (l = 0; l < P->nlevel[k]; l++) {
        mt = 0;
#ifdef _OPENMP
#pragma omp for nowait schedule(dynamic)
#endif
        for (i = 0; i < P->n - 1; i++) {
          s = i * P->n - i * (i + 1) / 2; /* index pair shift */
          fi = P->facdata[k * P->n + i];
          for (j = i + 1; j < P->n; j++) {
            fj = P->facdata[k * P->n + j];
            sijk = (P->facfun[k]) ? (fi == l) * (fj == l) :
              (fi == l) + (fj == l);
            mt += sijk * mu[s + j - i - 1]; /* note: relative to i */
          }
        }
#ifdef _OPENMP
#pragma omp atomic
#endif
        m[c + l] -= mt;
      }
      c += P->nlevel[k];
    }


    /* [ update variance ] */
    /* var x var block */
    for (k = 0; k < P->nvar; k++) {
      for (l = k; l < P->nvar; l++) {
        Vt = 0;
#ifdef _OPENMP
#pragma omp for nowait schedule(dynamic)
#endif
        for (i = 0; i < P->n - 1; i++) {
          s = i * P->n - i * (i + 1) / 2; /* index pair shift */
          xi = P->vardata[k * P->n + i];
          zi = P->vardata[l * P->n + i];
          for (j = i + 1; j < P->n; j++) {
            xj = P->vardata[k * P->n + j];
            sijk = (P->varfun[k]) ? xi * xj : xi + xj;
            zj = P->vardata[l * P->n + j];
            sijl = (P->varfun[l]) ? zi * zj : zi + zj;
            Vt += mueta[s + j - i - 1] * sijk * sijl;
          }
        }
#ifdef _OPENMP
#pragma omp atomic
#endif
        V[k + nt * l] += Vt;
      }
    }

    /* var x fac block */
    for (k = 0; k < P->nvar; k++) {
      c = P->nvar;
      for (f = 0; f < P->nfac; f++) {
        for (l = 0; l < P->nlevel[f]; l++) {
          Vt = 0;
#ifdef _OPENMP
#pragma omp for nowait schedule(dynamic)
#endif
          for (i = 0; i < P->n - 1; i++) {
            s = i * P->n - i * (i + 1) / 2; /* index pair shift */
            xi = P->vardata[k * P->n + i];
            fi = P->facdata[f * P->n + i];
            for (j = i + 1; j < P->n; j++) {
              xj = P->vardata[k * P->n + j];
              sijk = (P->varfun[k]) ? xi * xj : xi + xj;
              fj = P->facdata[f * P->n + j];
              sijl = (P->facfun[f]) ? (fi == l) * (fj == l) :
                                      (fi == l) + (fj == l);
              Vt += mueta[s + j - i - 1] * sijk * sijl;
            }
          }
#ifdef _OPENMP
#pragma omp atomic
#endif
          V[k + nt * (c + l)] += Vt;
        }
        c += P->nlevel[f];
      }
    }

    /* fac x fac block */
    c = P->nvar;
    for (k = 0; k < P->nfac; k++) {
      for (l = 0; l < P->nlevel[k]; l++) {
        d = c;
        for (f = k; f < P->nfac; f++) {
          for (r = 0; r < P->nlevel[f]; r++) {
            Vt = 0;
#ifdef _OPENMP
#pragma omp for nowait schedule(dynamic)
#endif
            for (i = 0; i < P->n - 1; i++) {
              s = i * P->n - i * (i + 1) / 2; /* index pair shift */
              fi = P->facdata[k * P->n + i];
              gi = P->facdata[f * P->n + i];
              for (j = i + 1; j < P->n; j++) {
                fj = P->facdata[k * P->n + j];
                sijk = (P->facfun[k]) ? (fi == l) * (fj == l) :
                                        (fi == l) + (fj == l);
                gj = P->facdata[f * P->n + j];
                sijl = (P->facfun[f]) ? (gi == r) * (gj == r) :
                                        (gi == r) + (gj == r);
                Vt += mueta[s + j - i - 1] * sijk * sijl;
              }
            }
#ifdef _OPENMP
#pragma omp atomic
#endif
            V[(c + l) + nt * (d + r)] += Vt;
          }
          d += P->nlevel[f];
        }
      }
      c += P->nlevel[k];
    }
  }

  free(edgei); free(edgej);
  if (initial_stats) free(edgeij);
  free(mu); free(mueta);
  return sm;
}
#else
SEXP stats_pairs_gglmparamset (SEXP sp, SEXP sbeta, SEXP soffset, SEXP sinit,
    SEXP sm, SEXP sV) {
  GGLMParamSet *P = GGLMPARAMSET(sp);
  double *beta = REAL(sbeta);
  double *offset = REAL(soffset);
  int loffset = length(soffset);
  int initial_stats = asInteger(sinit);
  double wij;
  int e, k, i, j, l, c, s;
  double *m = REAL(sm);
  double *V = REAL(sV);
  int *edgei, *edgej;
  int *edgeij = NULL; /* if edgeij[p] > 0 for pair p, edgeij holds eid */
  int ne = igraph_ecount(&P->g);
  int np = P->n * (P->n - 1) / 2;
  int nt = P->nvar + P->slevel;

  /* TODO: more general weights */
  wij = (P->family == 0 && P->wresponse > 0) ? P->wresponse : 1.;

  /* [ init ] */
  edgei = (int *) malloc(ne * sizeof(int));
  edgej = (int *) malloc(ne * sizeof(int));
  for (e = 0; e < ne; e++)
    igraph_edge(&P->g, e, &edgei[e], &edgej[e]);
  
  if (initial_stats) {
    edgeij = (int *) malloc(np * sizeof(int));
    for (i = 0; i < np; i++) edgeij[i] = NOEDGE;
    for (e = 0; e < ne; e++) {
      i = edgei[e]; j = edgej[e];
      if (i > j) { s = i; i = j; j = s; } /* (i, j) = sort(i, j) */
      s = i * P->n - i * (i + 1) / 2; /* index pair shift */
      edgeij[s + j - i - 1] = e;
    }
  }

  /* TODO: move to R, at gglm.fit */
  /* [ init ] */
  for (k = 0; k < P->nvar; k++) {
    m[k] = 0;
    for (j = k; j < nt; j++) V[k + nt * j] = 0;
  }
  c = P->nvar; /* factor column baseline */
  for (k = 0; k < P->nfac; k++) {
    for (l = 0; l < P->nlevel[k]; l++) {
      m[c + l] = 0;
      for (j = c + l; j < nt; j++) V[(c + l) + nt * j] = 0;
    }
    c += P->nlevel[k];
  }


  {
    int i, j, k, l, c, d, e, f, r, s;
    int fi, fj, gi, gj;
    double yij, xi, xj, sijk, sijl, etaij, muij, muetaij;

    /* [ init: m <- m + X' * y ] */
    /* FIXME: move to R, at gglm.fit; remove dependency on edgei, edgej */
    for (e = 0; e < ne; e++) { /* for all edges */
      yij = (P->nresponse > 0) ? VECTOR(P->response)[e] : 1.;
      yij *= wij; /* restore counts in binomial case */
      i = edgei[e]; j = edgej[e];
      /* var block */
      for (k = 0; k < P->nvar; k++) {
        xi = P->vardata[k * P->n + i];
        xj = P->vardata[k * P->n + j];
        sijk = (P->varfun[k]) ? xi * xj : xi + xj;
        m[k] += sijk * yij;
      }
      /* fac block */
      c = P->nvar;
      for (k = 0; k < P->nfac; k++) {
        fi = P->facdata[k * P->n + i];
        fj = P->facdata[k * P->n + j];
        for (l = 0; l < P->nlevel[k]; l++) {
          sijk = (P->facfun[k]) ? (fi == l) * (fj == l) :
                                  (fi == l) + (fj == l);
          m[c + l] += sijk * yij;
        }
        c += P->nlevel[k];
      }
    }

    /* [ update ] */
    for (i = 0; i < P->n - 1; i++) {
      s = i * P->n - i * (i + 1) / 2; /* index pair shift */
      for (j = i + 1; j < P->n; j++) {
        /* compute mu and mueta */
        if (initial_stats) {
          e = edgeij[s + j - i - 1];
          if (e != NOEDGE)
            yij = (P->nresponse > 0) ? VECTOR(P->response)[e] : 1.;
          else
            yij = 0.;
          muij = initialize[P->family](yij);
          muetaij = variance[P->family](muij);
          muij -= linkfun[P->family](muij) * muetaij;
        }
        else {
          etaij = compute_etaij(P, i, j, beta, offset, loffset);
          muij = linkinv[P->family](etaij);
          muetaij = variance[P->family](muij);
        }
        muij *= wij; muetaij *= wij;

        /* update mean and variance */
        /* var block */
        for (k = 0; k < P->nvar; k++) {
          xi = P->vardata[k * P->n + i];
          xj = P->vardata[k * P->n + j];
          sijk = (P->varfun[k]) ? xi * xj : xi + xj;
          m[k] -= sijk * muij;

          /* var x var block */
          for (l = k; l < P->nvar; l++) {
            xi = P->vardata[l * P->n + i];
            xj = P->vardata[l * P->n + j];
            sijl = (P->varfun[l]) ? xi * xj : xi + xj;
            V[k + nt * l] += muetaij * sijk * sijl;
          }
          /* var x fac block */
          c = P->nvar;
          for (f = 0; f < P->nfac; f++) {
            fi = P->facdata[f * P->n + i];
            fj = P->facdata[f * P->n + j];
            for (l = 0; l < P->nlevel[f]; l++) {
              sijl = (P->facfun[f]) ? (fi == l) * (fj == l) :
                                      (fi == l) + (fj == l);
              V[k + nt * (c + l)] += muetaij * sijk * sijl;
            }
            c += P->nlevel[f];
          }
        }

        /* fac block */
        c = P->nvar;
        for (k = 0; k < P->nfac; k++) {
          fi = P->facdata[k * P->n + i];
          fj = P->facdata[k * P->n + j];
          for (l = 0; l < P->nlevel[k]; l++) {
            sijk = (P->facfun[k]) ? (fi == l) * (fj == l) :
                                    (fi == l) + (fj == l);
            m[c + l] -= sijk * muij;

            /* fac x fac block */
            d = c;
            for (f = k; f < P->nfac; f++) {
              gi = P->facdata[f * P->n + i];
              gj = P->facdata[f * P->n + j];
              for (r = 0; r < P->nlevel[f]; r++) {
                sijl = (P->facfun[f]) ? (gi == r) * (gj == r) :
                                        (gi == r) + (gj == r);
                V[(c + l) + nt * (d + r)] += muetaij * sijk * sijl;
              }
              d += P->nlevel[f];
            }
          }
          c += P->nlevel[k];
        }
      }
    }
  }

  free(edgei); free(edgej);
  if (initial_stats) free(edgeij);
  return sm;
}
#endif



SEXP lhood_gglmparamset (SEXP sp, SEXP sbeta, SEXP soffset) {
  GGLMParamSet *P = GGLMPARAMSET(sp);
  double *beta = REAL(sbeta);
  double *offset = REAL(soffset);
  int loffset = length(soffset);
  int i, j, e;
  double yij, etaij;
  double lhood = 0;
  double wij = (P->family == 0 && P->wresponse > 0) ? P->wresponse : 1.;
#ifdef _OPENMP
  /*int nthreads = (R_num_math_threads > 0) ? R_num_math_threads : 1;*/
  int nthreads = omp_get_num_procs();
#endif
 
  if (ignore_non_edges) { /* for edges only? */
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) default(none) \
    reduction(+:lhood) \
    private(yij, etaij, i, j) \
    shared(P, wij, ne, cumulant, beta, offset, loffset)
#endif
    for (e = 0; e < ne; e++) {
      yij = (P->nresponse > 0) ? VECTOR(P->response)[e] : 1.;
      igraph_edge(&P->g, e, &i, &j);
      etaij = compute_etaij(P, i, j, beta, offset, loffset);
      lhood += yij * etaij - wij * cumulant[P->family](etaij);
    }

  }
  else { /* for all pairs of vertices (i, j) */
    igraph_vector_t adj;
    int s, k, np = P->n * (P->n - 1) / 2;
    int *isneighbor = (int *) malloc(np * sizeof(int));
    igraph_vector_init(&adj, P->n);

    for (i = 0; i < np; i++) isneighbor[i] = 0; /* init */
    s = 0; /* index pair shift for isneighbor */
    for (i = 0; i < P->n - 1; i++) {
      igraph_neighbors(&P->g, &adj, i, IGRAPH_ALL);
      for (k = 0; k < igraph_vector_size(&adj); k++) {
        j = ((int) VECTOR(adj)[k]) - i - 1; /* relative to i */
        if (j >= 0) isneighbor[s + j] = 1;
      }
      s += P->n - i - 1; /* update shift */
    }
    igraph_vector_destroy(&adj);

#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) default(none) \
    reduction(+:lhood) \
    private(yij, etaij, j, e, s) \
    shared(P, wij, isneighbor, cumulant, beta, offset, loffset)
#endif
    for (i = 0; i < P->n - 1; i++) {
      s = i * P->n - i * (i + 1) / 2; /* index pair shift */
      for (j = i + 1; j < P->n; j++) {
        etaij = compute_etaij(P, i, j, beta, offset, loffset);
        if (isneighbor[s + j - i - 1]) {
          igraph_get_eid(&P->g, &e, i, j, IGRAPH_UNDIRECTED, 1);
          yij = (P->nresponse > 0) ? VECTOR(P->response)[e] : 1.;
        }
        else
          yij = 0.;
        lhood += yij * etaij - wij * cumulant[P->family](etaij);
      }
    }

    free(isneighbor);
  }
  return ScalarReal(lhood);
}

/* TODO: deviance */

SEXP setfactor_gglmparamset (SEXP sp, SEXP sfac, SEXP sz) {
  GGLMParamSet *P = GGLMPARAMSET(sp);
  int *z = INTEGER(sz);
  int *f = P->facdata + asInteger(sfac) * P->n;
  int i;
  for (i = 0; i < P->n; i++) f[i] = z[i];
  return sz;
}

SEXP setvariable_gglmparamset (SEXP sp, SEXP svar, SEXP sz) {
  GGLMParamSet *P = GGLMPARAMSET(sp);
  double *z = REAL(sz);
  double *v = P->vardata + asInteger(svar) * P->n;
  int i;
  for (i = 0; i < P->n; i++) v[i] = z[i];
  return sz;
}

/* compute log-likelihood factor statistics; 'isneighbor' is cache of size
 * P->n; 'l' holds conditional likelihoods: l_ki for label 'k' and node 'i',
 * where k = P->nlevel[fac] */
SEXP labelstats_gglmparamset (SEXP sp, SEXP sfac,
    SEXP sbeta, SEXP soffset, SEXP slogpi, SEXP signore,
    SEXP sisneighbor, SEXP sl) {
  GGLMParamSet *P = GGLMPARAMSET(sp);
  int fac = asInteger(sfac);
  double *beta = REAL(sbeta);
  double *offset = REAL(soffset);
  int loffset = length(soffset);
  double *logpi = REAL(slogpi);
  int ignore_non_edges = asInteger(signore);
  int *isneighbor = INTEGER(sisneighbor);
  double *l = REAL(sl);
  int i, j, k, e;
  int nlevel = P->nlevel[fac];
  int *sigma = P->facdata + fac * P->n;
  igraph_vector_t adj;
  double yij, eta, etaij;
  double *gamma, *li;
  double wij = (P->family == 0 && P->wresponse > 0) ? P->wresponse : 1.;
  
  gamma = beta + P->nvar; /* fac baseline */
  for (k = 0; k < fac; k++) gamma += P->nlevel[k];
  igraph_vector_init(&adj, P->n);

  li = l;
  for (i = 0; i < P->n; i++) {
    /* init */
    for (k = 0; k < nlevel; k++) li[k] = logpi[k];
    igraph_neighbors(&P->g, &adj, i, IGRAPH_ALL);
    for (j = 0; j < P->n; j++) isneighbor[j] = 0;
    for (j = 0; j < igraph_vector_size(&adj); j++)
      isneighbor[(int) VECTOR(adj)[j]] = 1;

    /* compute stats */
    for (j = 0; j < P->n; j++) {
      if (i == j || (ignore_non_edges && !isneighbor[j])) continue;
      eta = compute_etaij_factor(P, i, j, beta, offset, loffset, fac);

      /* [ update l ] */
      yij = (isneighbor[j]) ? 1. : 0.;
      if (P->nresponse > 0 && isneighbor[j]) {
        igraph_get_eid(&P->g, &e, i, j, 0, 0);
        yij = VECTOR(P->response)[e];
      }
      for (k = 0; k < nlevel; k++) {
        etaij = eta;
        if (P->facfun[fac]) /* '*' ? */
          etaij += (sigma[j] == k) * gamma[k];
        else /* '+' */
          etaij += gamma[k] + gamma[sigma[j]];
        li[k] += yij * etaij - wij * cumulant[P->family](etaij);
      }
    }
    li += nlevel;
  }
  igraph_vector_destroy(&adj);
  return sl;
}

/* given a 'k'-by'n' matrix 'l', find sigma[i] = argmax_k l_{ki},
 * i = 1, ..., n, such that card[j] >= nconstraint for each j = 1, ..., k;
 * 'sigma' contains initial configuration */
SEXP labeloptim_gglmparamset (SEXP sk, SEXP sn, SEXP sl,
    SEXP snconstraint, SEXP scard, SEXP ssigma, SEXP slhood) {
  int k = asInteger(sk);
  int n = asInteger(sn);
  double *l = REAL(sl);
  int nconstraint = asInteger(snconstraint);
  int *card = INTEGER(scard);
  int *sigma = INTEGER(ssigma);
  double *lhood = REAL(slhood);
  int i, j, m;
  double u, *li;

  /* init card */
  for (i = 0; i < k; i++) card[i] = 0;
  for (i = 0; i < n; i++) card[sigma[i]]++;

  /* update each position at a time */
  li = l;
  for (i = 0; i < n; i++) {
    if (card[sigma[i]] <= nconstraint) /* cannot change? */
      *lhood += li[sigma[i]];
    else {
      m = 0; u = li[m];
      for (j = 1; j < k; j++)
        if (li[j] > u) { m = j; u = li[m]; }
      if (m != sigma[i]) { /* update card? */
        card[sigma[i]]--;
        card[m]++;
      }
      sigma[i] = m;
      *lhood += u;
    }
    li += k;
  }
  return ssigma;
}


SEXP remap_gglmparamset (SEXP sp, SEXP sfac, SEXP six, SEXP srevind) {
  GGLMParamSet *P = GGLMPARAMSET(sp);
  int fac = asInteger(sfac);
  int *sigma = P->facdata + fac * P->n;
  remap(P->n, P->nlevel[fac], sigma, INTEGER(six), INTEGER(srevind));
  return six;
}


SEXP iter_gglmparamset (SEXP sp, SEXP si, SEXP sj, SEXP sx) {
  GGLMParamSet *P = GGLMPARAMSET(sp);
  int i = asInteger(si);
  int j = asInteger(sj);
  double *x = REAL(sx);
  double xi, xj;
  int fi, fj, c, k;

  for (k = 0; k < P->nvar; k++) {
    xi = P->vardata[k * P->n + i];
    xj = P->vardata[k * P->n + j];
    x[k] = (P->varfun[k]) ? xi * xj : xi + xj;
  }
  c = P->nvar; /* factor column baseline in beta */
  for (k = 0; k < P->nfac; k++) {
    fi = P->facdata[k * P->n + i];
    fj = P->facdata[k * P->n + j];
    if (P->facfun[k]) /* '*' ? */
      x[c + fi] = (fi == fj);
    else /* '+' */
      x[c + fi] = x[c + fj] = 1;
    c += P->nlevel[k];
  }
  return sx;
}




/* [===[ Interface ]===] */

static const R_CallMethodDef callMethods[] = {
  {"sweep", (DL_FUNC) &sweep, 4},
  {"new_gglmparamset", (DL_FUNC) &new_gglmparamset, 12},
  {"setfactor_gglmparamset", (DL_FUNC) &setfactor_gglmparamset, 3},
  {"setvariable_gglmparamset", (DL_FUNC) &setvariable_gglmparamset, 3},
  {"predict_gglmparamset", (DL_FUNC) &predict_gglmparamset, 6},
  {"stats_edges_gglmparamset", (DL_FUNC) &stats_edges_gglmparamset, 6},
  {"stats_pairs_gglmparamset", (DL_FUNC) &stats_pairs_gglmparamset, 6},
  {"lhood_gglmparamset", (DL_FUNC) &lhood_gglmparamset, 4},
  {"labelstats_gglmparamset", (DL_FUNC) &labelstats_gglmparamset, 8},
  {"labeloptim_gglmparamset", (DL_FUNC) &labeloptim_gglmparamset, 7},
  {"remap_gglmparamset", (DL_FUNC) &remap_gglmparamset, 4},
  {"ppl_gglmparamset", (DL_FUNC) &ppl_gglmparamset, 4},
  {"iter_gglmparamset", (DL_FUNC) &iter_gglmparamset, 4},
  {NULL, NULL, 0}
};

void R_init_gglm (DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

