#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
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


static double compute_etaij (GGLMParamSet *P, int i, int j,
    double *beta, double *offset, int loffset) {
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
    fi = P->facdata[k * P->n + i];
    fj = P->facdata[k * P->n + j];
    if (P->facfun[k]) /* '*' ? */
      etaij += (fi == fj) * beta[c + fi]; /*facfun 0 or 1, change this if we include other options.*/
    else /* '+' */
      etaij += beta[c + fi] + beta[c + fj];
    c += P->nlevel[k];
  }
  /*if (loffset > 0)
    etaij += (loffset == 1) ? (*offset) : offset[i * P->n + j];
  return etaij;
  */
  if (loffset > 0) {
    if (loffset == 1) etaij += *offset;
    else {
      int s = i * P->n - i * (i + 1) / 2; /* index pair shift */
      etaij += offset[s + j - i - 1];
    }
  }
  return etaij;
}


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
    SEXP seta) {
  GGLMParamSet *P = GGLMPARAMSET(sp);
  double *beta = REAL(sbeta);
  double *offset = REAL(soffset);
  int type = asInteger(stype);
  double *eta = REAL(seta);
  int loffset = length(soffset);
  int i, j;
  double etaij;

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
  return seta;
}


SEXP stats_pairs_gglmparamset (SEXP sp, SEXP sbeta, SEXP soffset, SEXP sinit,
    SEXP sm, SEXP sV) {
  GGLMParamSet *P = GGLMPARAMSET(sp); /*defines pointer*/
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


SEXP lhood_gglmparamset (SEXP sp, SEXP sbeta, SEXP soffset) {
  GGLMParamSet *P = GGLMPARAMSET(sp);
  double *beta = REAL(sbeta);
  double *offset = REAL(soffset);
  int loffset = length(soffset);
  int i, j, e;
  double yij, etaij;
  double lhood = 0;
  double wij = (P->family == 0 && P->wresponse > 0) ? P->wresponse : 1.;
 
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
  return ScalarReal(lhood);
}


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




/* [===[ Interface ]===] */

static const R_CallMethodDef callMethods[] = {
  {"new_gglmparamset", (DL_FUNC) &new_gglmparamset, 12},
  {"setfactor_gglmparamset", (DL_FUNC) &setfactor_gglmparamset, 3},
  {"setvariable_gglmparamset", (DL_FUNC) &setvariable_gglmparamset, 3},
  {"predict_gglmparamset", (DL_FUNC) &predict_gglmparamset, 5},
  {"stats_pairs_gglmparamset", (DL_FUNC) &stats_pairs_gglmparamset, 6},
  {"lhood_gglmparamset", (DL_FUNC) &lhood_gglmparamset, 3},
  {NULL, NULL, 0}
};

void R_init_gglm (DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

