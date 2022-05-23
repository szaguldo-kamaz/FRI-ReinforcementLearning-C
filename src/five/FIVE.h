/*
** Fuzzy Inference by Interpolation in Vague Environment toolbox for ANSI C
**
** https://github.com/szaguldo-kamaz/FRI-ReinforcementLearning-C
**
** ANSI C / AVX port and FixRes optimization by David Vincze <david.vincze@webcode.hu>
** Various contibutions by Daniel Palko <palko.daniel@uni-miskolc.hu>
** mex port by Sandor Zsuga (zsuga@iit.uni-miskolc.hu)
** Original FIVE MATLAB code by Szilveszter Kovacs <szkovacs@iit.uni-miskolc.hu>
*/

///
///\file FIVE.h
///

#ifndef FIVE_H
#define FIVE_H

#define FIVE_MAX_NUM_OF_UNIVERSES  8

#include "config.h"

///common FIVE struct
struct FIVERB {
	double *u;         ///<all universe vectors
	double *ve;        ///<all vague environment vectors
//	double *rb;        ///<rule-base (all rule vectors)
	double *psc;       ///<scaling points
	double *scf;       ///<scaling functions
// mu volt
	int numofunivs;    ///<dimension count of the universe array (number of antecedents + conclusion in case it is also fuzzy, usually onbly no., of antec)
// nu volt
	int univlength;    ///<length of one universe (resolution of one dimension, number of elements in one "u")
// mr volt
	int numofrules;    ///<number of rules (number of vectors in *rant and *rconc)
	int maxnumofrules;  ///<maximum number of rules
// nr volt
	int rulelength;    ///<number of rule dimensions (number of elements in each vector (one rule) in *rb)
	int numofantecedents; ///<number of rule antecedents
	int p;             ///<Shepard interpolation parameter

	// cache, precalc
	double *valvagp;   ///<ValVag - value of VE point
	double *valvagu;   ///<ValVag - pointer to the consequent dim U
	double *valvagve;  ///<ValVag - pointer to the consequent dim VE
	double valvagdims; ///<number of dimensions used in ValVag - normally it is only used for the consequent dim if it has a VE - so normally this should be 1
	double *ruledists; ///<rule distances
	double *rant;      ///<rule-base antecedents (without conclusions) - values for each dimension following each other sequentially
	double **rseqant;  ///<array of rule-base antecedents (without conclusions) where values are stored in sequential order for each dimension
	double *ract;      ///<a pointer to the action dimension in rseqant (rseqant[last]) - for FRIRL
	unsigned int *rant_uindex;      ///<rule-base antecedents (without conclusions) - values for each dimension following each other sequentially - aligned to possible points in U, index in U is stored
	unsigned int *rant_veval;       ///<VE values of rule-base antecedents (without conclusions) - values for each dimension following each other sequentially - aligned to possible points in U, VE value for index in U is stored
	unsigned int **rseqant_uindex;  ///<array of rule-base antecedents (without conclusions) where values are stored in sequential order for each dimension - aligned to possible points in U, index in U is stored
	double **rseqant_veval;         ///<array of VE values of rule-base antecedents (without conclusions) where values are stored in sequential order for each dimension - aligned to possible points in U, VE value for index in U is stored
	unsigned int *ract_uindex;      ///<a pointer to the action dimension in rseqant_uindex (rseqant_uindex[last]) - for FRIRL
	double *ract_veval;             ///<a pointer to the action dimension in rseqant_veval (rseqant_veval[last]) - for FRIRL
	double *rconc;       ///<rule-base consequents (without antecedents)
	double *weights;     ///<VagConclWeight result for all rules
	unsigned int uksize; ///<resolution of universe
	double *ukdomains  ; ///<lengths of the universe domains: u[last] - u[first]
	double *udivs;       ///<step sizes of the universes ukdomain/univlength ( univlength=size(u) )
	double **uk;         ///<shortcut pointers to u by dims (uk[0..k])
	double **vek;        ///<shortcut pointers to ve by dims (vek[0..k])
	double *wi;          ///<to be used in VagConclWeight
	double *frd_dists;   ///<to be used in RuleDist
	double *fvc_vagdist; ///<to be used in VagConcl
	// TODO RB removed
	//	double *newrule;     ///<to be used in addrule(), points to the rule place after the last rule
	double *newrant;     ///<to be used in addrule(), points to the rule ant place after the last rule ant
	double *newrconc;
// TODO remove debug
unsigned int epno; // temp
#ifdef BUILD_AVX2
	unsigned int avx2_rbsize; ///<numofrules/4 + ((numofrules mod 4) > 0)
#endif
};

// API rc5 (obsoleted, kept for backward compatibility)
struct FIVERB *FIVEInit(double *u, double *ve, int p, int numofunivs, int univlength, int numofrules, int maxnumofrules, int rulelength, double *rant, double *rconc);

double *FIVEGScFunc(double *u, int numofunivs, int univlength, double *psc, int mp, int np, double nls);
double *FIVEGVagEnv(double *u, int numofunivs, int univlength, double *scf);

int FIVEValVag(struct FIVERB *frb, double *vp);

double FIVEVagConcl(struct FIVERB *frb, double *x);
unsigned int FIVEVagConclWeight(struct FIVERB *frb, double *x);
double FIVEVagConcl_FRIRL_BestAct(struct FIVERB *frb, double *ruledists);
int     FIVEAddRule(struct FIVERB *frb, double *newrule);

// old API (obsoleted, to be removed)
int five_vague_distance(struct FIVERB *frb, fri_float *p1, fri_float *p2, fri_float *d);
int five_vague_distance_parallel(struct FIVERB *frb, fri_float *p1, int p1_offset, fri_float *p2, fri_float *d);
int five_rule_distance(struct FIVERB *frb, fri_float *x);
int five_add_rule(struct FIVERB *frb, fri_float *ruletoadd);
int five_remove_rule(struct FIVERB *frb, unsigned int rulenotoremove);
void five_deinit(struct FIVERB *frb);

// API (current)
int FIVE_add_rule(struct FIVERB *frb, fri_float *rant, fri_float rconc);
unsigned int FIVE_vag_concl_weight(struct FIVERB *frb, double *ant, double *weights);
unsigned int FIVE_vag_concl(struct FIVERB *frb, double *ant, double *conc);
int FIVE_GSc_func(double *u, int numofunivs, int univlength, double *psc, int mp, int np, double nls, double *scf);

#endif //FIVE_H
