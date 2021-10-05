/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__Ca
#define _nrn_initial _nrn_initial__Ca
#define nrn_cur _nrn_cur__Ca
#define _nrn_current _nrn_current__Ca
#define nrn_jacob _nrn_jacob__Ca
#define nrn_state _nrn_state__Ca
#define _net_receive _net_receive__Ca 
#define state state__Ca 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gtcabar _p[0]
#define gncabar _p[1]
#define glcabar _p[2]
#define e_ca _p[3]
#define gtca _p[4]
#define gnca _p[5]
#define glca _p[6]
#define ca_i _p[7]
#define a _p[8]
#define b _p[9]
#define c _p[10]
#define d _p[11]
#define e _p[12]
#define ica _p[13]
#define Dca_i _p[14]
#define Da _p[15]
#define Db _p[16]
#define Dc _p[17]
#define Dd _p[18]
#define De _p[19]
#define v _p[20]
#define _g _p[21]
#define _ion_ica	*_ppvar[0]._pval
#define _ion_dicadv	*_ppvar[1]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static void _hoc_alphae(void);
 static void _hoc_alphad(void);
 static void _hoc_alphac(void);
 static void _hoc_alphab(void);
 static void _hoc_alphaa(void);
 static void _hoc_betae(void);
 static void _hoc_betad(void);
 static void _hoc_betac(void);
 static void _hoc_betab(void);
 static void _hoc_betaa(void);
 static void _hoc_exponential(void);
 static void _hoc_f(void);
 static void _hoc_logistic(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_Ca", _hoc_setdata,
 "alphae_Ca", _hoc_alphae,
 "alphad_Ca", _hoc_alphad,
 "alphac_Ca", _hoc_alphac,
 "alphab_Ca", _hoc_alphab,
 "alphaa_Ca", _hoc_alphaa,
 "betae_Ca", _hoc_betae,
 "betad_Ca", _hoc_betad,
 "betac_Ca", _hoc_betac,
 "betab_Ca", _hoc_betab,
 "betaa_Ca", _hoc_betaa,
 "exponential_Ca", _hoc_exponential,
 "f_Ca", _hoc_f,
 "logistic_Ca", _hoc_logistic,
 0, 0
};
#define alphae alphae_Ca
#define alphad alphad_Ca
#define alphac alphac_Ca
#define alphab alphab_Ca
#define alphaa alphaa_Ca
#define betae betae_Ca
#define betad betad_Ca
#define betac betac_Ca
#define betab betab_Ca
#define betaa betaa_Ca
#define exponential exponential_Ca
#define f f_Ca
#define logistic logistic_Ca
 extern double alphae( _threadargsprotocomma_ double );
 extern double alphad( _threadargsprotocomma_ double );
 extern double alphac( _threadargsprotocomma_ double );
 extern double alphab( _threadargsprotocomma_ double );
 extern double alphaa( _threadargsprotocomma_ double );
 extern double betae( _threadargsprotocomma_ double );
 extern double betad( _threadargsprotocomma_ double );
 extern double betac( _threadargsprotocomma_ double );
 extern double betab( _threadargsprotocomma_ double );
 extern double betaa( _threadargsprotocomma_ double );
 extern double exponential( _threadargsprotocomma_ double , double , double , double );
 extern double f( _threadargsprotocomma_ double , double , double , double );
 extern double logistic( _threadargsprotocomma_ double , double , double , double );
 /* declare global and static user variables */
#define cao cao_Ca
 double cao = 1.2;
#define ca0 ca0_Ca
 double ca0 = 7e-005;
#define tauc tauc_Ca
 double tauc = 9;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "ca0_Ca", "mM",
 "cao_Ca", "mM",
 "tauc_Ca", "ms",
 "gtcabar_Ca", "S/cm2",
 "gncabar_Ca", "S/cm2",
 "glcabar_Ca", "S/cm2",
 "ca_i_Ca", "mM",
 "e_ca_Ca", "mV",
 "gtca_Ca", "S/cm2",
 "gnca_Ca", "S/cm2",
 "glca_Ca", "S/cm2",
 0,0
};
 static double a0 = 0;
 static double b0 = 0;
 static double c0 = 0;
 static double ca_i0 = 0;
 static double delta_t = 0.01;
 static double d0 = 0;
 static double e0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "ca0_Ca", &ca0_Ca,
 "cao_Ca", &cao_Ca,
 "tauc_Ca", &tauc_Ca,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[2]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"Ca",
 "gtcabar_Ca",
 "gncabar_Ca",
 "glcabar_Ca",
 0,
 "e_ca_Ca",
 "gtca_Ca",
 "gnca_Ca",
 "glca_Ca",
 0,
 "ca_i_Ca",
 "a_Ca",
 "b_Ca",
 "c_Ca",
 "d_Ca",
 "e_Ca",
 0,
 0};
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 22, _prop);
 	/*initialize range parameters*/
 	gtcabar = 0;
 	gncabar = 0;
 	glcabar = 0;
 	_prop->param = _p;
 	_prop->param_size = 22;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 	_ppvar[0]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[1]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _ca_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("ca", -10000.);
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 22, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 Ca C:/Users/Viktor/Desktop/Emory/PV_IN/Bartos2002/ca.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double B = .26;
 static double F = 96485.3;
 static double R = 8.3145;
 static double TEMP = 25.0;
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[6], _dlist1[6];
 static int state(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   Dca_i = - B * ica - ( ca_i - ca0 ) / tauc ;
   Da = alphaa ( _threadargscomma_ v ) * ( 1.0 - a ) - betaa ( _threadargscomma_ v ) * a ;
   Db = alphab ( _threadargscomma_ v ) * ( 1.0 - b ) - betab ( _threadargscomma_ v ) * b ;
   Dc = alphac ( _threadargscomma_ v ) * ( 1.0 - c ) - betac ( _threadargscomma_ v ) * c ;
   Dd = alphad ( _threadargscomma_ v ) * ( 1.0 - d ) - betad ( _threadargscomma_ v ) * d ;
   De = alphae ( _threadargscomma_ v ) * ( 1.0 - e ) - betae ( _threadargscomma_ v ) * e ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 Dca_i = Dca_i  / (1. - dt*( ( - ( ( 1.0 ) ) / tauc ) )) ;
 Da = Da  / (1. - dt*( ( alphaa ( _threadargscomma_ v ) )*( ( ( - 1.0 ) ) ) - ( betaa ( _threadargscomma_ v ) )*( 1.0 ) )) ;
 Db = Db  / (1. - dt*( ( alphab ( _threadargscomma_ v ) )*( ( ( - 1.0 ) ) ) - ( betab ( _threadargscomma_ v ) )*( 1.0 ) )) ;
 Dc = Dc  / (1. - dt*( ( alphac ( _threadargscomma_ v ) )*( ( ( - 1.0 ) ) ) - ( betac ( _threadargscomma_ v ) )*( 1.0 ) )) ;
 Dd = Dd  / (1. - dt*( ( alphad ( _threadargscomma_ v ) )*( ( ( - 1.0 ) ) ) - ( betad ( _threadargscomma_ v ) )*( 1.0 ) )) ;
 De = De  / (1. - dt*( ( alphae ( _threadargscomma_ v ) )*( ( ( - 1.0 ) ) ) - ( betae ( _threadargscomma_ v ) )*( 1.0 ) )) ;
  return 0;
}
 /*END CVODE*/
 static int state (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
    ca_i = ca_i + (1. - exp(dt*(( - ( ( 1.0 ) ) / tauc ))))*(- ( ( - B )*( ica ) - ( ( ( - ca0 ) ) ) / tauc ) / ( ( - ( ( 1.0 ) ) / tauc ) ) - ca_i) ;
    a = a + (1. - exp(dt*(( alphaa ( _threadargscomma_ v ) )*( ( ( - 1.0 ) ) ) - ( betaa ( _threadargscomma_ v ) )*( 1.0 ))))*(- ( ( alphaa ( _threadargscomma_ v ) )*( ( 1.0 ) ) ) / ( ( alphaa ( _threadargscomma_ v ) )*( ( ( - 1.0 ) ) ) - ( betaa ( _threadargscomma_ v ) )*( 1.0 ) ) - a) ;
    b = b + (1. - exp(dt*(( alphab ( _threadargscomma_ v ) )*( ( ( - 1.0 ) ) ) - ( betab ( _threadargscomma_ v ) )*( 1.0 ))))*(- ( ( alphab ( _threadargscomma_ v ) )*( ( 1.0 ) ) ) / ( ( alphab ( _threadargscomma_ v ) )*( ( ( - 1.0 ) ) ) - ( betab ( _threadargscomma_ v ) )*( 1.0 ) ) - b) ;
    c = c + (1. - exp(dt*(( alphac ( _threadargscomma_ v ) )*( ( ( - 1.0 ) ) ) - ( betac ( _threadargscomma_ v ) )*( 1.0 ))))*(- ( ( alphac ( _threadargscomma_ v ) )*( ( 1.0 ) ) ) / ( ( alphac ( _threadargscomma_ v ) )*( ( ( - 1.0 ) ) ) - ( betac ( _threadargscomma_ v ) )*( 1.0 ) ) - c) ;
    d = d + (1. - exp(dt*(( alphad ( _threadargscomma_ v ) )*( ( ( - 1.0 ) ) ) - ( betad ( _threadargscomma_ v ) )*( 1.0 ))))*(- ( ( alphad ( _threadargscomma_ v ) )*( ( 1.0 ) ) ) / ( ( alphad ( _threadargscomma_ v ) )*( ( ( - 1.0 ) ) ) - ( betad ( _threadargscomma_ v ) )*( 1.0 ) ) - d) ;
    e = e + (1. - exp(dt*(( alphae ( _threadargscomma_ v ) )*( ( ( - 1.0 ) ) ) - ( betae ( _threadargscomma_ v ) )*( 1.0 ))))*(- ( ( alphae ( _threadargscomma_ v ) )*( ( 1.0 ) ) ) / ( ( alphae ( _threadargscomma_ v ) )*( ( ( - 1.0 ) ) ) - ( betae ( _threadargscomma_ v ) )*( 1.0 ) ) - e) ;
   }
  return 0;
}
 
double alphaa ( _threadargsprotocomma_ double _lv ) {
   double _lalphaa;
 _lalphaa = f ( _threadargscomma_ 2.0 , 0.1 , _lv , 19.26 ) ;
   
return _lalphaa;
 }
 
static void _hoc_alphaa(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  alphaa ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double betaa ( _threadargsprotocomma_ double _lv ) {
   double _lbetaa;
 _lbetaa = exponential ( _threadargscomma_ 0.009 , - 0.045393 , _lv , 0.0 ) ;
   
return _lbetaa;
 }
 
static void _hoc_betaa(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  betaa ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double alphab ( _threadargsprotocomma_ double _lv ) {
   double _lalphab;
 _lalphab = exponential ( _threadargscomma_ 1e-6 , - 0.061501 , _lv , 0.0 ) ;
   
return _lalphab;
 }
 
static void _hoc_alphab(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  alphab ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double betab ( _threadargsprotocomma_ double _lv ) {
   double _lbetab;
 _lbetab = logistic ( _threadargscomma_ 1.0 , - 0.1 , _lv , 29.79 ) ;
   
return _lbetab;
 }
 
static void _hoc_betab(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  betab ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double alphac ( _threadargsprotocomma_ double _lv ) {
   double _lalphac;
 _lalphac = f ( _threadargscomma_ 1.9 , 0.1 , _lv , 19.88 ) ;
   
return _lalphac;
 }
 
static void _hoc_alphac(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  alphac ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double betac ( _threadargsprotocomma_ double _lv ) {
   double _lbetac;
 _lbetac = exponential ( _threadargscomma_ 0.046 , - 0.048239 , _lv , 0.0 ) ;
   
return _lbetac;
 }
 
static void _hoc_betac(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  betac ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double alphad ( _threadargsprotocomma_ double _lv ) {
   double _lalphad;
 _lalphad = exponential ( _threadargscomma_ 1.6e-4 , - 0.020661 , _lv , 0.0 ) ;
   
return _lalphad;
 }
 
static void _hoc_alphad(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  alphad ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double betad ( _threadargsprotocomma_ double _lv ) {
   double _lbetad;
 _lbetad = logistic ( _threadargscomma_ 1.0 , - 0.1 , _lv , 39.0 ) ;
   
return _lbetad;
 }
 
static void _hoc_betad(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  betad ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double alphae ( _threadargsprotocomma_ double _lv ) {
   double _lalphae;
 _lalphae = f ( _threadargscomma_ 156.9 , 0.1 , _lv , 81.5 ) ;
   
return _lalphae;
 }
 
static void _hoc_alphae(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  alphae ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double betae ( _threadargsprotocomma_ double _lv ) {
   double _lbetae;
 _lbetae = exponential ( _threadargscomma_ 0.29 , - 0.092081 , _lv , 0.0 ) ;
   
return _lbetae;
 }
 
static void _hoc_betae(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  betae ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double f ( _threadargsprotocomma_ double _lA , double _lk , double _lv , double _lD ) {
   double _lf;
 double _lx ;
  _lx = _lk * ( _lv - _lD ) ;
   if ( fabs ( _lx ) > 1e-6 ) {
     _lf = _lA * _lx / ( 1.0 - exp ( - _lx ) ) ;
     }
   else {
     _lf = _lA / ( 1.0 - 0.5 * _lx ) ;
     }
    
return _lf;
 }
 
static void _hoc_f(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  f ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) );
 hoc_retpushx(_r);
}
 
double logistic ( _threadargsprotocomma_ double _lA , double _lk , double _lv , double _lD ) {
   double _llogistic;
  _llogistic = _lA / ( 1.0 + exp ( _lk * ( _lv - _lD ) ) ) ;
    
return _llogistic;
 }
 
static void _hoc_logistic(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  logistic ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) );
 hoc_retpushx(_r);
}
 
double exponential ( _threadargsprotocomma_ double _lA , double _lk , double _lv , double _lD ) {
   double _lexponential;
  _lexponential = _lA * exp ( _lk * ( _lv - _lD ) ) ;
    
return _lexponential;
 }
 
static void _hoc_exponential(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  exponential ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 6;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 6; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  a = a0;
  b = b0;
  c = c0;
  ca_i = ca_i0;
  d = d0;
  e = e0;
 {
   ca_i = ca0 ;
   a = alphaa ( _threadargscomma_ v ) / ( alphaa ( _threadargscomma_ v ) + betaa ( _threadargscomma_ v ) ) ;
   b = alphab ( _threadargscomma_ v ) / ( alphab ( _threadargscomma_ v ) + betab ( _threadargscomma_ v ) ) ;
   c = alphac ( _threadargscomma_ v ) / ( alphac ( _threadargscomma_ v ) + betac ( _threadargscomma_ v ) ) ;
   d = alphad ( _threadargscomma_ v ) / ( alphad ( _threadargscomma_ v ) + betad ( _threadargscomma_ v ) ) ;
   e = alphae ( _threadargscomma_ v ) / ( alphae ( _threadargscomma_ v ) + betae ( _threadargscomma_ v ) ) ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   e_ca = ( 1000.0 ) * ( TEMP + 273.15 ) * R / ( 2.0 * F ) * log ( cao / ca_i ) ;
   gtca = gtcabar * a * a * b ;
   gnca = gncabar * c * c * d ;
   glca = glcabar * e * e ;
   ica = ( gtca + gnca + glca ) * ( v - e_ca ) ;
   }
 _current += ica;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dica;
  _dica = ica;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dicadv += (_dica - ica)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ica += ica ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 {   state(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(ca_i) - _p;  _dlist1[0] = &(Dca_i) - _p;
 _slist1[1] = &(a) - _p;  _dlist1[1] = &(Da) - _p;
 _slist1[2] = &(b) - _p;  _dlist1[2] = &(Db) - _p;
 _slist1[3] = &(c) - _p;  _dlist1[3] = &(Dc) - _p;
 _slist1[4] = &(d) - _p;  _dlist1[4] = &(Dd) - _p;
 _slist1[5] = &(e) - _p;  _dlist1[5] = &(De) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "ca.mod";
static const char* nmodl_file_text = 
  ": Ca channels (T,N,L-type)\n"
  "\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX Ca\n"
  "	USEION ca WRITE ica\n"
  "	RANGE gtcabar, gncabar, glcabar, gtca, gnca, glca, e_ca\n"
  "	GLOBAL ca0, cao, tauc\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(molar) = (1/liter)\n"
  "	(mM) = (millimolar)\n"
  "	(mV) = (millivolt)\n"
  "	(mA) = (milliamp)\n"
  "	(S) = (siemens)\n"
  "	B = .26 (mM-cm2/mA-ms)\n"
  "	F = (faraday) (coulomb)\n"
  "	R = (k-mole) (joule/degC)\n"
  "	TEMP = 25 (degC)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	ca0 = .00007	(mM)		: initial calcium concentration inside\n"
  "	cao = 1.2		(mM)		: calcium concentration outside\n"
  "	tauc = 9		(ms)\n"
  "	gtcabar = .00	(S/cm2)	: maximum permeability\n"
  "	gncabar = .00	(S/cm2)\n"
  "	glcabar = .00	(S/cm2)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v		(mV)\n"
  "	e_ca	(mV)\n"
  "	ica		(mA/cm2)\n"
  "	gtca		(S/cm2)\n"
  "	gnca		(S/cm2)\n"
  "	glca		(S/cm2)\n"
  "}\n"
  "\n"
  "STATE { ca_i (mM) a b c d e}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE state METHOD cnexp\n"
  "	e_ca = (1000)*(TEMP+273.15)*R/(2*F)*log(cao/ca_i)\n"
  "	gtca = gtcabar*a*a*b\n"
  "	gnca = gncabar*c*c*d\n"
  "	glca = glcabar*e*e\n"
  "	ica = (gtca+gnca+glca)*(v - e_ca)\n"
  "}\n"
  "\n"
  "DERIVATIVE state {	: exact when v held constant; integrates over dt step\n"
  "	ca_i' = -B*ica-(ca_i-ca0)/tauc\n"
  "	a' = alphaa(v)*(1-a)-betaa(v)*a\n"
  "	b' = alphab(v)*(1-b)-betab(v)*b\n"
  "	c' = alphac(v)*(1-c)-betac(v)*c\n"
  "	d' = alphad(v)*(1-d)-betad(v)*d\n"
  "	e' = alphae(v)*(1-e)-betae(v)*e\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	ca_i = ca0\n"
  "	a = alphaa(v)/(alphaa(v)+betaa(v))\n"
  "	b = alphab(v)/(alphab(v)+betab(v))\n"
  "	c = alphac(v)/(alphac(v)+betac(v))\n"
  "	d = alphad(v)/(alphad(v)+betad(v))\n"
  "	e = alphae(v)/(alphae(v)+betae(v))\n"
  "}\n"
  "\n"
  "FUNCTION alphaa(v (mV)) (/ms) {\n"
  "	alphaa = f(2,0.1,v,19.26)\n"
  "}\n"
  "\n"
  "FUNCTION betaa(v (mV)) (/ms) {\n"
  "	betaa = exponential(0.009,-0.045393,v,0)\n"
  "}\n"
  "\n"
  "FUNCTION alphab(v (mV)) (/ms) {\n"
  "	alphab = exponential(1e-6,-0.061501,v,0)\n"
  "}\n"
  "\n"
  "FUNCTION betab(v (mV)) (/ms) {\n"
  "	betab = logistic(1,-0.1,v,29.79)\n"
  "}\n"
  "\n"
  "FUNCTION alphac(v (mV)) (/ms) {\n"
  "	alphac = f(1.9,0.1,v,19.88)\n"
  "}\n"
  "\n"
  "FUNCTION betac(v (mV)) (/ms) {\n"
  "	betac = exponential(0.046,-0.048239,v,0)\n"
  "}\n"
  "\n"
  "FUNCTION alphad(v (mV)) (/ms) {\n"
  "	alphad = exponential(1.6e-4,-0.020661,v,0)\n"
  "}\n"
  "\n"
  "FUNCTION betad(v (mV)) (/ms) {\n"
  "	betad = logistic(1,-0.1,v,39)\n"
  "}\n"
  "\n"
  "FUNCTION alphae(v (mV)) (/ms) {\n"
  "	alphae = f(156.9,0.1,v,81.5)\n"
  "}\n"
  "\n"
  "FUNCTION betae(v (mV)) (/ms) {\n"
  "	betae = exponential(0.29,-0.092081,v,0)\n"
  "}\n"
  "\n"
  "FUNCTION f(A, k, v (mV), D) (/ms) {\n"
  "	LOCAL x\n"
  "	UNITSOFF\n"
  "	x = k*(v-D)\n"
  "	if (fabs(x) > 1e-6) {\n"
  "		f = A*x/(1-exp(-x))\n"
  "	}else{\n"
  "		f = A/(1-0.5*x)\n"
  "	}\n"
  "	UNITSON\n"
  "}\n"
  "\n"
  "FUNCTION logistic(A, k, v (mV), D) (/ms) {\n"
  "	UNITSOFF\n"
  "	logistic = A/(1+exp(k*(v-D)))\n"
  "	UNITSON\n"
  "}\n"
  "\n"
  "FUNCTION exponential(A, k, v (mV), D) (/ms) {\n"
  "	UNITSOFF\n"
  "	exponential = A*exp(k*(v-D))\n"
  "	UNITSON\n"
  "}\n"
  ;
#endif
