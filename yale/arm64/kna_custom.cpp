/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#define _pval pval
// clang-format off
#include "md1redef.h"
#include "section_fwd.hpp"
#include "nrniv_mf.h"
#include "md2redef.h"
#include "nrnconf.h"
// clang-format on
#include "neuron/cache/mechanism_range.hpp"
static constexpr auto number_of_datum_variables = 6;
static constexpr auto number_of_floating_point_variables = 22;
namespace {
template <typename T>
using _nrn_mechanism_std_vector = std::vector<T>;
using _nrn_model_sorted_token = neuron::model_sorted_token;
using _nrn_mechanism_cache_range = neuron::cache::MechanismRange<number_of_floating_point_variables, number_of_datum_variables>;
using _nrn_mechanism_cache_instance = neuron::cache::MechanismInstance<number_of_floating_point_variables, number_of_datum_variables>;
using _nrn_non_owning_id_without_container = neuron::container::non_owning_identifier_without_container;
template <typename T>
using _nrn_mechanism_field = neuron::mechanism::field<T>;
template <typename... Args>
void _nrn_mechanism_register_data_fields(Args&&... args) {
  neuron::mechanism::register_data_fields(std::forward<Args>(args)...);
}
}
 
#if !NRNGPU
#undef exp
#define exp hoc_Exp
#if NRN_ENABLE_ARCH_INDEP_EXP_POW
#undef pow
#define pow hoc_pow
#endif
#endif
 
#define nrn_init _nrn_init__kna_custom
#define _nrn_initial _nrn_initial__kna_custom
#define nrn_cur _nrn_cur__kna_custom
#define _nrn_current _nrn_current__kna_custom
#define nrn_jacob _nrn_jacob__kna_custom
#define nrn_state _nrn_state__kna_custom
#define _net_receive _net_receive__kna_custom 
#define rates rates__kna_custom 
#define states states__kna_custom 
 
#define _threadargscomma_ _ml, _iml, _ppvar, _thread, _globals, _nt,
#define _threadargsprotocomma_ Memb_list* _ml, size_t _iml, Datum* _ppvar, Datum* _thread, double* _globals, NrnThread* _nt,
#define _internalthreadargsprotocomma_ _nrn_mechanism_cache_range* _ml, size_t _iml, Datum* _ppvar, Datum* _thread, double* _globals, NrnThread* _nt,
#define _threadargs_ _ml, _iml, _ppvar, _thread, _globals, _nt
#define _threadargsproto_ Memb_list* _ml, size_t _iml, Datum* _ppvar, Datum* _thread, double* _globals, NrnThread* _nt
#define _internalthreadargsproto_ _nrn_mechanism_cache_range* _ml, size_t _iml, Datum* _ppvar, Datum* _thread, double* _globals, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *hoc_getarg(int);
 
#define t _nt->_t
#define dt _nt->_dt
#define gkbar _ml->template fpfield<0>(_iml)
#define gkbar_columnindex 0
#define gnabar _ml->template fpfield<1>(_iml)
#define gnabar_columnindex 1
#define gk _ml->template fpfield<2>(_iml)
#define gk_columnindex 2
#define gna _ml->template fpfield<3>(_iml)
#define gna_columnindex 3
#define n _ml->template fpfield<4>(_iml)
#define n_columnindex 4
#define m _ml->template fpfield<5>(_iml)
#define m_columnindex 5
#define h _ml->template fpfield<6>(_iml)
#define h_columnindex 6
#define Dn _ml->template fpfield<7>(_iml)
#define Dn_columnindex 7
#define Dm _ml->template fpfield<8>(_iml)
#define Dm_columnindex 8
#define Dh _ml->template fpfield<9>(_iml)
#define Dh_columnindex 9
#define ek _ml->template fpfield<10>(_iml)
#define ek_columnindex 10
#define ena _ml->template fpfield<11>(_iml)
#define ena_columnindex 11
#define ik _ml->template fpfield<12>(_iml)
#define ik_columnindex 12
#define ina _ml->template fpfield<13>(_iml)
#define ina_columnindex 13
#define alpha_n _ml->template fpfield<14>(_iml)
#define alpha_n_columnindex 14
#define beta_n _ml->template fpfield<15>(_iml)
#define beta_n_columnindex 15
#define alpha_m _ml->template fpfield<16>(_iml)
#define alpha_m_columnindex 16
#define beta_m _ml->template fpfield<17>(_iml)
#define beta_m_columnindex 17
#define alpha_h _ml->template fpfield<18>(_iml)
#define alpha_h_columnindex 18
#define beta_h _ml->template fpfield<19>(_iml)
#define beta_h_columnindex 19
#define v _ml->template fpfield<20>(_iml)
#define v_columnindex 20
#define _g _ml->template fpfield<21>(_iml)
#define _g_columnindex 21
#define _ion_ek *(_ml->dptr_field<0>(_iml))
#define _p_ion_ek static_cast<neuron::container::data_handle<double>>(_ppvar[0])
#define _ion_ik *(_ml->dptr_field<1>(_iml))
#define _p_ion_ik static_cast<neuron::container::data_handle<double>>(_ppvar[1])
#define _ion_dikdv *(_ml->dptr_field<2>(_iml))
#define _ion_ena *(_ml->dptr_field<3>(_iml))
#define _p_ion_ena static_cast<neuron::container::data_handle<double>>(_ppvar[3])
#define _ion_ina *(_ml->dptr_field<4>(_iml))
#define _p_ion_ina static_cast<neuron::container::data_handle<double>>(_ppvar[4])
#define _ion_dinadv *(_ml->dptr_field<5>(_iml))
 /* Thread safe. No static _ml, _iml or _ppvar. */
 static int hoc_nrnpointerindex =  -1;
 static _nrn_mechanism_std_vector<Datum> _extcall_thread;
 static Prop* _extcall_prop;
 /* _prop_id kind of shadows _extcall_prop to allow validity checking. */
 static _nrn_non_owning_id_without_container _prop_id{};
 /* external NEURON variables */
 /* declaration of user functions */
 static void _hoc_rates(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mechtype);
#endif
 static void _hoc_setdata();
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 {"setdata_kna_custom", _hoc_setdata},
 {"rates_kna_custom", _hoc_rates},
 {0, 0}
};
 
/* Direct Python call wrappers to density mechanism functions.*/
 static double _npy_rates(Prop*);
 
static NPyDirectMechFunc npy_direct_func_proc[] = {
 {"rates", _npy_rates},
 {0, 0}
};
 /* declare global and static user variables */
 #define gind 0
 #define _gth 0
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 {0, 0, 0}
};
 static HocParmUnits _hoc_parm_units[] = {
 {"gkbar_kna_custom", "S/cm2"},
 {"gnabar_kna_custom", "S/cm2"},
 {"gk_kna_custom", "S/cm2"},
 {"gna_kna_custom", "S/cm2"},
 {0, 0}
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 static double n0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 {0, 0}
};
 static DoubVec hoc_vdoub[] = {
 {0, 0, 0}
};
 static double _sav_indep;
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 _prop_id = _nrn_get_prop_id(_prop);
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 static void nrn_alloc(Prop*);
static void nrn_init(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
static void nrn_state(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
 static void nrn_cur(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
static void nrn_jacob(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(Prop*, int, neuron::container::data_handle<double>*, neuron::container::data_handle<double>*, double*, int);
static void _ode_spec(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
static void _ode_matsol(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
 
#define _cvode_ieq _ppvar[6].literal_value<int>()
 static void _ode_matsol_instance1(_internalthreadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"kna_custom",
 "gkbar_kna_custom",
 "gnabar_kna_custom",
 0,
 "gk_kna_custom",
 "gna_kna_custom",
 0,
 "n_kna_custom",
 "m_kna_custom",
 "h_kna_custom",
 0,
 0};
 static Symbol* _k_sym;
 static Symbol* _na_sym;
 
 /* Used by NrnProperty */
 static _nrn_mechanism_std_vector<double> _parm_default{
     0.005, /* gkbar */
     0.06, /* gnabar */
 }; 
 
 
extern Prop* need_memb(Symbol*);
static void nrn_alloc(Prop* _prop) {
  Prop *prop_ion{};
  Datum *_ppvar{};
   _ppvar = nrn_prop_datum_alloc(_mechtype, 7, _prop);
    _nrn_mechanism_access_dparam(_prop) = _ppvar;
     _nrn_mechanism_cache_instance _ml_real{_prop};
    auto* const _ml = &_ml_real;
    size_t const _iml{};
    assert(_nrn_mechanism_get_num_vars(_prop) == 22);
 	/*initialize range parameters*/
 	gkbar = _parm_default[0]; /* 0.005 */
 	gnabar = _parm_default[1]; /* 0.06 */
 	 assert(_nrn_mechanism_get_num_vars(_prop) == 22);
 	_nrn_mechanism_access_dparam(_prop) = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0] = _nrn_mechanism_get_param_handle(prop_ion, 0); /* ek */
 	_ppvar[1] = _nrn_mechanism_get_param_handle(prop_ion, 3); /* ik */
 	_ppvar[2] = _nrn_mechanism_get_param_handle(prop_ion, 4); /* _ion_dikdv */
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[3] = _nrn_mechanism_get_param_handle(prop_ion, 0); /* ena */
 	_ppvar[4] = _nrn_mechanism_get_param_handle(prop_ion, 3); /* ina */
 	_ppvar[5] = _nrn_mechanism_get_param_handle(prop_ion, 4); /* _ion_dinadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 {0, 0}
};
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
void _nrn_thread_table_reg(int, nrn_thread_table_check_t);
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 extern "C" void _kna_custom_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("k", -10000.);
 	ion_reg("na", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
 hoc_register_parm_default(_mechtype, &_parm_default);
         hoc_register_npy_direct(_mechtype, npy_direct_func_proc);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  register_nmodl_text_and_filename(_mechtype);
#endif
   _nrn_mechanism_register_data_fields(_mechtype,
                                       _nrn_mechanism_field<double>{"gkbar"} /* 0 */,
                                       _nrn_mechanism_field<double>{"gnabar"} /* 1 */,
                                       _nrn_mechanism_field<double>{"gk"} /* 2 */,
                                       _nrn_mechanism_field<double>{"gna"} /* 3 */,
                                       _nrn_mechanism_field<double>{"n"} /* 4 */,
                                       _nrn_mechanism_field<double>{"m"} /* 5 */,
                                       _nrn_mechanism_field<double>{"h"} /* 6 */,
                                       _nrn_mechanism_field<double>{"Dn"} /* 7 */,
                                       _nrn_mechanism_field<double>{"Dm"} /* 8 */,
                                       _nrn_mechanism_field<double>{"Dh"} /* 9 */,
                                       _nrn_mechanism_field<double>{"ek"} /* 10 */,
                                       _nrn_mechanism_field<double>{"ena"} /* 11 */,
                                       _nrn_mechanism_field<double>{"ik"} /* 12 */,
                                       _nrn_mechanism_field<double>{"ina"} /* 13 */,
                                       _nrn_mechanism_field<double>{"alpha_n"} /* 14 */,
                                       _nrn_mechanism_field<double>{"beta_n"} /* 15 */,
                                       _nrn_mechanism_field<double>{"alpha_m"} /* 16 */,
                                       _nrn_mechanism_field<double>{"beta_m"} /* 17 */,
                                       _nrn_mechanism_field<double>{"alpha_h"} /* 18 */,
                                       _nrn_mechanism_field<double>{"beta_h"} /* 19 */,
                                       _nrn_mechanism_field<double>{"v"} /* 20 */,
                                       _nrn_mechanism_field<double>{"_g"} /* 21 */,
                                       _nrn_mechanism_field<double*>{"_ion_ek", "k_ion"} /* 0 */,
                                       _nrn_mechanism_field<double*>{"_ion_ik", "k_ion"} /* 1 */,
                                       _nrn_mechanism_field<double*>{"_ion_dikdv", "k_ion"} /* 2 */,
                                       _nrn_mechanism_field<double*>{"_ion_ena", "na_ion"} /* 3 */,
                                       _nrn_mechanism_field<double*>{"_ion_ina", "na_ion"} /* 4 */,
                                       _nrn_mechanism_field<double*>{"_ion_dinadv", "na_ion"} /* 5 */,
                                       _nrn_mechanism_field<int>{"_cvode_ieq", "cvodeieq"} /* 6 */);
  hoc_register_prop_size(_mechtype, 22, 7);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 6, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 
    hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 kna_custom /Users/carlos/Desktop/neuron/numerical_stability_in_neurovisor/yale/kna_custom.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static const char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_internalthreadargsprotocomma_ double);
 
static int _ode_spec1(_internalthreadargsproto_);
/*static int _ode_matsol1(_internalthreadargsproto_);*/
 static neuron::container::field_index _slist1[3], _dlist1[3];
 static int states(_internalthreadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (_internalthreadargsproto_) {int _reset = 0; {
   rates ( _threadargscomma_ v ) ;
   Dn = alpha_n * ( 1.0 - n ) - beta_n * n ;
   Dm = alpha_m * ( 1.0 - m ) - beta_m * m ;
   Dh = alpha_h * ( 1.0 - h ) - beta_h * h ;
   }
 return _reset;
}
 static int _ode_matsol1 (_internalthreadargsproto_) {
 rates ( _threadargscomma_ v ) ;
 Dn = Dn  / (1. - dt*( ( alpha_n )*( ( ( - 1.0 ) ) ) - ( beta_n )*( 1.0 ) )) ;
 Dm = Dm  / (1. - dt*( ( alpha_m )*( ( ( - 1.0 ) ) ) - ( beta_m )*( 1.0 ) )) ;
 Dh = Dh  / (1. - dt*( ( alpha_h )*( ( ( - 1.0 ) ) ) - ( beta_h )*( 1.0 ) )) ;
  return 0;
}
 /*END CVODE*/
 static int states (_internalthreadargsproto_) { {
   rates ( _threadargscomma_ v ) ;
    n = n + (1. - exp(dt*(( alpha_n )*( ( ( - 1.0 ) ) ) - ( beta_n )*( 1.0 ))))*(- ( ( alpha_n )*( ( 1.0 ) ) ) / ( ( alpha_n )*( ( ( - 1.0 ) ) ) - ( beta_n )*( 1.0 ) ) - n) ;
    m = m + (1. - exp(dt*(( alpha_m )*( ( ( - 1.0 ) ) ) - ( beta_m )*( 1.0 ))))*(- ( ( alpha_m )*( ( 1.0 ) ) ) / ( ( alpha_m )*( ( ( - 1.0 ) ) ) - ( beta_m )*( 1.0 ) ) - m) ;
    h = h + (1. - exp(dt*(( alpha_h )*( ( ( - 1.0 ) ) ) - ( beta_h )*( 1.0 ))))*(- ( ( alpha_h )*( ( 1.0 ) ) ) / ( ( alpha_h )*( ( ( - 1.0 ) ) ) - ( beta_h )*( 1.0 ) ) - h) ;
   }
  return 0;
}
 
static int  rates ( _internalthreadargsprotocomma_ double _lv ) {
   double _lvm ;
 _lvm = _lv ;
   alpha_n = 1.0 * 0.032 * ( 15.0 - _lvm ) / ( exp ( ( 15.0 - _lvm ) / 5.0 ) - 1.0 ) ;
   beta_n = 1.0 * 0.5 * exp ( ( 10.0 - _lvm ) / 40.0 ) ;
   alpha_m = 1.0 * 0.32 * ( 13.0 - _lvm ) / ( exp ( ( 13.0 - _lvm ) / 4.0 ) - 1.0 ) ;
   beta_m = 1.0 * 0.28 * ( _lvm - 40.0 ) / ( exp ( ( _lvm - 40.0 ) / 5.0 ) - 1.0 ) ;
   alpha_h = 1.0 * 0.128 * exp ( ( 17.0 - _lvm ) / 18.0 ) ;
   beta_h = 1.0 * 4.0 / ( exp ( ( 40.0 - _lvm ) / 5.0 ) + 1.0 ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 
  Prop* _local_prop = _prop_id ? _extcall_prop : nullptr;
  _nrn_mechanism_cache_instance _ml_real{_local_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _local_prop ? _nrn_mechanism_access_dparam(_local_prop) : nullptr;
_thread = _extcall_thread.data();
double* _globals = nullptr;
if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
_nt = nrn_threads;
 _r = 1.;
 rates ( _threadargscomma_ *getarg(1) );
 hoc_retpushx(_r);
}
 
static double _npy_rates(Prop* _prop) {
    double _r{0.0};
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 _nrn_mechanism_cache_instance _ml_real{_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _nrn_mechanism_access_dparam(_prop);
_thread = _extcall_thread.data();
double* _globals = nullptr;
if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
_nt = nrn_threads;
 _r = 1.;
 rates ( _threadargscomma_ *getarg(1) );
 return(_r);
}
 
static int _ode_count(int _type){ return 3;}
 
static void _ode_spec(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
   Datum* _ppvar;
   size_t _iml;   _nrn_mechanism_cache_range* _ml;   Node* _nd{};
  double _v{};
  int _cntml;
  _nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
  _ml = &_lmr;
  _cntml = _ml_arg->_nodecount;
  Datum *_thread{_ml_arg->_thread};
  double* _globals = nullptr;
  if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _ppvar = _ml_arg->_pdata[_iml];
    _nd = _ml_arg->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
  ena = _ion_ena;
     _ode_spec1 (_threadargs_);
   }}
 
static void _ode_map(Prop* _prop, int _ieq, neuron::container::data_handle<double>* _pv, neuron::container::data_handle<double>* _pvdot, double* _atol, int _type) { 
  Datum* _ppvar;
  _ppvar = _nrn_mechanism_access_dparam(_prop);
  _cvode_ieq = _ieq;
  for (int _i=0; _i < 3; ++_i) {
    _pv[_i] = _nrn_mechanism_get_param_handle(_prop, _slist1[_i]);
    _pvdot[_i] = _nrn_mechanism_get_param_handle(_prop, _dlist1[_i]);
    _cvode_abstol(_atollist, _atol, _i);
  }
 }
 
static void _ode_matsol_instance1(_internalthreadargsproto_) {
 _ode_matsol1 (_threadargs_);
 }
 
static void _ode_matsol(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
   Datum* _ppvar;
   size_t _iml;   _nrn_mechanism_cache_range* _ml;   Node* _nd{};
  double _v{};
  int _cntml;
  _nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
  _ml = &_lmr;
  _cntml = _ml_arg->_nodecount;
  Datum *_thread{_ml_arg->_thread};
  double* _globals = nullptr;
  if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _ppvar = _ml_arg->_pdata[_iml];
    _nd = _ml_arg->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
  ena = _ion_ena;
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel(_internalthreadargsproto_) {
  int _i; double _save;{
  h = h0;
  m = m0;
  n = n0;
 
}
}

static void nrn_init(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type){
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto* const _vec_v = _nt->node_voltage_storage();
auto* const _ml = &_lmr;
Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
_thread = _ml_arg->_thread;
double* _globals = nullptr;
if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
for (_iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
   _v = _vec_v[_ni[_iml]];
 v = _v;
  ek = _ion_ek;
  ena = _ion_ena;
 initmodel(_threadargs_);
  }
}

static double _nrn_current(_internalthreadargsprotocomma_ double _v) {
double _current=0.; v=_v;
{ {
   gk = gkbar * pow( n , 4.0 ) ;
   gna = gnabar * pow( m , 3.0 ) * h ;
   ik = gk * ( v - ek ) ;
   ina = gna * ( v - ena ) ;
   }
 _current += ik;
 _current += ina;

} return _current;
}

static void nrn_cur(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto const _vec_rhs = _nt->node_rhs_storage();
auto const _vec_sav_rhs = _nt->node_sav_rhs_storage();
auto const _vec_v = _nt->node_voltage_storage();
auto* const _ml = &_lmr;
Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
_thread = _ml_arg->_thread;
double* _globals = nullptr;
if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
for (_iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
   _v = _vec_v[_ni[_iml]];
  ek = _ion_ek;
  ena = _ion_ena;
 auto const _g_local = _nrn_current(_threadargscomma_ _v + .001);
 	{ double _dina;
 double _dik;
  _dik = ik;
  _dina = ina;
 _rhs = _nrn_current(_threadargscomma_ _v);
  _ion_dikdv += (_dik - ik)/.001 ;
  _ion_dinadv += (_dina - ina)/.001 ;
 	}
 _g = (_g_local - _rhs)/.001;
  _ion_ik += ik ;
  _ion_ina += ina ;
	 _vec_rhs[_ni[_iml]] -= _rhs;
 
}
 
}

static void nrn_jacob(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto const _vec_d = _nt->node_d_storage();
auto const _vec_sav_d = _nt->node_sav_d_storage();
auto* const _ml = &_lmr;
Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
_thread = _ml_arg->_thread;
double* _globals = nullptr;
if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
for (_iml = 0; _iml < _cntml; ++_iml) {
  _vec_d[_ni[_iml]] += _g;
 
}
 
}

static void nrn_state(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto* const _vec_v = _nt->node_voltage_storage();
auto* const _ml = &_lmr;
Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni;
_ni = _ml_arg->_nodeindices;
size_t _cntml = _ml_arg->_nodecount;
_thread = _ml_arg->_thread;
double* _globals = nullptr;
if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
for (size_t _iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
 _nd = _ml_arg->_nodelist[_iml];
   _v = _vec_v[_ni[_iml]];
 v=_v;
{
  ek = _ion_ek;
  ena = _ion_ena;
 {   states(_threadargs_);
  }  }}

}

static void terminal(){}

static void _initlists(){
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = {n_columnindex, 0};  _dlist1[0] = {Dn_columnindex, 0};
 _slist1[1] = {m_columnindex, 0};  _dlist1[1] = {Dm_columnindex, 0};
 _slist1[2] = {h_columnindex, 0};  _dlist1[2] = {Dh_columnindex, 0};
_first = 0;
}

#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mech_type) {
    const char* nmodl_filename = "/Users/carlos/Desktop/neuron/numerical_stability_in_neurovisor/yale/kna_custom.mod";
    const char* nmodl_file_text = 
  "COMMENT\n"
  "  kna_custom.mod\n"
  "  A Hodgkin-Huxley style channel containing:\n"
  "    - Potassium gate n with alpha_n, beta_n from PDE code\n"
  "    - Sodium gates m, h with alpha_m, beta_m, alpha_h, beta_h from PDE code\n"
  "  No leak current included. Reversal potentials are read from NEURON's ek, ena.\n"
  "  Units are consistent with NEURON (mV, ms, mA/cm2, S/cm2).\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "  SUFFIX kna_custom\n"
  "  USEION k  READ ek  WRITE ik\n"
  "  USEION na READ ena WRITE ina\n"
  "  RANGE gkbar, gk, gna, gnabar\n"
  "  RANGE n, m, h\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "  gkbar  = 0.005 (S/cm2)  : K conductance (convert from PDE  50 S/m2 => 0.005 S/cm2)\n"
  "  gnabar = 0.06  (S/cm2)  : Na conductance (convert from PDE 600 S/m2 => 0.06 S/cm2)\n"
  "}\n"
  "\n"
  "STATE {\n"
  "  n\n"
  "  m\n"
  "  h\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "  v   (mV)\n"
  "  ek  (mV)\n"
  "  ena (mV)\n"
  "  ik  (mA/cm2)\n"
  "  ina (mA/cm2)\n"
  "  alpha_n (1/ms)\n"
  "  beta_n  (1/ms)\n"
  "  alpha_m (1/ms)\n"
  "  beta_m  (1/ms)\n"
  "  alpha_h (1/ms)\n"
  "  beta_h  (1/ms)\n"
  "  gk  (S/cm2)\n"
  "  gna (S/cm2)\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "  SOLVE states METHOD cnexp\n"
  "  gk  = gkbar  * n^4\n"
  "  gna = gnabar * m^3 * h\n"
  "  ik  = gk  * (v - ek)\n"
  "  ina = gna * (v - ena)\n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "  rates(v)\n"
  "  n' = alpha_n*(1 - n) - beta_n*n\n"
  "  m' = alpha_m*(1 - m) - beta_m*m\n"
  "  h' = alpha_h*(1 - h) - beta_h*h\n"
  "}\n"
  "\n"
  "PROCEDURE rates(v (mV)) {\n"
  "  LOCAL vm  : v in mV\n"
  "  vm = v\n"
  "\n"
  "  : alpha_n = 1.0e3 * 0.032*(15 - Vin) / ( exp((15 - Vin)/5) - 1 )\n"
  "  alpha_n = 1(e3)*0.032 * (15 - vm)/( exp((15 - vm)/5) - 1 )\n"
  "  beta_n  = 1(e3)*0.5   * exp((10 - vm)/40)\n"
  "\n"
  "  : alpha_m = 1.0e3 * 0.32*(13 - Vin) / ( exp((13 - Vin)/4) - 1 )\n"
  "  alpha_m = 1(e3)*0.32 * (13 - vm)/( exp((13 - vm)/4) - 1 )\n"
  "  beta_m  = 1(e3)*0.28 * (vm - 40)/( exp((vm - 40)/5) - 1 )\n"
  "\n"
  "  : alpha_h = 1.0e3 * 0.128 * exp((17 - Vin)/18)\n"
  "  alpha_h = 1(e3)*0.128 * exp( (17 - vm)/18 )\n"
  "  beta_h  = 1(e3)*4.0   / (exp((40 - vm)/5) + 1 )\n"
  "}\n"
  ;
    hoc_reg_nmodl_filename(mech_type, nmodl_filename);
    hoc_reg_nmodl_text(mech_type, nmodl_file_text);
}
#endif
