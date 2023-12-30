/* This file was automatically generated by CasADi.
   The CasADi copyright holders make no ownership claim of its contents. */
#ifdef __cplusplus
extern "C" {
#endif

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int long long int
#endif

int EKF_Func(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem);
int EKF_Func_alloc_mem(void);
int EKF_Func_init_mem(int mem);
void EKF_Func_free_mem(int mem);
int EKF_Func_checkout(void);
void EKF_Func_release(int mem);
void EKF_Func_incref(void);
void EKF_Func_decref(void);
casadi_int EKF_Func_n_out(void);
casadi_int EKF_Func_n_in(void);
casadi_real EKF_Func_default_in(casadi_int i);
const char* EKF_Func_name_in(casadi_int i);
const char* EKF_Func_name_out(casadi_int i);
const casadi_int* EKF_Func_sparsity_in(casadi_int i);
const casadi_int* EKF_Func_sparsity_out(casadi_int i);
int EKF_Func_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w);
#ifdef __cplusplus
} /* extern "C" */
#endif
