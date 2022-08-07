#if !defined(LINALG_H)
#define LINALG_H

extern int my_dgtsv_b(const int n, double *tdm_l, double *tdm_c, double *tdm_u, double *tdm_q);
extern int my_dgtsv_p(const int n, double *tdm_l, double *tdm_c, double *tdm_u, double *tdm_q);

#endif // LINALG_H
