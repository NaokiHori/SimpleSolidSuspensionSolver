#if !defined(ARRAYS_SUSPENSIONS_H)
#define ARRAYS_SUSPENSIONS_H
/* This file is automatically generated by define_array.py */

/* dux (double, 2D) [0:itot+1] x [0:jsize+1] */
#define DUX_MIN_I (0)
#define DUX_MAX_I (itot+1)
#define DUX_LEN_I (DUX_MAX_I-DUX_MIN_I+1)
#define DUX_MIN_J (0)
#define DUX_MAX_J (jsize+1)
#define DUX_LEN_J (DUX_MAX_J-DUX_MIN_J+1)
#define DUX_MEMSIZE (sizeof(double)*(DUX_LEN_I)*(DUX_LEN_J))
#define DUX(I, J) (dux[(J)*(itot+2)+(I)])

/* duy (double, 2D) [0:itot+1] x [0:jsize+1] */
#define DUY_MIN_I (0)
#define DUY_MAX_I (itot+1)
#define DUY_LEN_I (DUY_MAX_I-DUY_MIN_I+1)
#define DUY_MIN_J (0)
#define DUY_MAX_J (jsize+1)
#define DUY_LEN_J (DUY_MAX_J-DUY_MIN_J+1)
#define DUY_MEMSIZE (sizeof(double)*(DUY_LEN_I)*(DUY_LEN_J))
#define DUY(I, J) (duy[(J)*(itot+2)+(I)])

#endif // ARRAYS_SUSPENSIONS_H
