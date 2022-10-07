//#ifndef CIO_H
//#define CIO_H
#include <stdio.h>
#include <stdlib.h>
typedef unsigned char T;
typedef long S;
void cwritebin(T *r,S* row, S*col );
//void cprint(T *x,S* row, S * col);
void getsizebin(S* row, S * col, char * fname);
void creadbin(T *r, S* len, S* col, char * fname);
//#endif
