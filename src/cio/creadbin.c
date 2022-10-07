#include "cio.h"
void creadbin(T *r,S * row, S * col, char * fname)
{
    FILE * f = fopen(fname, "rb");
    char str[13];
    fread(str, sizeof(str), 1, f);
    fread(str, sizeof(str), 1, f);
    // assuming r is in the correct shape. 
    // if not, it needs to be done before this.
    // i.e., call getsizebin and then allocate
    // your array.
    fread(r, sizeof(T), *row * *col, f);
    fclose(f);
}
