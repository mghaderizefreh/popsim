#include "cio.h"
void cwritebin(T *r, S* row, S* col )
{
    FILE * f = fopen("genotype.bin", "wb");
    char str[13];
    sprintf(str, "%ld\n",(long) *row);
    fwrite(str, sizeof(str), 1, f);
    sprintf(str, "%ld\n",(long) *col);
    fwrite(str, sizeof(str), 1, f);
    fwrite(r, sizeof(T), *row * *col, f);
    fclose(f);
}
