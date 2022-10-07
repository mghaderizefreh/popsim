#include "cio.h"

void getsizebin(S* row, S* col, char * fname){
    FILE * f = fopen(fname, "rb");
    char str[13];
    fread(str, sizeof(str), 1, f);
    *row = atoi(str);
    fread(str, sizeof(str), 1, f);
    *col = atoi(str);
    fclose(f);
}
