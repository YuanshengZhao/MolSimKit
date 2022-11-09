#include "dump.h"
#include "global.h"
#include <stdio.h>

void dumpXYZ(const char *fname, const char *elm[])
{
    FILE *fp=fopen(fname,"w");
    fprintf(fp,"%d\nx y z w 0 %.17e\n",natom,bl);
    for(int i=0;i<natom;++i)
        fprintf(fp,"%s %.17e %.17e %.17e\n",elm[typ[i]],x[i][0],x[i][1],x[i][2]);
    fclose(fp);
}
