// tst_mat.c - test MAT file save & load

#include <stdio.h>
#include "sigpro.h"

void
print_var(char *fnam, VAR *vl)
{
    float *f;
    int i;

    printf("%s:", fnam);
    for (i = 0; i < SP_MAXVAR; i++) {
	if (vl[i].name && vl[i].data) {
	    f = (float *) vl[i].data;
	    printf(" %s=%.1f", vl[i].name, *f);
	}
	if (vl[i].last) {
	    break;
	}
    }
    printf("\n");
}

int
main()
{
    
    float a = 1.1F, b = 2.2F, c = 3.3F, d = 4.4F;
    VAR *vl1, *vl2;
    static char *fnam = "tst_mat.mat";
    static char *frmt = "f4";

    vl1 = sp_var_alloc(4);
    sp_var_set(vl1 + 0, "a", &a, 1, 1, frmt);
    sp_var_set(vl1 + 1, "b", &b, 1, 1, frmt);
    sp_var_set(vl1 + 2, "c", &c, 1, 1, frmt);
    sp_var_set(vl1 + 3, "d", &d, 1, 1, frmt);
    sp_mat_save(fnam, vl1);

    vl2 = sp_mat_load(fnam);
    print_var(fnam, vl2);

    sp_var_clear_all();

    return (0);
}

