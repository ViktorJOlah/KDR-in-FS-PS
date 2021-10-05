#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _hh_bj_reg();
extern void _kdr_reg();
extern void _na8st_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," hh_bj.mod");
fprintf(stderr," kdr.mod");
fprintf(stderr," na8st.mod");
fprintf(stderr, "\n");
    }
_hh_bj_reg();
_kdr_reg();
_na8st_reg();
}
