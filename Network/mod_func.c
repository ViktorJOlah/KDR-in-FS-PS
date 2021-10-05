#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _ca_reg();
extern void _gap_reg();
extern void _hh_bj_reg();
extern void _kdr_reg();
extern void _km_reg();
extern void _nethhwbm_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," ca.mod");
fprintf(stderr," gap.mod");
fprintf(stderr," hh_bj.mod");
fprintf(stderr," kdr.mod");
fprintf(stderr," km.mod");
fprintf(stderr," nethhwbm.mod");
fprintf(stderr, "\n");
    }
_ca_reg();
_gap_reg();
_hh_bj_reg();
_kdr_reg();
_km_reg();
_nethhwbm_reg();
}
