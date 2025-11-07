#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _cal_custom_reg();
extern void _kna_custom_reg();
extern void _low_calcium_reg();
extern void _original_ca_reg();
extern void _original_k_reg();
extern void _original_na_reg();
extern void _slow_potassium_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," cal_custom.mod");
fprintf(stderr," kna_custom.mod");
fprintf(stderr," low_calcium.mod");
fprintf(stderr," original_ca.mod");
fprintf(stderr," original_k.mod");
fprintf(stderr," original_na.mod");
fprintf(stderr," slow_potassium.mod");
fprintf(stderr, "\n");
    }
_cal_custom_reg();
_kna_custom_reg();
_low_calcium_reg();
_original_ca_reg();
_original_k_reg();
_original_na_reg();
_slow_potassium_reg();
}
