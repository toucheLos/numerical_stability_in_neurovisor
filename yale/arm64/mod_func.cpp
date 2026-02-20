#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern "C" void _cal_custom_reg(void);
extern "C" void _kna_custom_reg(void);
extern "C" void _low_calcium_reg(void);
extern "C" void _original_ca_reg(void);
extern "C" void _original_k_reg(void);
extern "C" void _original_na_reg(void);
extern "C" void _slow_potassium_reg(void);

extern "C" void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"cal_custom.mod\"");
    fprintf(stderr, " \"kna_custom.mod\"");
    fprintf(stderr, " \"low_calcium.mod\"");
    fprintf(stderr, " \"original_ca.mod\"");
    fprintf(stderr, " \"original_k.mod\"");
    fprintf(stderr, " \"original_na.mod\"");
    fprintf(stderr, " \"slow_potassium.mod\"");
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
