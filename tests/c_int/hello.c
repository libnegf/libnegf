/* Basic example of libnegf call from C */

#include "libnegf.h"
#include <stdio.h>
#include <string.h>

void padRight(char *string, char pad, int padded_len) {
    int len = (int) strlen(string);
    if (len >= padded_len) {
        return;
    }
    int i;
    for (i = len; i < padded_len + 1; i++) {
      string[i] = pad;
    }
    return;
}


int main()
{
  int handler[NEGF_HSIZE];
  char realmat[] = "HR.dat";
  char imagmat[] = "HI.dat";
  int target_matrix = 0;
  int nrow = 100;

  //strcpy(realmat, "HR.dat");
  //strcpy(imagmat,"HI.dat");
  //padRight(realmat, ' ', NEGF_LC);
  //padRight(imagmat, ' ', NEGF_LC);
  
  printf("Initializing libNEGF \n");
  negf_init_session(handler);
  negf_init(handler);
  negf_read_hs(handler, &realmat[0], &imagmat[0], &target_matrix);
  negf_set_s_id(handler, &nrow);
  negf_read_input(handler);

  printf("Destroying libNEGF \n");
  negf_destruct_libnegf(handler);
  negf_destruct_session(handler);
  printf("Done \n");
}
