#include "power_measurement.h"
#include <stdio.h>

const char* event_names[] = { "rapl:::PACKAGE_ENERGY:PACKAGE%d",
                              "rapl:::DRAM_ENERGY:PACKAGE%d"};

const int N_EVTS = 2;
const int N_SOCK = 2;

void
add_event(int EventSet, int socket)
{
  for (int i = 0; i < N_EVTS; i++) {
    char buf[255];
    int code;
    int retval;
    PAPI_event_info_t info;

    sprintf(buf,  event_names[i], socket);
    retval = PAPI_event_name_to_code(buf, &code);

    retval = PAPI_get_event_info(code, &info);
    retval = PAPI_add_event(EventSet, code);
    /* if (retval != PAPI_OK) { */
    /*   DEBUGP("Activating multiplex\n"); */
    /*   retval = PAPI_set_multiplex(EventSet); */
    /*   PAPIDIE(retval); */
    /*   retval = PAPI_add_named_event(EventSet, buf); */
    /*   PAPIDIE(retval); */
    /* } */

  }
}
