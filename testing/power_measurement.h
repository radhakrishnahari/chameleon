#include <papi.h>

#ifndef _power_measurement_h_
#define _power_measurement_h_

extern const char* event_names[];

extern const int N_EVTS;
extern const int N_SOCK;

void
add_event(int EventSet, int socket);

#endif /* _power_measurement_h_ */
