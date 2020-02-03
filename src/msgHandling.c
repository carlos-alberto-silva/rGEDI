#include "R.h"
#include "Rinternals.h"


int errorf(const char* msg, ...) {
    va_list pl;
    va_start(pl, msg);
    REvprintf(msg,pl);
    return 0;
}

int msgf(const char* msg, ...) {
    va_list pl;
    va_start(pl, msg);
    Rvprintf(msg,pl);
    return 0;
}
