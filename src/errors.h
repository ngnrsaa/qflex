#ifndef ERROR__H
#define ERROR__H

#include "global.h"
#include "utils.h"

#define ERROR_MSG(...)                                                       \
  qflex::utils::concat("ERROR (" __FILE__ ": ", __FUNCTION__, ":", __LINE__, \
                       ") --> ", __VA_ARGS__)

#endif
