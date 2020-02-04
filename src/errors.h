#ifndef ERROR__H
#define ERROR__H

#include "global.h"
#include "utils.h"

#define ERROR_MSG(...)                                                       \
  qflex::utils::concat("ERROR (" __FILE__ ": ", __FUNCTION__, ":", __LINE__, \
                       ") --> ", __VA_ARGS__)

#define WARN_MSG(...) \
  qflex::utils::concat("[", qflex::utils::get_date_time(), "] ", __VA_ARGS__)

#endif
