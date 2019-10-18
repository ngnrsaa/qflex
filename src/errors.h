#ifndef ERROR__H
#define ERROR__H

#include <stdexcept>
#include <string>

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define ERROR_MSG(MSG) "ERROR (" __FILE__ ":" TOSTRING(__LINE__) ") --> " + std::string(MSG)

#endif
