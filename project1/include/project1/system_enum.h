#pragma once

#include <string>

enum SystemEnum {
    SIMPLE_ANALYTICAL,
    SIMPLE_ANALYTICAL_IMPORTANCE,
    SIMPLE_NUMERICAL,
    INTERACTIVE,
    NO_JASTROW,
};

std::string system_enum_to_string(SystemEnum system_enum);
