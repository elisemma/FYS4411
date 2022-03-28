#include "project1/system_enum.h"

std::string system_enum_to_string(SystemEnum system_enum) {
    switch (system_enum) {
        case SIMPLE_ANALYTICAL:
            return "simple_analytical_brute_force";
        case SIMPLE_ANALYTICAL_IMPORTANCE:
            return "simple_analytical_importance";
        case SIMPLE_NUMERICAL:
            return "simple_numerical_brute_force";
        case INTERACTIVE:
            return "interactive";
        case NO_JASTROW:
            return "no_jastrow";
        default:
            return "unkown";
    }
}
