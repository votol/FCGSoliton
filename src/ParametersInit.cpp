#include "ParametersInit.h"
#include "schema.h"
#include "ParameterDefines.h"
ParametersHolder::ParametersHolder(YAML::Node&& in)
{
    parameters.resize(5);
    Nfibs = in[FCGSolitonSchema::PARAMETER_Nfibs].as<double>();
    Nfibs_calc = in[FCGSolitonSchema::PARAMETER_Nfibs_calc].as<double>();
    L = in[FCGSolitonSchema::PARAMETER_L].as<double>();
    gamma = in[FCGSolitonSchema::PARAMETER_gamma].as<double>();
}
