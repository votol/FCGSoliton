#include <iostream>
#include "yaml-cpp/yaml.h"
#include "schema.h"
#include "NetCdfWriter.h"
#include "CLmanager.h"
#include "PolynomialOperator.h"
#include "ParameterDefines.h"
#include "ParametersInit.h"

#define alph(n) n
#define norm(n, m) static_cast<unsigned int>(Nfibs) + n * static_cast<unsigned int>(Nfibs) + m
#define anom(n, m) static_cast<unsigned int>(Nfibs) + \
    static_cast<unsigned int>(Nfibs) * static_cast<unsigned int>(Nfibs) + \
    n * static_cast<unsigned int>(Nfibs) + m

using namespace clde;

void buildPolynomial(Polynomial& poly, const std::vector<double>& parameters)
{
    //alpha
    for(unsigned int ind = 0 ; ind < Nfibs; ++ind)
    {
        poly.push_back(Monomial());
        poly.back().coe = -gamma;
        poly.back().inInds.push_back(alph(ind));
        poly.back().outInd = alph(ind);

        poly.push_back(Monomial());
        poly.back().coe = 1.0;
        poly.back().inInds.push_back(alph(ind));
        poly.back().inInds.push_back(alph(ind));
        poly.back().inInds.push_back(alph(ind));
        poly.back().outInd = alph(ind);

        poly.push_back(Monomial());
        poly.back().coe = 2.0;
        poly.back().inInds.push_back(alph(ind));
        poly.back().inInds.push_back(norm(ind,ind));
        poly.back().outInd = alph(ind);

        poly.push_back(Monomial());
        poly.back().coe = 1.0;
        poly.back().inInds.push_back(alph(ind));
        poly.back().inInds.push_back(anom(ind,ind));
        poly.back().outInd = alph(ind);

        if(ind != 0)
        {
            poly.push_back(Monomial());
            poly.back().coe = 1.0;
            poly.back().inInds.push_back(alph(ind - 1));
            poly.back().outInd = alph(ind);
        }

        if(ind != static_cast<unsigned int>(Nfibs) - 1 )
        {
            poly.push_back(Monomial());
            poly.back().coe = 1.0;
            poly.back().inInds.push_back(alph(ind + 1));
            poly.back().outInd = alph(ind);
        }
    }
    //normal
    for(unsigned int ind1 = 0 ; ind1 < Nfibs; ++ind1)
    {
        for(unsigned int ind2 = 0 ; ind2 < Nfibs; ++ind2)
        {
            if(ind1 != 0 )
            {
                poly.push_back(Monomial());
                poly.back().coe = -1.0;
                poly.back().inInds.push_back(norm(ind1 - 1, ind2));
                poly.back().outInd = norm(ind1, ind2);
            }

            if(ind1 != static_cast<unsigned int>(Nfibs) - 1 )
            {
                poly.push_back(Monomial());
                poly.back().coe = -1.0;
                poly.back().inInds.push_back(norm(ind1 + 1, ind2));
                poly.back().outInd = norm(ind1, ind2);
            }
            if(ind2 != 0 )
            {
                poly.push_back(Monomial());
                poly.back().coe = 1.0;
                poly.back().inInds.push_back(norm(ind1, ind2 - 1));
                poly.back().outInd = norm(ind1, ind2);
            }

            if(ind2 != static_cast<unsigned int>(Nfibs) - 1 )
            {
                poly.push_back(Monomial());
                poly.back().coe = 1.0;
                poly.back().inInds.push_back(norm(ind1, ind2 + 1));
                poly.back().outInd = norm(ind1, ind2);
            }

            poly.push_back(Monomial());
            poly.back().coe = 2.0;
            poly.back().inInds.push_back(alph(ind2));
            poly.back().inInds.push_back(alph(ind2));
            poly.back().inInds.push_back(norm(ind1, ind2));
            poly.back().outInd = norm(ind1, ind2);

            poly.push_back(Monomial());
            poly.back().coe = -2.0;
            poly.back().inInds.push_back(alph(ind1));
            poly.back().inInds.push_back(alph(ind1));
            poly.back().inInds.push_back(norm(ind1, ind2));
            poly.back().outInd = norm(ind1, ind2);

            poly.push_back(Monomial());
            poly.back().coe = 1.0;
            poly.back().inInds.push_back(alph(ind2));
            poly.back().inInds.push_back(alph(ind2));
            poly.back().inInds.push_back(anom(ind1, ind2));
            poly.back().outInd = norm(ind1, ind2);

            poly.push_back(Monomial());
            poly.back().coe = -1.0;
            poly.back().inInds.push_back(alph(ind1));
            poly.back().inInds.push_back(alph(ind1));
            poly.back().inInds.push_back(anom(ind1, ind2));
            poly.back().outInd = norm(ind1, ind2);

            poly.push_back(Monomial());
            poly.back().coe = 2.0;
            poly.back().inInds.push_back(norm(ind2, ind2));
            poly.back().inInds.push_back(norm(ind1, ind2));
            poly.back().outInd = norm(ind1, ind2);

            poly.push_back(Monomial());
            poly.back().coe = -2.0;
            poly.back().inInds.push_back(norm(ind1, ind1));
            poly.back().inInds.push_back(norm(ind1, ind2));
            poly.back().outInd = norm(ind1, ind2);

            poly.push_back(Monomial());
            poly.back().coe = 1.0;
            poly.back().inInds.push_back(anom(ind2, ind2));
            poly.back().inInds.push_back(anom(ind1, ind2));
            poly.back().outInd = norm(ind1, ind2);

            poly.push_back(Monomial());
            poly.back().coe = -1.0;
            poly.back().inInds.push_back(anom(ind1, ind1));
            poly.back().inInds.push_back(anom(ind1, ind2));
            poly.back().outInd = norm(ind1, ind2);
        }
    }
    //anomal
    for(unsigned int ind1 = 0 ; ind1 < Nfibs; ++ind1)
    {
        for(unsigned int ind2 = 0 ; ind2 < Nfibs; ++ind2)
        {
            if(ind1 != 0 )
            {
                poly.push_back(Monomial());
                poly.back().coe = 1.0;
                poly.back().inInds.push_back(anom(ind1 - 1, ind2));
                poly.back().outInd = anom(ind1, ind2);
            }

            if(ind1 != static_cast<unsigned int>(Nfibs) - 1 )
            {
                poly.push_back(Monomial());
                poly.back().coe = 1.0;
                poly.back().inInds.push_back(anom(ind1 + 1, ind2));
                poly.back().outInd = anom(ind1, ind2);
            }
            if(ind2 != 0 )
            {
                poly.push_back(Monomial());
                poly.back().coe = 1.0;
                poly.back().inInds.push_back(anom(ind1, ind2 - 1));
                poly.back().outInd = anom(ind1, ind2);
            }

            if(ind2 != static_cast<unsigned int>(Nfibs) - 1 )
            {
                poly.push_back(Monomial());
                poly.back().coe = 1.0;
                poly.back().inInds.push_back(anom(ind1, ind2 + 1));
                poly.back().outInd = anom(ind1, ind2);
            }

            poly.push_back(Monomial());
            poly.back().coe = -2.0 * gamma;
            poly.back().inInds.push_back(anom(ind1, ind2));
            poly.back().outInd = anom(ind1, ind2);

            poly.push_back(Monomial());
            poly.back().coe = 2.0;
            poly.back().inInds.push_back(alph(ind2));
            poly.back().inInds.push_back(alph(ind2));
            poly.back().inInds.push_back(anom(ind1, ind2));
            poly.back().outInd = anom(ind1, ind2);

            poly.push_back(Monomial());
            poly.back().coe = 2.0;
            poly.back().inInds.push_back(alph(ind1));
            poly.back().inInds.push_back(alph(ind1));
            poly.back().inInds.push_back(anom(ind1, ind2));
            poly.back().outInd = anom(ind1, ind2);

            poly.push_back(Monomial());
            poly.back().coe = 1.0;
            poly.back().inInds.push_back(alph(ind2));
            poly.back().inInds.push_back(alph(ind2));
            poly.back().inInds.push_back(norm(ind1, ind2));
            poly.back().outInd = anom(ind1, ind2);

            poly.push_back(Monomial());
            poly.back().coe = 1.0;
            poly.back().inInds.push_back(alph(ind1));
            poly.back().inInds.push_back(alph(ind1));
            poly.back().inInds.push_back(norm(ind1, ind2));
            poly.back().outInd = anom(ind1, ind2);

            poly.push_back(Monomial());
            poly.back().coe = 2.0;
            poly.back().inInds.push_back(norm(ind2, ind2));
            poly.back().inInds.push_back(anom(ind1, ind2));
            poly.back().outInd = anom(ind1, ind2);

            poly.push_back(Monomial());
            poly.back().coe = 2.0;
            poly.back().inInds.push_back(norm(ind1, ind1));
            poly.back().inInds.push_back(anom(ind1, ind2));
            poly.back().outInd = anom(ind1, ind2);

            poly.push_back(Monomial());
            poly.back().coe = 1.0;
            poly.back().inInds.push_back(anom(ind2, ind2));
            poly.back().inInds.push_back(norm(ind1, ind2));
            poly.back().outInd = anom(ind1, ind2);

            poly.push_back(Monomial());
            poly.back().coe = 1.0;
            poly.back().inInds.push_back(anom(ind1, ind1));
            poly.back().inInds.push_back(norm(ind1, ind2));
            poly.back().outInd = anom(ind1, ind2);

            if(ind1 == ind2)
            {
                poly.push_back(Monomial());
                poly.back().coe = 0.5 * L;
                poly.back().inInds.push_back(alph(ind1));
                poly.back().inInds.push_back(alph(ind1));
                poly.back().outInd = anom(ind1, ind2);

                poly.push_back(Monomial());
                poly.back().coe = 0.5 * L;
                poly.back().inInds.push_back(alph(ind2));
                poly.back().inInds.push_back(alph(ind2));
                poly.back().outInd = anom(ind1, ind2);

                poly.push_back(Monomial());
                poly.back().coe = 0.5 * L;
                poly.back().inInds.push_back(anom(ind1, ind1));
                poly.back().outInd = anom(ind1, ind2);

                poly.push_back(Monomial());
                poly.back().coe = 0.5 * L;
                poly.back().inInds.push_back(anom(ind2, ind2));
                poly.back().outInd = anom(ind1, ind2);
            }
        }
    }
}

int main(int argc, char **argv)
{
    YAML::Node config = YAML::LoadFile(argv[1]);
    std::string output_dir = config["properties"]
								  [FCGSolitonSchema::PROPERTY_output_path].as<std::string>();

    ParametersHolder parameters_holder_instance(config["parameters"]);
    std::shared_ptr<ICLmanager> manag = std::make_shared<CLmanager>(config["properties"]);

    Polynomial poly;
    buildPolynomial(poly, parameters_holder_instance.GetParameters());
    OperatorDimension operDim = PolynomialOperator::calculateDimension(poly.begin(), poly.end(), true);
    auto grad = polynomialGradient(poly, operDim.in_dim);
    for(auto it = poly.begin(); it != poly.end(); ++it)
        std::cout<< *it << std::endl;


    return 0;
}
