#include <fstream>
#include <iostream>
#include <string>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <unistd.h>
#include <netcdf>
#include <algorithm>
#include "yaml-cpp/yaml.h"
#include "schema.h"
#include "NetCdfWriter.h"
#include "CLmanager.h"
#include "PolynomialOperator.h"
#include "ParameterDefines.h"
#include "ParametersInit.h"


#define alph(n) (n)
#define norm(n, m) static_cast<unsigned int>(Nfibs) + (n) * static_cast<unsigned int>(Nfibs) + (m)
#define anom(n, m) static_cast<unsigned int>(Nfibs) + \
    static_cast<unsigned int>(Nfibs) * static_cast<unsigned int>(Nfibs) + \
    (n) * static_cast<unsigned int>(Nfibs) + (m)

using namespace clde;
using namespace netCDF;
using namespace netCDF::exceptions;

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

void loadVector(std::vector<double>& vec, std::string& path)
{
    std::ifstream ifs;

    ifs.open (path + "/init.bin", std::ifstream::in | std::ifstream::binary);

    ifs.read((char *)&(vec.data()[1]), 8 * (vec.size() - 1));

    vec[0] = 1.0;
    ifs.close();
}

void printPoly(const Polynomial& in, const unsigned int& dim)
{
    std::vector<Polynomial> tmp(dim);
    for(auto it = in.begin(); it != in.end(); ++it)
    {
        tmp[it->outInd].push_back(*it);
    }
    std::cout << "------------------------" << std::endl;
    for(auto it= tmp.begin(); it != tmp.end(); ++it)
    {
        for( auto it1 = it->begin(); it1 != it->end(); ++it1)
        {
            std::cout<< *it1;
        }
        std::cout<<std::endl;
    }
    std::cout << "------------------------" << std::endl;
}

class Processor
{
    size_t m_dim;
    CLDataStorage<double> m_input_dev;
    CLDataStorage<double> m_output_dev;
    std::unique_ptr<PolynomialOperator> m_base_oper;
    std::list<PolynomialOperator> m_grad_oper;
    std::vector<double> m_base_result;
    std::vector<double> m_grad_result;
public:
    Processor(Polynomial& poly, OperatorDimension dim, std::shared_ptr<ICLmanager>& m):m_dim(dim.in_dim),
        m_input_dev(dim.in_dim + 1,m), m_output_dev(dim.in_dim + 1, m)
    {
        auto grad = polynomialGradient(poly, dim.in_dim);
        m_base_oper = std::unique_ptr<PolynomialOperator>(new
                           PolynomialOperator(poly.begin(), poly.end(), 1, dim, m));
        for(auto it = grad.begin(); it != grad.end(); ++it)
        {
            m_grad_oper.emplace_back(it->begin(),it->end(), 1,dim, m);
        }
        m_base_result.resize(dim.in_dim);
        m_grad_result.resize(dim.in_dim * dim.in_dim);
    }

    void process(const std::vector<double>& in)
    {
        m_input_dev = in;
        m_base_oper->apply(m_input_dev, m_output_dev, std::vector<double>());

        auto tmp = m_output_dev.read();

        std::copy(++tmp.begin(), tmp.end(), m_base_result.begin());
        auto ins_it = m_grad_result.begin();
        for(auto it = m_grad_oper.begin(); it != m_grad_oper.end(); ++it)
        {
            it->apply(m_input_dev, m_output_dev, std::vector<double>());
            tmp = m_output_dev.read();
            std::copy(++tmp.begin(), tmp.end(), ins_it);
            ins_it += m_dim;
        }
    }
    void save(const std::string& file_name)
    {
        NcFile NcFile_instatnce(file_name, NcFile::replace);
        NcDim dim0 = NcFile_instatnce.addDim("dim0", m_dim);
        NcDim dim1 = NcFile_instatnce.addDim("dim1", m_dim);
        NcVar  var = NcFile_instatnce.addVar("vector", ncDouble, dim0);
        var.putVar(m_base_result.data());
        std::vector<NcDim> NcDims;
        NcDims.push_back(dim0);
        NcDims.push_back(dim1);
        var = NcFile_instatnce.addVar("matrix", ncDouble, NcDims);
        var.putVar(m_grad_result.data());
    }
};

int main(int argc, char **argv)
{
    YAML::Node config = YAML::LoadFile(argv[1]);
    std::string output_dir = config["properties"]
								  [FCGSolitonSchema::PROPERTY_output_path].as<std::string>();
    std::string tmp_dir = config["properties"]
                         [FCGSolitonSchema::PROPERTY_tmp_path].as<std::string>();

    ParametersHolder parameters_holder_instance(config["parameters"]);
    auto parameters = parameters_holder_instance.GetParameters();
    std::shared_ptr<ICLmanager> manag = std::make_shared<CLmanager>(config["properties"]);

    Polynomial poly;
    buildPolynomial(poly, parameters_holder_instance.GetParameters());
    OperatorDimension operDim = PolynomialOperator::calculateDimension(poly.begin(), poly.end(), true);
    //printPoly(poly, operDim.in_dim);
    Processor proc(poly, operDim, manag);
    std::vector<double> input(operDim.in_dim + 1);


    std::string output_file = tmp_dir + "/tmp.nc";
    std::string socket_path = tmp_dir + "/socket";
    std::string exit_command = "exit";
    std::string calc_command = "calc";
    char buf[100];
    int fd,cl,rc;
    if ( (fd = socket(AF_UNIX, SOCK_STREAM, 0)) == -1) {
        perror("socket error");
        exit(-1);
      }

    struct sockaddr_un addr;
    memset(&addr, 0, sizeof(addr));
    addr.sun_family = AF_UNIX;
    strncpy(addr.sun_path, socket_path.c_str(), sizeof(addr.sun_path)-1);

    if (bind(fd, (struct sockaddr*)&addr, sizeof(addr)) == -1) {
        perror("bind error");
        exit(-1);
    }

    if (listen(fd, 5) == -1) {
        perror("listen error");
        exit(-1);
    }

    while (1) {
        if ( (cl = accept(fd, NULL, NULL)) == -1) {
            perror("accept error");
            continue;
        }
        break;
    }
    while ( (rc=read(cl,buf,sizeof(buf))) > 0) {
        loadVector(input, tmp_dir);
        proc.process(input);
        proc.save(output_file);
        write(cl, exit_command.c_str(), 1);
    }
    return 0;
}
