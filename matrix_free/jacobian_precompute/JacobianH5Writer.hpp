#ifndef JACOBIANH5WRITER_HPP
#define JACOBIANH5WRITER_HPP

#include <string>
#include <vector>

class JacobianH5Writer
{
    public:
        void Write(const std::string &filename,
            const int &p,
            const int &q,
            const int &nlocalelemx,
            const int &nlocalelemy,
            const int &nqp1,
            const int &nqp2,
            const std::vector<double> &elem_size1,
            const std::vector<double> &elem_size2,
            const std::vector<double> &qp1,
            const std::vector<double> &qp2,
            const std::vector<double> &w1,
            const std::vector<double> &w2,
            const std::vector<double> &inv_jacobian,
            const std::vector<double> &det_jacobian) const;
};

#endif
