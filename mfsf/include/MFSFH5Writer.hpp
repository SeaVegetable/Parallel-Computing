#ifndef MFSF_H5_WRITER_HPP
#define MFSF_H5_WRITER_HPP

#include <string>

#include "MFSFPrecompute.hpp"

class MFSFH5Writer
{
    public:
        void Write(const std::string &filename,
            const MFSFPrecomputeData<2> &data) const;

        void Write(const std::string &filename,
            const MFSFPrecomputeData<3> &data) const;
};

#endif
