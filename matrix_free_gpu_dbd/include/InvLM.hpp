#ifndef INVLM_HPP
#define INVLM_HPP

#include <vector>

class InvLM
{
public:
    InvLM(const int &nLocBas,
    const int &nLocElem,
    const int &nlocfunc,
    const std::vector<int> ID,
    const std::vector<int> IEN);

private:
    std::vector<int> elemNum{};
    std::vector<int> elemIdx{};
    std::vector<int> baseIdx{};
};

#endif