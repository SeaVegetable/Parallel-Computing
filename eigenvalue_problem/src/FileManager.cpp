#include "FileManager.hpp"

void FileManager::WritePartition(const std::string &filename, const std::vector<int> &local_to_global,
    const std::vector<double> &CP, const std::vector<int> &ID, const std::vector<int> &IEN,
    const std::vector<double> &NURBSExtraction1, const std::vector<double> &NURBSExtraction2) const
{
    std::ofstream file(filename.c_str());
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    file << "local_to_global" << std::endl;
    std::copy(local_to_global.begin(), local_to_global.end(), std::ostream_iterator<int>(file, " "));
    file << std::endl;

    file << "CP" << std::endl;
    std::copy(CP.begin(), CP.end(), std::ostream_iterator<int>(file, " "));
    file << std::endl;

    file << "ID" << std::endl;
    std::copy(ID.begin(), ID.end(), std::ostream_iterator<int>(file, " "));
    file << std::endl;

    file << "IEN" << std::endl;
    std::copy(IEN.begin(), IEN.end(), std::ostream_iterator<int>(file, " "));
    file << std::endl;

    file << "NURBSExtraction1" << std::endl;
    std::copy(NURBSExtraction1.begin(), NURBSExtraction1.end(), std::ostream_iterator<double>(file, " "));
    file << std::endl;

    file << "NURBSExtraction2" << std::endl;
    std::copy(NURBSExtraction2.begin(), NURBSExtraction2.end(), std::ostream_iterator<double>(file, " "));
    file << std::endl;

    file.close();
}

void FileManager::ReadPartition(const std::string &filename, std::vector<int> &local_to_global,
    std::vector<double> &CP, std::vector<int> &ID, std::vector<int> &IEN,
    std::vector<double> &NURBSExtraction1, std::vector<double> &NURBSExtraction2) const
{
    std::ifstream file(filename.c_str());
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(file, line))
    {
        if (line == "local_to_global")
        {
            std::string local_to_global_str;
            std::getline(file, local_to_global_str);
            std::istringstream local_to_global_ss(local_to_global_str);
            local_to_global.clear();
            int i;
            while (local_to_global_ss >> i)
            {
                local_to_global.push_back(i);
                if (local_to_global_ss.peek() == ' ')
                    local_to_global_ss.ignore();
            }
        }
        else if (line == "CP")
        {
            std::string CP_str;
            std::getline(file, CP_str);
            std::istringstream CP_ss(CP_str);
            CP.clear();
            double d;
            while (CP_ss >> d)
            {
                CP.push_back(d);
                if (CP_ss.peek() == ' ')
                    CP_ss.ignore();
            }
        }
        else if (line == "ID")
        {
            std::string ID_str;
            std::getline(file, ID_str);
            std::istringstream ID_ss(ID_str);
            ID.clear();
            int i;
            while (ID_ss >> i)
            {
                ID.push_back(i);
                if (ID_ss.peek() == ' ')
                    ID_ss.ignore();
            }
        }
        else if (line == "IEN")
        {
            std::string IEN_str;
            std::getline(file, IEN_str);
            std::istringstream IEN_ss(IEN_str);
            IEN.clear();
            int i;
            while (IEN_ss >> i)
            {
                IEN.push_back(i);
                if (IEN_ss.peek() == ' ')
                    IEN_ss.ignore();
            }
        }
        else if (line == "NURBSExtraction1")
        {
            std::string NURBSExtraction1_str;
            std::getline(file, NURBSExtraction1_str);
            std::istringstream NURBSExtraction1_ss(NURBSExtraction1_str);
            NURBSExtraction1.clear();
            double d;
            while (NURBSExtraction1_ss >> d)
            {
                NURBSExtraction1.push_back(d);
                if (NURBSExtraction1_ss.peek() == ' ')
                    NURBSExtraction1_ss.ignore();
            }
        }
        else if (line == "NURBSExtraction2")
        {
            std::string NURBSExtraction2_str;
            std::getline(file, NURBSExtraction2_str);
            std::istringstream NURBSExtraction2_ss(NURBSExtraction2_str);
            NURBSExtraction2.clear();
            double d;
            while (NURBSExtraction2_ss >> d)
            {
                NURBSExtraction2.push_back(d);
                if (NURBSExtraction2_ss.peek() == ' ')
                    NURBSExtraction2_ss.ignore();
            }
        }
    }
}

void FileManager::WritePreprocessInfo(const std::string &filename, const int &p, const int &q, const double &Lx, const double &Ly,
    const int &nElemX, const int &nElemY, const int &part_num_1d, const int &dim, const std::string &base_name) const
{
    std::ofstream file(filename.c_str());
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    file << "p: " << p << std::endl;
    file << "q: " << q << std::endl;
    file << "Lx: " << Lx << std::endl;
    file << "Ly: " << Ly << std::endl;
    file << "nElemX: " << nElemX << std::endl;
    file << "nElemY: " << nElemY << std::endl;
    file << "part_num_1d: " << part_num_1d << std::endl;
    file << "dim: " << dim << std::endl;
    file << "base_name: " << base_name << std::endl;

    file.close();
}

void FileManager::ReadPreprocessInfo(const std::string &filename, int &p, int &q, double &Lx, double &Ly,
    int &nElemX, int &nElemY, int &part_num_1d, int &dim, std::string &base_name) const
{
    std::ifstream file(filename.c_str());
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(file, line))
    {
        if (line.find("p: ") != std::string::npos)
        {
            p = std::stoi(line.substr(3));
        }
        else if (line.find("q: ") != std::string::npos)
        {
            q = std::stoi(line.substr(3));
        }
        else if (line.find("Lx: ") != std::string::npos)
        {
            Lx = std::stod(line.substr(4));
        }
        else if (line.find("Ly: ") != std::string::npos)
        {
            Ly = std::stod(line.substr(4));
        }
        else if (line.find("nElemX: ") != std::string::npos)
        {
            nElemX = std::stoi(line.substr(8));
        }
        else if (line.find("nElemY: ") != std::string::npos)
        {
            nElemY = std::stoi(line.substr(8));
        }
        else if (line.find("part_num_1d: ") != std::string::npos)
        {
            part_num_1d = std::stoi(line.substr(13));
        }
        else if (line.find("dim: ") != std::string::npos)
        {
            dim = std::stoi(line.substr(5));
        }
        else if (line.find("base_name: ") != std::string::npos)
        {
            base_name = line.substr(11);
        }
    }
}

std::string FileManager::GetPartitionFilename(const std::string &base_name, const int &rank) const
{
    return base_name + "_" + std::to_string(rank) + ".txt";
}