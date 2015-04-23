#ifndef MOLECULE
#define MOLECULE

#include "utils.h"
#include "sysutils.h"

namespace phys {

// map: first - name, second - atom mass
std::map<std::string, int> _mendel;
// map: first - name, second - valence
std::map<std::string, int> _valences;
// map: first - bond, second - num;
std::map<std::string, int> _bonds;
// map: first - bond, second - num;
std::map<std::string, float> _radii;

const std::string resFile_HO = "resFile_HO";
const std::string resFile_HCa = "resFile_HCa";
const std::string resFile_HP = "resFile_HP";
const std::string resFile_OO = "resFile_OO";
const std::string resFile_OP = "resFile_OP";
const std::string resFile_OCa = "resFile_OCa";
const std::string resFile_PP = "resFile_PP";
const std::string resFile_PCa = "resFile_PCa";
const std::string resFile_CaCa = "resFile_CaCa";

void fillMendelAndValences()
{
    _mendel["H"] = 1;
    _mendel["O"] = 16;
    _mendel["P"] = 31;
    _mendel["Ca"] = 40;

    _valences["H"] = 1;
    _valences["O"] = 2;
    _valences["P"] = 5;
    _valences["Ca"] = 2;

    _radii["H"] = 0.31f;
    _radii["O"] = 0.66f;
    _radii["P"] = 1.07f;
    _radii["Ca"] = 1.76f;

    int i = 0;
    _bonds["HO"] = ++i;
    _bonds["HCa"] = ++i;
    _bonds["HP"] = ++i;
    _bonds["OO"] = ++i;
    _bonds["OP"] = ++i;
    _bonds["OCa"] = ++i;
    _bonds["PP"] = ++i;
    _bonds["PCa"] = ++i;
    _bonds["CaCa"] = ++i;

    std::ofstream outfile(resFile_HO.c_str());
    outfile << resFile_HO << std::endl;
    outfile.close();

    std::ofstream outfile1(resFile_HCa.c_str());
    outfile1 << resFile_HCa << std::endl;
    outfile1.close();

//    std::ofstream outfile2(resFile_HO.c_str());
//    outfile2 << resFile_HO << std::endl;
//    outfile2.close();

    std::ofstream outfile3(resFile_HP.c_str());
    outfile3 << resFile_HP << std::endl;
    outfile3.close();

    std::ofstream outfile4(resFile_OO.c_str());
    outfile4 << resFile_OO << std::endl;
    outfile4.close();

    std::ofstream  outfile5(resFile_OP.c_str());
    outfile5 << resFile_OP << std::endl;
    outfile5.close();
    
    std::ofstream outfile6(resFile_CaCa.c_str());
    outfile6 << resFile_CaCa << std::endl;
    outfile6.close();

    std::ofstream outfile7(resFile_PP.c_str());
    outfile << resFile_PP << std::endl;
    outfile.close();

    std::ofstream outfile8(resFile_PCa.c_str());
    outfile << resFile_PCa << std::endl;
    outfile.close();

    std::ofstream outfile9(resFile_CaCa.c_str());
    outfile << resFile_CaCa << std::endl;
    outfile.close();
}

struct Atom
{
    std::string mName;
    Vector3D mCoord;
    int mValence;
};

struct Molecule
{
    std::vector<Atom> mAtoms;
    std::string mName;
    std::string mPath;
    std::map<std::string, int> mBonds;
    unsigned int hash;
    float mEnergy;

    void getHz(std::vector<uint>& distances) const
    {
        distances.clear();
        if (!mAtoms.empty())
        {
            const int size = mAtoms.size();
            for (int i = 0; i < size; ++i)
                for (int j = i+1; j < size; ++j)
                {
                    const int dist = 10000.f * (mAtoms[i].mCoord - mAtoms[j].mCoord).Magnitude2();
                    distances.push_back(dist);
                }
            std::sort(distances.begin(), distances.end());
        }
    }
    
    bool equal(const Molecule& mol) const
    {
        bool res = false;
        std::vector<uint> first;
        std::vector<uint> second;
        getHz(first);
        mol.getHz(second);
        
        if (!first.empty() && !second.empty() &&
            first.size() == second.size())
        {
            res = true;
            for (int i = 0, size = first.size();
                 i < size && res; ++i)
                res = first[i] == second[i];
        }
        return res;
    }
    
    uint getHash()
    {
        hash = 0;
        if (!mBonds.empty() && !mAtoms.empty())
        {
            for (std::map<std::string, int>::iterator it = _bonds.begin();
                 it != _bonds.end(); ++it)
            {
                unsigned int bond = mBonds[it->first];
                hash += bond;
                hash *= 10;
            }
        }
        return hash;
    }

    uint mustBeBonds() const
    {
        uint num = 0;
        for (int i = 0, size = mAtoms.size(); i < size; ++i)
            num += _valences[mAtoms[i].mName];
        return num / 2;
    }

    uint getBondsNum()
    {
        uint num = 0;
        if (!mBonds.empty())
        {
            for (std::map<std::string, int>::iterator it = _bonds.begin();
                 it != _bonds.end(); ++it)
                num += mBonds[it->first];
        }
        return num;
    }

    void fillDistancesAndBonds1()
    {
        mBonds.clear();
        const int size = mAtoms.size();
    
        std::map<int, std::vector<std::pair<int, int> > > distances;
    
        uint numSemiBonds = 0;
        for (int i = 0; i < size; ++i)
        {
            Atom atom1 = mAtoms[i];
            numSemiBonds += atom1.mValence;
            for (int j = i + 1; j < size; ++j)
            {
                Atom atom2 = mAtoms[j];
                int magnitude = int(1000 * Magnitude(atom1.mCoord - atom2.mCoord));
                distances[magnitude].push_back(std::make_pair(i, j));
            }
        }
    
    //    int num = 0;
    //    for (std::map<int, std::vector<std::pair<int, int> > >::iterator it = distances.begin();
    //         num < numSemiBonds / 2; ++it)
    //    {
    //        int sizeSize = it->second.size();
    //        for (int i = 0; i < sizeSize; ++i)
    //        {
    //            std::pair<int, int> thisPair = it->second[i];
    //            const std::string name1 = mol.mAtoms[thisPair.first].mName;
    //            const std::string name2 = mol.mAtoms[thisPair.second].mName;
    //            const std::string res = _mendel[name1] < _mendel[name2] ? name1 + name2 : name2 + name1;
    //            ++mol.mBonds[res];
    //            ++num;
    //        }
    //    }
    
    
        int num = 0;
        for (std::map<int, std::vector<std::pair<int, int> > >::iterator it = distances.begin();
             num < numSemiBonds / 2 && it != distances.end(); ++it)
        {
            std::vector<std::pair<int, int> > tmpBonds = it->second;
            const int tmpBondsSize = tmpBonds.size();
            for (int i = 0; i < tmpBondsSize && num < numSemiBonds / 2; ++i)
            {
                std::pair<int, int> bond = tmpBonds[i];
                const std::string name1 = mAtoms[bond.first].mName;
                const std::string name2 = mAtoms[bond.second].mName;
    
                const float equilibriumDist = _radii[name1] + _radii[name2];
                const float curDist = it->first / 1000.f;
                const float ratioDist = fabs(curDist - equilibriumDist) / equilibriumDist;
    //            if ((name1 == "Ca" || name2 == "Ca") && ratioDist < 0.4f)
    //            {
    //                const std::string res = _mendel[name1] < _mendel[name2] ? name1 + name2 : name2 + name1;
    //                ++mol.mBonds[res];
    //                ++num;
    //            } else
                if (ratioDist < 0.2f)
                {
                    const std::string res = _mendel[name1] < _mendel[name2] ? name1 + name2 : name2 + name1;
                    ++mBonds[res];
                    ++num;
                }
            }
        }
    }
    
    bool readFromFile()
    {
        if (!mPath.empty())
        {
            FILE* inFile = fopen(mPath.c_str(), "r");
            if (inFile)
            {
                int numAtoms = 0;
                fscanf(inFile, "%d", &numAtoms);
                
                mAtoms.assign(numAtoms, Atom());
                
                for (int j = 0; j < numAtoms; ++j)
                {
                    Atom* atom = &(mAtoms[j]);
                    char name[4]; // max name of atom have 4 symbols
                    fscanf(inFile, "%s", &name);
                    atom->mName = std::string(name);
                    fscanf(inFile, "%f %f %f", &(atom->mCoord.x), &(atom->mCoord.y), &(atom->mCoord.z));
                    atom->mValence = _valences[name];
                }
                
                fclose(inFile);
                return true;
            }
        }
        
        return false;
    }
    
    uint findDoubleOBond()
    {
        assert(!mAtoms.empty());
    
        std::ofstream ofs("resultDoubleBonds", std::ios_base::app);
    
        uint numDoubleBonds = 0;
        const uint size = mAtoms.size();
        for (int i = 0; i < size; ++i)
        {
            const Atom* atom1 = &(mAtoms[i]);
            if (0 == atom1->mName.compare("O"))
            {
                int endAtom = -1;
                uint numNeighborAtoms = 0;
                for (int j = 0; j < size; ++j)
                {
                    const Atom* atom2 = &(mAtoms[j]);
    
                    const float equilibriumDist = _radii["O"] + _radii[atom2->mName];
                    const float dist = (atom1->mCoord - atom2->mCoord).Magnitude();
                    const float ratio = fabs(dist - equilibriumDist) / equilibriumDist;
                    if (ratio < 0.2f)
                    {
                        ++numNeighborAtoms;
                        endAtom = j;
                    }
                }
    
                if (numNeighborAtoms == 1)
                {
                    ++numDoubleBonds;
                    ofs << mName << "\t" << i << "\t" << endAtom << "\t" << mEnergy << std::endl;
                    std::cout << mName << "\t" << i << "\t" << endAtom << "\t" << mEnergy << std::endl;
                }
            }
        }
    
        ofs.close();
        return numDoubleBonds;
    }
    
};

}

#endif // MOLECULE

