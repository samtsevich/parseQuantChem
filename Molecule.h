#ifndef MOLECULE
#define MOLECULE

#include "utils.h"
#include "sysutils.h"
#include <sstream>
#include <set>

namespace phys {

const float E_ex = 14.4f;
const float E_eV = 23.0609f;        // eV to kcal/mol
const float a_Bohr = 0.529177;      // Angstrom

const char* energyHeader = "%nproc=2\n# b3lyp/6-31G(d,p) pop=(ReadRadii,MK)\n\n";

// map: first - name, second - atom mass
std::map<std::string, uint> _mendel;
// map: first - name, second - valence
std::map<std::string, int> _valences;
// set of names of bonds;
std::set<std::string> _bonds;
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

    _bonds.insert("H1O");
    _bonds.insert("H1Ca");
    _bonds.insert("H1P");
    _bonds.insert("O1O");
    _bonds.insert("O1P");
    _bonds.insert("O1Ca");
    //_bonds.insert("P1P");
    _bonds.insert("P1Ca");
    _bonds.insert("Ca1Ca");
    _bonds.insert("O2O");
    _bonds.insert("O2P");
    _bonds.insert("O2Ca");
    _bonds.insert("P2P");
    //_bonds.insert("P2Ca");
    _bonds.insert("Ca2Ca");
    _bonds.insert("ionCa");

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
    outfile7 << resFile_PP << std::endl;
    outfile.close();

    std::ofstream outfile8(resFile_PCa.c_str());
    outfile8 << resFile_PCa << std::endl;
    outfile.close();

    std::ofstream outfile9(resFile_CaCa.c_str());
    outfile9 << resFile_CaCa << std::endl;
    outfile.close();
}

struct Atom
{
    std::string mName;
    Vector3D mCoord;
    uint mValence;
    float mCharge;
    uint Z;
    float radius;
    bool valid;

    Atom() : valid(true)
    {}
};

inline float coloumb(const Atom& atom1, const Atom& atom2)
{
    const float dist = (atom1.mCoord - atom2.mCoord).Magnitude();
    return E_eV * E_ex * atom1.mCharge * atom2.mCharge / dist;
}

struct Molecule
{
    int mNumAtoms;
    bool coherent;
    std::vector<Atom> mAtoms;
    std::string mName;
    std::string mPath;
    std::map<std::string, int> mBonds;
    int** bondsGrid;
    const char* chargeheader;
    const char* dipoleheader;

    uint param;
    float mEnergy;
    bool mExcludeCa;
    float mIonEnergy;

    Molecule() : mNumAtoms(0), coherent(true), param(0), mIonEnergy(0.f), mEnergy(0.f),
                 chargeheader("Charge="), dipoleheader("Dipole="), mExcludeCa(false)
    {
    }

    ~Molecule()
    {
        if (mNumAtoms > 0)
        {
            for (int i = 0; i < mNumAtoms; ++i)
                delete bondsGrid[i];
            delete bondsGrid;
        }
    }

    void setNumAtoms(const int& numAtoms)
    {
        if (0 < (mNumAtoms = numAtoms))
        {
            if (mAtoms.empty())
                mAtoms.assign(mNumAtoms, Atom());
            assert((int)mAtoms.size() == numAtoms);
            bondsGrid = new int*[numAtoms];

            for (int i = 0; i < numAtoms; ++i)
            {
                bondsGrid[i] = new int[numAtoms];
                for (int j = 0; j < numAtoms; ++j)
                    bondsGrid[i][j] = 0;
            }
        }
    }

    void setNumAtoms()
    {
        setNumAtoms((int)mAtoms.size());
    }

    bool isNormBond(int i, int j)
    {
        return i >= 0 && i < mNumAtoms && j >= 0 && j < mNumAtoms;
    }

    int bond(int i, int j)
    {
        return isNormBond(i, j) ? bondsGrid[i][j] : -1;
    }

    void getHz(std::vector<uint>& distances) const
    {
        distances.clear();
        if (!mAtoms.empty())
        {
            for (int i = 0; i < mNumAtoms; ++i)
                for (int j = i+1; j < mNumAtoms; ++j)
                {
                    const int dist = 50.f * (mAtoms[i].mCoord - mAtoms[j].mCoord).Magnitude();
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
            first.size() == second.size() &&
            0 != mName.compare(mol.mName))
        {
            res = true;
            for (int i = 0, size = first.size();
                 i < size && res; ++i)
                res = first[i] == second[i];
        }
        return res;
    }

    ull getHash()
    {
        ull hash = 0;
        if (!mBonds.empty() && !mAtoms.empty())
        {
            for (std::set<std::string>::iterator it = _bonds.begin();
                 it != _bonds.end(); ++it)
                if (*it != "ionCa")
                {
                    uint bond = mBonds[*it];
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
            if (mAtoms[i].valid)
                num += mAtoms[i].mValence;
        return num / 2;
    }

    uint getBondsNum()
    {
        uint num = 0;
        if (!mBonds.empty())
        {
            for (std::set<std::string>::iterator it = _bonds.begin();
                 it != _bonds.end(); ++it)
                num += mBonds[*it];
        }
        return num;
    }

    uint findDoubleOBond()
    {
        assert(!mAtoms.empty());
        std::ofstream ofs("resultDoubleBonds", std::ios_base::app);
        uint numDoubleBonds = 0;
        const uint size = mAtoms.size();
        for (uint i = 0; i < size; ++i)
        {
            const Atom* atom1 = &(mAtoms[i]);
            if (0 == atom1->mName.compare("O"))
            {
                int endAtom = -1;
                uint numNeighborAtoms = 0;
                const Atom* atom2 = 0;
                for (uint j = 0; j < size; ++j)
                {
                    atom2 = &(mAtoms[j]);

                    const float equilibriumDist = _radii["O"] + atom2->radius;
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
                    bondsGrid[i][endAtom] = bondsGrid[endAtom][i] = 2;
                    std::string name = _mendel["O"] < atom2->Z ?
                                std::string("O2") + atom2->mName : atom2->mName + std::string("2O");    // here 2-Oxygen or Oxygen-2
                    mBonds[name] += 2;

                    ++numDoubleBonds;
                    ofs << mName << "\t" << i << "\t" << endAtom << "\t" << mEnergy << std::endl;
                    std::cout << mName << "\t" << i << "\t" << endAtom << "\t" << mEnergy << std::endl;
                }
            }
        }

        ofs.close();
        return numDoubleBonds;
    }

    // Fill distances files and bonds
    void fillDistancesAndBonds(std::ofstream& ofs)
    {

        std::map<std::string, int> bonds;       // num bonds in molecule
        const int size = mAtoms.size();

        std::vector< std::vector<char> > matrixBonds;
        std::vector<char> tmp(size, 0);
        for (int i = 0; i < size; ++i)
            matrixBonds.push_back(tmp);


        std::map<std::pair<std::string,int>, int> binding;     // respond to current binding in molecule

        for (int i = 0; i < size; ++i)
        {
            Atom atom1 = mAtoms[i];
            std::string s1 = atom1.mName;
            std::map<float, std::pair<std::string, int> > distances; // to another atoms
            for (int j = i + 1; j < size; ++j)
            {
                Atom atom2 = mAtoms[j];
                std::string s2 = atom2.mName;

                float magnitude = int(1000 * Magnitude(atom1.mCoord - atom2.mCoord)) / 1000.f;
                distances[magnitude] = std::make_pair(s2, j);

            }

            int v = 0,
                valence = _valences[s1];
            for (std::map<float, std::pair<std::string, int> >::iterator it = distances.begin();
                 it != distances.end() && v < valence; ++it)
            {
                std::pair<std::string, int> endAtom = it->second;
                if (!matrixBonds[i][endAtom.second] && !matrixBonds[endAtom.second][i] &&
                    binding[std::make_pair(s1, i)] < _valences[s1] &&
                    binding[endAtom] < _valences[endAtom.first])
                {
                    ++binding[std::make_pair(s1, i)];
                    ++binding[endAtom];
                    std::string bond = _mendel[s1] < _mendel[endAtom.first] ? s1 + endAtom.first : endAtom.first + s1;
                    ++bonds[bond];
                    matrixBonds[i][endAtom.second] = matrixBonds[endAtom.second][i] = 1;
                    ++v;
                }
            }
        }

        for (std::set<std::string>::iterator it = _bonds.begin();
             it != _bonds.end(); ++it)
            ofs << *it << "\t";
        ofs << std::endl;

        for (std::set<std::string>::iterator it = _bonds.begin();
             it != _bonds.end(); ++it)
            ofs << bonds[*it] << "\t";
        ofs << std::endl;
    }

    void fillDistancesAndBonds1()
    {
        mBonds.clear();
        std::map<int, std::vector<std::pair<int, int> > > distances;

        uint numSemiBonds = -2 * findDoubleOBond();
        //uint numSemiBonds = 0;
        for (int i = 0; i < mNumAtoms; ++i)
        {
            Atom* atom1 = &mAtoms[i];
            numSemiBonds += atom1->mValence;
            for (int j = i + 1; j < mNumAtoms; ++j)
            {
                Atom* atom2 = &mAtoms[j]; 
                if (0 == bondsGrid[i][j])
                {
                    int magnitude = int(1000 * (atom1->mCoord - atom2->mCoord).Magnitude());
                    distances[magnitude].push_back(std::make_pair(i, j));
                }
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

        uint num = 0;
        const std::string singleBond("1");
        for (std::map<int, std::vector<std::pair<int, int> > >::iterator it = distances.begin();
             num < numSemiBonds / 2 && it != distances.end(); ++it)
        {
            std::vector<std::pair<int, int> > tmpBonds = it->second;
            const uint tmpBondsSize = tmpBonds.size();
            for (uint i = 0; i < tmpBondsSize && num < numSemiBonds / 2; ++i)
            {
                std::pair<int, int> bond = tmpBonds[i];
                const Atom* atom1 = &mAtoms[bond.first];
                const Atom* atom2 = &mAtoms[bond.second];
                const std::string name1 = atom1->mName;
                const std::string name2 = atom2->mName;

                const float equilibriumDist = atom1->radius + atom2->radius;
                const float curDist = it->first / 1000.f;
                const float ratioDist = fabs(curDist - equilibriumDist) / equilibriumDist;
                if (ratioDist < 0.2f)
                {
                    const std::string name = atom1->Z < atom2->Z ? name1 + singleBond + name2 : name2 + singleBond + name1;
                    ++mBonds[name];
                    bondsGrid[bond.first][bond.second] = bondsGrid[bond.second][bond.first] = 1;
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

                setNumAtoms(numAtoms);

                Atom* atomH = 0;
                for (int j = 0; j < numAtoms; ++j)
                {
                    Atom* atom = &(mAtoms[j]);
                    char name[4]; // max name of atom have 4 symbols
                    fscanf(inFile, "%s", &name);
                    atom->mName = std::string(name);
                    atomH = 0 == atom->mName.compare("H") ? atom : 0;
                    atom->Z = _mendel[atom->mName];
                    atom->radius = _radii[atom->mName];
                    fscanf(inFile, "%f %f %f", &(atom->mCoord.x), &(atom->mCoord.y), &(atom->mCoord.z));
                    atom->mValence = _valences[name];
                }

                fclose(inFile);

                Vector3D sum;
                if (atomH)
                {
                    for (int i = 0; i < mNumAtoms; ++i) 
                        sum += mAtoms[i].mCoord - atomH->mCoord;
                }

                param = uint(sum.Magnitude2() * 100);
                return true;
            }
        }
        return false;
    }

    bool isCoherent()
    {
        if (!coherent)
            return false;
        vertexes.clear();
        vertexes.assign(mNumAtoms, 0);
        dfs(0);

        for (int i = 0; i < mNumAtoms && coherent; ++i)
            coherent = 1 == vertexes[i];

        return coherent;
    }

    void writeGaussianFile()
    {
        assert(!mAtoms.empty() && coherent);

        std::ofstream ofs((mPath + ".gjf").c_str());
        ofs << energyHeader;
        ofs << mName << "\n\n 0 1\n";

        for (int i = 0; i < mNumAtoms; ++i)
            ofs << " " << mAtoms[i].mName << "\t" <<
                   mAtoms[i].mCoord.x << "\t" <<
                   mAtoms[i].mCoord.y << "\t" <<
                   mAtoms[i].mCoord.z << "\n";
        ofs << "\n Ca 1.188\n\n\n";
        ofs.close();
    }

    void readCharges(const std::string& path)
    {
        FILE* inFile = fopen(path.c_str(), "r");
        if (inFile)
        {
            const int maxSymbols = 200;
            char str[maxSymbols];

            while (inFile != NULL && fgets(str, maxSymbols, inFile))
            {
                std::istringstream iss(str);
                std::string charge, dipole;
                float totalChargeValue;
                iss >> charge >> totalChargeValue >> dipole;
                if (0 == charge.compare(chargeheader) &&
                    0 == dipole.compare(dipoleheader))
                {
                    fgets(str, maxSymbols, inFile);
                    for (int i = 0; i < mNumAtoms; ++i)
                    {
                        fgets(str, maxSymbols, inFile);
                        std::istringstream iss1(str);
                        int n = -1;
                        std::string name;
                        float chargeValue;
                        iss1 >> n >> name >> chargeValue;
                        if (0 == name.compare(mAtoms[i].mName.c_str()))
                            mAtoms[i].mCharge = chargeValue;
                    }
                }
            }
        }

        for (int i = 0; i < mNumAtoms; ++i)
            for (int j = i+1; j < mNumAtoms; ++j)
                mIonEnergy += coloumb(mAtoms[i], mAtoms[j]);
    }

    void excludeCa()
    {
        mExcludeCa = true;
        float energy = 0.f;
        const std::string Ca = "Ca";
        for (int i = 0; i < mNumAtoms; ++i)
        {
            Atom* atom1 = &(mAtoms[i]);
            for (int j = i+1; j < mNumAtoms; ++j)
            {
                Atom* atom2 = &(mAtoms[j]);
                if ((0 == Ca.compare(atom1->mName) ||
                     0 == Ca.compare(atom2->mName)) && 
                     bondsGrid[i][j] && bondsGrid[j][i])
                    energy += coloumb(*atom1, *atom2);
            }
        }
        mBonds["ionCa"] = energy;
        //mEnergy -= energy;

//        mBonds["H1Ca"] = 0;
//        mBonds["O1Ca"] = 0;
//        mBonds["P1Ca"] = 0;
//        mBonds["Ca1Ca"] = 0;
//        mBonds["O2Ca"] = 0;
//        mBonds["P2Ca"] = 0;
//        mBonds["Ca2Ca"] = 0;
    }

    void fillChargesStat(std::map<int, std::pair<int, float> >& allCharges)
    {
        for (int i = 0; i < mNumAtoms; ++i)
        {
            Atom* atom = &mAtoms[i];
            if (0 == atom->Z)
                std::cout << mName << std::endl;
            if (allCharges.end() == allCharges.find(atom->Z))
                allCharges[atom->Z] = std::make_pair(1, atom->mCharge);
            else
            {
                std::pair<int, float>* p = &(allCharges[atom->Z]);
                ++p->first;
                p->second += atom->mCharge;
            }
        }
    }

private:
    std::vector<int> vertexes;

    void dfs(int v)
    {
        if (vertexes[v])
            return;
        vertexes[v] = 1;

        const int size = vertexes.size();
        for (int i = 0; i < size; ++i)
        {
            int bond = bondsGrid[v][i];
            if (!vertexes[i] && bond)
                dfs(i);
        }
    }
};
}

#endif // MOLECULE

