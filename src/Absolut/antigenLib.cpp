#include "antigenLib.h"
#include "../Ymir/proteins.h"
#include "toml.hpp"
#include <fstream>
#include <sstream>
using namespace std;

static bool loaded = true;
static vector<string> listLoaded = vector<string>();

vector<string> listIDs()
{
    if (listLoaded.size() == 0)
    {
        loaded = false;
        std::pair<superProtein *, vector<int>> whatever = getAntigen("NotToBeFoundID");
        loaded = true;
    }
    return listLoaded;
}

void showListIDs()
{
    cout << printVector(listIDs()) << endl;
}

int recstr2int(const char *str, int h = 0)
{
    return !str[h] ? 5381 : (recstr2int(str, h + 1) * 33) ^ str[h];
}

int str2int(const char *str) // actually, unsigned int
{
    if (!loaded)
    {
        if (string(str).compare(string("NotToBeFoundID")))
        {
            listLoaded.push_back(string(str));
        }
    }
    return recstr2int(str);
}

string IDshortcut(int ID)
{
    vector<string> list = listIDs();

    if ((ID < 1) || (ID > static_cast<int>(list.size())))
    {
        cerr << "NO Antigen known with integer ID " << ID << endl;
        cerr << "For your information, antigens are: \n";
        showListIDs();
        return string("");
    }
    return list[ID];
}

string IDshortcut(string justPDB)
{
    vector<string> listAG = listIDs();
    string res = "";
    for (size_t i = 0; i < listAG.size(); ++i)
    {
        std::size_t found = listAG[i].find(justPDB);
        if (found != std::string::npos)
        {
            if (res.size() > 0)
                res = res + " ";
            res = res + listAG[i];
        }
    }
    return res;
}

void testIDshortcut()
{
    cout << IDshortcut("1NCA_N");
    cout << IDshortcut("1NCA");
    cout << IDshortcut("1CZ8");
    cout << IDshortcut("1CZ8_VW_S");
}

static vector<int> lastSearchGlycan;
static vector<vector<int>> lastSearchHotspotsCore;
static vector<vector<int>> lastSearchHotspotsLarge;
static string lastAntibodyChains;

antigenInfo getAntigenInfos(std::string ID)
{

    std::pair<superProtein *, vector<int>> AG = getAntigen(ID);

    antigenInfo res;
    res.ID = ID;
    res.first = AG.first;
    res.second = AG.second;
    res.glycans = lastSearchGlycan;
    res.hotspotsCore = lastSearchHotspotsCore;
    res.hotspotsLarge = lastSearchHotspotsLarge;
    res.antibodyChains = lastAntibodyChains;
    return res;
}

antigenInfo::antigenInfo()
{
    first = nullptr;
    ID = "?";
}

string antigenInfo::print()
{
    if (!first)
        return string("ID") + ID + string(", Empty Antigen + Infos");
    stringstream res;
    res << printProtein(*first) << endl;
    res << "Forbidden positions: ";
    res << printVector(second) << endl;
    res << "Glycans=" << printVector(glycans) << endl;
    size_t nHotspots = hotspotsCore.size();
    if (hotspotsLarge.size() != nHotspots)
    {
        cerr << "ERR: antigenInfo::print, for antigen " << ID << ", different number of nb of hotspots for Core and Large positions" << endl;
        return string("ERR");
    }
    for (size_t i = 0; i < nHotspots; ++i)
    {
        res << "Hotspot " << ID << "_H" << i + 1 << "\tCore=\t" << printVector(hotspotsCore[i]) << "\tBound 100%:\t" << printVector(hotspotsLarge[i]) << endl;
    }
    res << "Known thresholds:" << printVector(thresholds) << endl;
    return res.str();
}

std::pair<superProtein *, vector<int>> getAntigen(std::string ID, bool display)
{

    string agStruct;
    string agSeq;
    bool multichain = false;
    int startPos = -1;
    vector<int> blockV;
    vector<int> glycans; // this is just for info, they should be included inside blockV as well.
    vector<vector<int>> hotspotsCore;
    vector<vector<int>> hotspotsLarge;
    string antibodyChains = "LH"; // most of them are. If different, will be rewritten
    // multi-chains proteins
    vector<int> startPosChains;
    vector<string> structChains;

    lastSearchGlycan.clear();
    lastSearchHotspotsCore.clear();
    lastSearchHotspotsLarge.clear();
    lastAntibodyChains.clear();

    int codeID = str2int(ID.c_str());

    try
    {
        auto data = toml::parse("antigenlib.toml");
        auto antigens = toml::find<toml::table>(data, "antigens");

        for (const auto &[key, value] : antigens)
        {
            if (str2int(key.c_str()) == codeID)
            {
                multichain = toml::find<bool>(value, "multichain");
                if (!multichain)
                {
                    agStruct = toml::find<string>(value, "agStruct");
                    agSeq = toml::find<string>(value, "agSeq");
                    startPos = toml::find<int>(value, "startPos");
                    blockV = toml::find<vector<int>>(value, "blockV");
                    if (value.contains("hotspotsCore"))
                        hotspotsCore = toml::find<vector<vector<int>>>(value, "hotspotsCore");
                    if (value.contains("hotspotsLarge"))
                        hotspotsLarge = toml::find<vector<vector<int>>>(value, "hotspotsLarge");
                    if (value.contains("antibodyChains"))
                        antibodyChains = toml::find<string>(value, "antibodyChains");

                    if ((loaded == true) && (agStruct.size() > 0))
                    {
                        if (display)
                            displayLigand(agStruct, startPos, blockV);
                        superProtein *sp = new superProtein(agStruct, startPos);
                        sp->setAAs(agSeq);

                        lastSearchGlycan = glycans;
                        lastSearchHotspotsCore = hotspotsCore;
                        lastSearchHotspotsLarge = hotspotsLarge;
                        lastAntibodyChains = antibodyChains;

                        return std::pair<superProtein *, vector<int>>(sp, blockV);
                    }
                }
                else
                {
                    startPosChains = toml::find<vector<int>>(value, "startPosChains");
                    structChains = toml::find<vector<string>>(value, "structChains");
                    agSeq = toml::find<string>(value, "agSeq");
                    blockV = toml::find<vector<int>>(value, "blockV");
                    if (value.contains("glycans"))
                    {
                        glycans = toml::find<vector<int>>(value, "glycans");
                    }
                    if (value.contains("hotspotsCore"))
                    {
                        hotspotsCore = toml::find<vector<vector<int>>>(value, "hotspotsCore");
                    }
                    if (value.contains("hotspotsLarge"))
                    {
                        hotspotsLarge = toml::find<vector<vector<int>>>(value, "hotspotsLarge");
                    }
                    if (value.contains("antibodyChains"))
                    {
                        antibodyChains = toml::find<string>(value, "antibodyChains");
                    }

                    size_t N = startPosChains.size();
                    if (N > 0)
                    {
                        if (structChains.size() != N)
                        {
                            cerr << "ERR: getAntigen(" << ID << "), got different number of starting positions and chains structures" << endl;
                        }
                        else
                        {
                            superProtein *P1 = new superProtein(structChains[0], startPosChains[0]);
                            for (size_t i = 1; i < N; ++i)
                            {
                                superProtein *P2 = new superProtein(insert(P1, structChains[i], startPosChains[i], 1000 * i));
                                delete P1;
                                P1 = P2;
                            }
                            if (agSeq.size() != P1->size())
                            {
                                cerr << "ERR:  getAntigen(" << ID << "), contains " << N << " chains with total number of AAs=" << P1->size() << " while " << agSeq.size() << " AAs are provided. Doesn't match." << endl;
                            }
                            else
                            {
                                P1->setAAs(agSeq);
                                if (display)
                                    displayLigand(P1, blockV);

                                lastSearchGlycan = glycans;
                                lastSearchHotspotsCore = hotspotsCore;
                                lastSearchHotspotsLarge = hotspotsLarge;
                                lastAntibodyChains = antibodyChains;

                                return std::pair<superProtein *, vector<int>>(P1, blockV);
                            }
                        }
                    }
                    else
                    {
                        if (ID.compare("NotToBeFoundID"))
                        {
                            cerr << "ERR: getAntigen(ID=" << ID << "), couldn't find this antigen. Call listIDs() for the list of existing antigen IDs.\n";
                        }
                    }
                }
                break;
            }
        }
    }
    catch (const toml::exception &e)
    {
        cerr << "Error parsing antigenlib.toml: " << e.what() << endl;
    }
    return std::pair<superProtein *, vector<int>>(nullptr, vector<int>());
}

affinityOneLigand *getAffinityAntigen(std::string ID, int sizeReceptors, int minInteract)
{

    std::pair<superProtein *, vector<int>> antigen = getAntigen(ID);
    if (antigen.first == nullptr)
        return nullptr;
    affinityOneLigand *letsGo = new affinityOneLigand(antigen.first, sizeReceptors, minInteract, -1, 1.0, antigen.second);
    return letsGo;
}
