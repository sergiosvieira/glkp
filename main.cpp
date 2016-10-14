#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <map>
#include <string>
#include <glpk.h>

//#define DEBUG

using namespace std;

/* Fowarding Class */
class Demand;
class Link;
class Element;

/* Constantes, Tipos e Enumerações */
static const int kMaxPriority = 9999;
static const char* kReservoir_RelZone = "Z";// Volume liberado pela zona de prioridade
//static const char* kReservoir_RecZone = "Z";// Volume liberado pela zona de prioridade
static const char* kDemand = "D"; // Demanda atendida
static const char* kMaximumVolume = "VMAX"; // Volume máximo suportado
static const char* kMinimumVolume = "VMIN"; // Volume mínimo suportado
static const char* kCurrentVolume = "VCUR"; // Voluma atual
static const char* kSilth = "V"; // Volume vertido pelo reservatório
static const char* kExcess = "E"; // Volume recebido a mais (excesso)
static const char* kLink = "L"; // Volume recebido pelo trecho
static const char* kDemandSubject = "r_demands";
static const char* kReservoirDemand = "R"; // Volume recebido por reservatório oriundo
                                           // de outro a montante
static const char* kExcessVar = "EV"; // Variável genérica de excesso

enum class Type
{
    RESERVOIR,
    LINK,
    DEMAND
};

using MapOfDemand = std::map<int, double>; // int = identificador do reservatorio, double = volume total a jusante do reservatório
using MapOfZones = std::map<int, double>;
using VectorOfElement = std::vector<Element*>;
using MapOfElement = std::map<int, Element*>;

/* Classes */
class Element
{
public:
    explicit Element(int id_): id(id_){};
    virtual ~Element(){};
    Element(){};
    int id = 0;
    Type type = Type::RESERVOIR;
    std::vector<Link> downLinks;
    std::vector<Link> upLinks;
};

class Link : public Element
{
public:
    explicit Link(int id): Element(id)
    {
        type = Type::LINK;
    };
    Element* from = nullptr;
    Element* to = nullptr;
    double max = 0.0;
    double min = 0.0;
};

class PriorityZone
{
public:
    PriorityZone(int p, double v): priority(p), upper_bound(v) {};
    int priority = kMaxPriority;
    double upper_bound = 0.0;
    friend std::ostream& operator<<(std::ostream& os, const PriorityZone& obj)
    {
        os << "["
           << obj.priority
           << "->"
           << obj.upper_bound
           << "]";
        return os;
    }
};

class Reservoir : public Element
{
public:
    Reservoir(int id, double max, double min, double cur): Element(id),
                                                           max_volume(max),
                                                           min_volume(min),
                                                           cur_volume(cur)
    {
        type = Type::RESERVOIR;
    };
    double max_volume = 0.0;
    double min_volume = 0.0;
    double cur_volume = 0.0;
    double silth = 0.0; // o que verte
    double excess = 0.0;
    double demand = 0.0;
    std::vector<PriorityZone> release_zones;
    std::vector<PriorityZone> receiving_zones;
};

class Demand: public Element
{
public:
    Demand(int id, int p, double v): Element(id),
                                     value(v),
                                     priority(p)
    {
        type = Type::DEMAND;
    }
    double value = 0.0;
    int priority = kMaxPriority;
};

/* Métodos de Apoio */
std::string subject(glp_prob* lp)
{
    std::string result = "r_" + std::to_string(glp_get_num_rows(lp));
    return result;
}

void debug(glp_prob* lp, const char* method, int rowIndex)
{
#ifdef DEBUG
    std::string name_ = std::string(method);
    name_ += subject(lp);
    glp_set_row_name(lp, rowIndex, name_.c_str());
#endif
}

std::string fromName(const char* variable, int x)
{
    std::string name = variable;
    name += "_";
    name += std::to_string(x);
    return name;
}

std::string fromToName(const char* variable, int x, int y)
{
    std::string name = fromName(variable, x);
    name += "_";
    name += std::to_string(y);
    return name;
}

double demand(Element * e)
{
    if (e == nullptr) return 0.0;
    if (e->type == Type::DEMAND)
    {
        Demand* d = dynamic_cast<Demand*>(e);
        return d->value;
    }
    double result = 0.0;
    for (const Link& l : e->downLinks)
    {
        Element* to = l.to;
        result += demand(to);
    }
    return result;
}

MapOfDemand makeMapOfDemand(MapOfElement& elements)
{
    MapOfDemand result;
    for (auto keyValue: elements)
    {
        int id = keyValue.first;
        Element* e = keyValue.second;
        if (e->type == Type::RESERVOIR)
        {
            double d = demand(e);
            result[id] = d;
        }
    }
    return result;
}

template <typename T>
void print(std::vector<T>& array)
{
    std::cout << "(";
    int counter = 0;
    for (const T& v: array)
    {
        std::cout << v;
        if (counter++ < array.size() - 1)
        {
            std::cout << ", ";
        }
    }
    std::cout << ")\n";
}

// Método de propósito geral para adicionar novas colunas no modelo
int colValue(glp_prob* lp, const std::string& name, int type, double lower = 0.0, double upper = 0.0)
{
    int result = 0;
    if (glp_get_num_cols(lp) != 0)
    {
        result = glp_find_col(lp, name.c_str());
    }
    // caso a variável não tenha sido criada, crie-a
    if (result == 0)
    {
        glp_create_index(lp);
        result = glp_add_cols(lp, 1);
        glp_set_col_name(lp, result, name.c_str());
        glp_set_col_bnds(lp, result, type, lower, upper);
    }
    return result;
}

int rowValue(glp_prob* lp, const std::string& name, int type, double lower = 0.0, double upper = 0.0)
{
    int result = 0;
    if (glp_get_num_rows(lp) != 0)
    {
        result = glp_find_row(lp, name.c_str());
    }
    // caso a variável não tenha sido criada, crie-a
    if (result == 0)
    {
        glp_create_index(lp);
        result = glp_add_rows(lp, 1);
        glp_set_row_name(lp, result, name.c_str());
        glp_set_row_bnds(lp, result, type, lower, upper);
    }
    return result;
}

// GLP_FX, value, value;
int assignColValue(glp_prob* lp, const std::string& name, double value)
{
    return colValue(lp, name, GLP_FX, value, value);
}

// GLP_FX, value, value;
int assignRowValue(glp_prob* lp, const std::string& name, double value)
{
    return rowValue(lp, name, GLP_FX, value, value);
}

// GLP_LO, 0.0, 0.0;
int lowColValue(glp_prob* lp, const std::string& name, double value = 0.0)
{
    return colValue(lp, name, GLP_LO, value, 0.0);
}

// GLP_LO, 0.0, 0.0;
int lowRowValue(glp_prob* lp, const std::string& name, double value = 0.0)
{
    return rowValue(lp, name, GLP_LO, value, 0.0);
}

// GLP_UP, 0.0, value;
int upRowValue(glp_prob* lp, const std::string& name, double value)
{
    return rowValue(lp, name, GLP_UP, 0.0, value);
}

// GLP_UP, 0.0, value;
int upColValue(glp_prob* lp, const std::string& name, double value)
{
    return colValue(lp, name, GLP_UP, 0.0, value);
}

// GLP_DB, lower, upper
int upDownRowValue(glp_prob* lp, const std::string& name, double lower, double upper)
{
    return rowValue(lp, name, GLP_DB, lower, upper);
}

// GLP_DB, lower, upper
int upDownColValue(glp_prob* lp, const std::string& name, double lower, double upper)
{
    return colValue(lp, name, GLP_DB, lower, upper);
}

void makeZonesSubject(glp_prob* lp, double total, int total_zones)
{
    int index = lowRowValue(lp, subject(lp), total);
    debug(lp, "MakeZonesSubject_", index);
    std::vector<int> indices(total_zones);
    std::iota(indices.begin(), indices.end(), 0);// 0, 1, 2, ..., n - 1
    std::vector<double> values(indices.size(), 1.0);
    values[0] = 0.;
    glp_set_mat_row(lp, index, indices.size() - 1, &indices[0], &values[0]);
}

void linkVolumeSubject(glp_prob* lp,
                       int rowIndex,
                       Link& l,
                       double volumeValue,
                       const std::string& name)
{
    // A variável name representa VMIN ou VMAX
    int colIndex = assignColValue(lp, name, volumeValue);
    int linkIndex = lowColValue(lp, fromToName(kLink, l.from->id, l.to->id));
    glp_set_mat_row(lp, rowIndex, 2, &std::vector<int>{0, colIndex, linkIndex}[0], &std::vector<double>{0.0, -1.0, 1.0}[0]);
    debug(lp, "LinkVolumeSubject_", rowIndex);
}


void makeLinksSubject(glp_prob* lp, Reservoir& r)
{
    for (Link& link: r.downLinks)
    {
        int from = link.from->id;
        int to = link.to->id;
        if (link.max > 0)
        {
            int maxIndex = upRowValue(lp, subject(lp), 0.0);
            linkVolumeSubject(lp,
                              maxIndex,
                              link,
                              link.max,
                              fromToName(kMaximumVolume, from, to));
        }
        if (link.min >= 0)
        {
            int minIndex = lowRowValue(lp, subject(lp));
            linkVolumeSubject(lp,
                              minIndex,
                              link,
                              link.min,
                              fromToName(kMinimumVolume, from, to));
        }
    }
}

int makeSilthSubject(glp_prob* lp, int id, int maxIndex, int curIndex)
{
    int index = lowRowValue(lp, subject(lp));
    int silthIndex = lowColValue(lp, fromName(kSilth, id));
    std::vector<int> indices = {0, maxIndex, curIndex, silthIndex};
    std::vector<double> values = {0, 1, -1, 1};
    glp_set_mat_row(lp, index, 3, &indices[0], &values[0]);
    debug(lp, "MakeSilthSubject", index);
    return silthIndex;
}

int makeExcessSubject(glp_prob* lp, int id, int silthIndex)
{
    int index = lowRowValue(lp, subject(lp));
    int excIndex = lowColValue(lp, fromName(kExcess, id));
    std::vector<int> indices = {0, silthIndex, excIndex};
    std::vector<double> values = {0, 1, -1};
    glp_set_mat_row(lp, index, 2, &indices[0], &values[0]);
    debug(lp, "MakeExcessSubject", index);
    return excIndex;
}

int makeLinkVolumeSubject(glp_prob* lp, Reservoir& r, std::vector<int>& indices)
{
    std::vector<int> aux(indices);
    std::vector<double> values(indices.size(), 1.0);
    values[0] = 0.0; // O primeiro valor sempre é zero (característica do API do GLPK, olhe o manual)
    int index = assignRowValue(lp, subject(lp), 0.0);
    for (Link& l: r.downLinks)
    {
        int colIndex = lowColValue(lp, fromToName(kLink, l.from->id, l.to->id));
        aux.push_back(colIndex);
        values.push_back(-1.0);
    }
    glp_set_mat_row(lp, index, aux.size() - 1, &aux[0], &values[0]);
    debug(lp, "MakeLinkVolumeSubject_", index);
}

int receivingPriority(Reservoir& r)
{
    int result = 0;
    auto item = std::find_if(r.receiving_zones.begin(), r.receiving_zones.end(), [=](PriorityZone const& a){
        return r.cur_volume < a.upper_bound;
    });
    if (item != r.receiving_zones.end())
    {
        result = (*item).priority;
    }
    return result;
}

std::vector<PriorityZone> updateZones(std::vector<PriorityZone> zones,
                                      double minimumVolume,
                                      double currentVolume)
{
    /*
    int counter = 0;
    std::vector<PriorityZone> aux(zones);
    std::sort(zones.begin(), zones.end(), [](PriorityZone const& a, PriorityZone const & b) -> bool {
        return a.upper_bound < b.upper_bound;
    });
    auto item = std::find_if(zones.begin(), zones.end(), [=](PriorityZone const& a) -> bool {
        return currentVolume <= a.upper_bound;
    });
    if (item == zones.end())
    {
        // Neste caso, o volume atual do reservatório é maior que o volume das zonas de prioridade configuradas
        // Portanto, os volumes serão subtraídos do volume anterior, a partir da segunda zona
        // Here, the current volume of the reservoir is greater than the volume of the configured priorities zones
        // So, all zone's value will be subtracted with the previous zone's value, from second zone to the last zone
        if (zones.size() > 0)
        {
            zones[0].upper_bound -= minimumVolume;
        }
        for (int i = 1; i < zones.size(); ++i)
        {
            zones[i].upper_bound -= aux[i - 1].upper_bound;
        }
    }
    else
    {
        // Neste caso existem duas coisas a serem feitas
        // 1. Zerar os valores das zonas com volumes acima do volume do reservatório inicial
        std::for_each(item + 1, zones.end(), [](PriorityZone & a) {
            a.upper_bound = 0.0;
        });
        // 2. Calcular os valores dos volumes de cada zona com volume menor ou igual ao volume inicial do reservatório
        int lastIndex = item - zones.begin();
        if (zones.size() > 0)
        {
            if (currentVolume > zones[0].upper_bound)
            {
                zones[0].upper_bound -= minimumVolume;
            }
            else
            {
                zones[0].upper_bound = currentVolume - minimumVolume;
            }
        }
        for (int i = 1; i <= lastIndex; ++i)
        {
            zones[i].upper_bound = currentVolume - aux[i - 1].upper_bound;
        }
    }
    */
    std::vector<PriorityZone> aux(zones);
    std::sort(zones.begin(), zones.end(), [](PriorityZone const& a, PriorityZone const & b) -> bool {
        return a.upper_bound < b.upper_bound;
    });
    for (int i = 1; i < zones.size(); ++i)
    {
        zones[i].upper_bound = aux[i].upper_bound - aux[i - 1].upper_bound;
    }
    return zones;
}

/** As vezes, quando um reservatório está ligado a um outro,
 *  o mais a jusante deve se comportar como uma demanda **/
void modelReservoirDemand(glp_prob* lp,
                          Reservoir& r,
                          double total_demands,
                          int priority)
{
    /**
     * TO DO:
     *  ModelReservoirDemand_r_30: + R_12 - VMAX_13_12 = 0 - Demanda é maior que o volume máximo
     *  ModelReservoirDemand_r_33: + R_12 - Z_12_2 - Z_12_1 = 0 - Demanda é menor que o volume máximo
     * Demanda deve ser o mínimo da vazão máxima e do volume util
    **/
    int reservoirIndex = lowColValue(lp, fromName(kReservoirDemand, r.id));
    glp_set_obj_coef(lp, reservoirIndex, priority);
    int rowIndex = lowRowValue(lp, subject(lp), total_demands);
    glp_set_mat_row(lp, rowIndex, 1, &std::vector<int>{0, reservoirIndex}[0], &std::vector<double>{0.0, 1.0}[0]);
    int excessIndex = glp_find_col(lp, fromName(kExcess, r.id).c_str());
    if (excessIndex != 0)
    {
        glp_set_obj_coef(lp, excessIndex, - 1.0);
    }
    debug(lp, "ModelReservoirDemand_", rowIndex);
}

void modelReservoir(glp_prob* lp,
                    Reservoir& r,
                    double total_demands)
{
    // Objective Function
    int counter = 1;
    std::vector<int> indices;
    indices.push_back(0);// o primeiro valor sempre é zero
    std::vector<PriorityZone> rel_zones = updateZones(r.release_zones, r.min_volume, r.cur_volume);
    double epsilon = 0.00001;
    for (const PriorityZone& pz: r.release_zones)
    {
        int index = lowColValue(lp, fromToName(kReservoir_RelZone, r.id, counter));
        indices.push_back(index);
        glp_set_obj_coef(lp, index, -pz.priority);
        // Tratar limites das zonas
        auto item = std::find_if(rel_zones.begin(), rel_zones.end(), [=](PriorityZone const& a){
            return a.priority == pz.priority;
        });
        double upperBound = (*item).upper_bound;
        // usado para testar se upperBound é igual a zero
        if (std::fabs(upperBound) <= epsilon * std::fabs(upperBound))
        {
            glp_set_col_bnds(lp, index, GLP_FX, 0.0, 0.0);
        }
        else
        {
            glp_set_col_bnds(lp, index, GLP_DB, 0.0, upperBound);
        }
        ++counter;
    }
    // Tratar valores liberados para cada Link, caso existam links a downsteam
    if (r.downLinks.size() > 0)
    {
        makeLinkVolumeSubject(lp, r, indices);
    }
    // Tratar Liberações considerando o total de demandas quando a demanda for maior que zero
//    if (std::fabs(total_demands) > epsilon * std::fabs(total_demands))
//    {
        //makeZonesSubject(lp, total_demands, r.release_zones.size() + 1);
    makeZonesSubject(lp, 0.0, r.release_zones.size() + 1);
//    }
    // Tratar Liberações para os links considerando seus volumes máximo e mínimos
    makeLinksSubject(lp, r);
    // Tratar vertimentos e excessos
    /* Volume Máximo */
    int maxIndex = assignColValue(lp, fromName(kMaximumVolume, r.id), r.max_volume);
    /* Volume Mínimo */
    int minIndex = assignColValue(lp, fromName(kMinimumVolume, r.id), r.min_volume);
    int minRowIndex = assignRowValue(lp, subject(lp), r.min_volume);
    // Tô fazendo isso só para retirar o volume mínimo da função objetivo
    glp_set_mat_row(lp, minRowIndex, 1, &std::vector<int>{0, minIndex}[0], &std::vector<double>{0, 1.0}[0]);
    debug(lp, "MinimumVolume_", minRowIndex);
    /* Volume Atual */
    int curIndex = assignColValue(lp, fromName(kCurrentVolume, r.id), r.cur_volume);
    /* Tratar o quer verte */
    int silthIndex = makeSilthSubject(lp, r.id, maxIndex, curIndex);
    /* Tratar o Excesso */
    makeExcessSubject(lp, r.id, silthIndex);
}

int modelDemand(Demand& d,
                glp_prob* lp,
                double min = 0.0)
{
    // Adicionar variável de demanda atendida na função objetivo (maximização)
    //double total = std::max(min, d.value);
    double total = d.value;
    int demandIndex = upDownColValue(lp, fromName(kDemand, d.id), 0.0, total);
    glp_set_obj_coef(lp, demandIndex, d.priority);
    // Adicionar excesso da demanda atendida na função objetivo (minimização)
    int excessIndex = lowColValue(lp, fromName(kExcess, d.id));
    glp_set_obj_coef(lp, excessIndex, -1.0);
    return demandIndex;
}

void modelLink(Link& l,
               glp_prob* lp)
{
    if (l.to != nullptr)
    {
        int from = l.from->id;
        int to = l.to->id;
        std::string rdName = fromName(kDemand, to);
        if (l.to->type == Type::RESERVOIR)
        {
            rdName = fromName(kReservoirDemand, to);
        }
        int linkIndex = lowColValue(lp, fromToName(kLink, from, to));
        int demandIndex = glp_find_col(lp, rdName.c_str());
        int excessIndex = glp_find_col(lp, fromName(kExcess, to).c_str());
        int index = assignRowValue(lp, subject(lp), 0.0);
        std::vector<int> indices = {0, linkIndex, demandIndex, excessIndex};
        std::vector<double> values = {0.0, 1.0, -1.0, -1.0};
        glp_set_mat_row(lp, index, indices.size() - 1, &indices[0], &values[0]);
        debug(lp, "ModelLink_", index);
    }
}

int main(int argc, char *argv[])
{
    glp_prob *lp = glp_create_prob();
    glp_set_prob_name(lp, "teste01");
    glp_set_obj_dir(lp, GLP_MAX);
    glp_add_rows(lp, 1);
    // Criando rede hidrológica
    /** Reservatório 12 **/
    Reservoir r12(12, 15000.0, 2000.0, 7000.0);
    r12.receiving_zones.push_back(PriorityZone(1, 15000.0));
    r12.receiving_zones.push_back(PriorityZone(4, 7000.0));
    r12.release_zones.push_back(PriorityZone(8, 6500));
    r12.release_zones.push_back(PriorityZone(4, 10000));
    r12.release_zones.push_back(PriorityZone(1, 15000));
    /** Demanda 1 **/
    Demand d1(1, 8, 26.784);
    /** Link do reservatório 12 a demanda 1 **/
    Link link12_1(121);
    link12_1.from = &r12;
    link12_1.to = &d1;
    link12_1.max = 5000.0;
    link12_1.min = 10.0;
    r12.downLinks.push_back(link12_1);
    /** Demanda 2 **/
    Demand d2(2, 4, 26.784);
    Link link12_2(122);
    link12_2.from = &r12;
    link12_2.to = &d2;
    link12_2.max = 5000.;
    link12_2.min = 10.0;
    r12.downLinks.push_back(link12_2);
    /** Demanda 3 **/
    Demand d3(3, 1, 26.784);
    Link link12_3(123);
    link12_3.from = &r12;
    link12_3.to = &d3;
    link12_3.max = 100.0;
    link12_3.min = 50.0;
    r12.downLinks.push_back(link12_3);
    /** Mapa de Elementos **/
    MapOfElement elements;
    elements[r12.id] = &r12;
    // Modelando rede hidrológica
    MapOfDemand demands = makeMapOfDemand(elements);
    for (auto kv : elements)
    {
        Element* e = kv.second;
        if (e->type == Type::RESERVOIR)
        {
            std::vector<int> demandIndices;
            demandIndices.push_back(0); // sempre começa com zero (vide manual GLPK)
            Reservoir& r = *dynamic_cast<Reservoir*>(e);
            double min = std::numeric_limits<double>::max();
            MapOfDemand::const_iterator it = demands.find(r.id);
            double total_demands = (it != demands.end()) ? (*it).second : 0.0;
            modelReservoir(lp, r, total_demands);
            for (Link& l: r.downLinks)
            {
                if (l.to->type == Type::DEMAND)
                {
                    Demand& d = *dynamic_cast<Demand*>(l.to);
                    min = std::min(min, d.value);
                    demandIndices.push_back(modelDemand(d, lp, l.min));
                }
                else if (l.to->type == Type::RESERVOIR)
                {
                    std::vector<PriorityZone> rel_zones = updateZones(r.release_zones, r.min_volume, r.cur_volume);
                    double util_volume = 0.0;
                    for (PriorityZone & pz : rel_zones)
                    {
                        util_volume += pz.upper_bound;
                    }
                    double total = std::min(min, std::min(util_volume, total_demands));
                    Reservoir& to_r = *dynamic_cast<Reservoir*>(l.to);
                    modelReservoir(lp, to_r, total_demands);
                    int priority = receivingPriority(to_r);
                    total_demands = std::min(total_demands, l.max);
                    total_demands = std::min(total_demands, total);
                    modelReservoirDemand(lp, to_r, total_demands, priority);
                }
                modelLink(l, lp);
            }
            // Adiciona restrição de atendimento da demanda
            if (demandIndices.size() > 1)
            {
                int excessIndex = lowColValue(lp, fromName(kExcessVar, e->id), 0.0);
                demandIndices.push_back(excessIndex);
                std::vector<PriorityZone> rel_zones = updateZones(r.release_zones, r.min_volume, r.cur_volume);
                double util_volume = 0.0;
                for (PriorityZone & pz : rel_zones)
                {
                    util_volume += pz.upper_bound;
                }
                //double total = std::min(min, std::min(util_volume, total_demands));
                double total = std::min(util_volume, total_demands);
                int index = assignRowValue(lp, kDemandSubject, total);
                std::vector<double> values(demandIndices.size(), 1.0);
                values[0] = 0.0; // sempre começa com zero (vide manual GLPK)
                values.push_back(1.0);
                glp_set_mat_row(lp, index, demandIndices.size() - 1, &demandIndices[0], &values[0]);
                debug(lp, "DemandReleaseSubject", index);
            }
        }
    }
    // Executando Modelo
//    glp_iocp parm;
//    glp_intopt(lp, nullptr);
    glp_write_lp(lp, nullptr, "./model.txt");
    int varIndex = glp_get_unbnd_ray(lp);

    std::cout << "Unbound: " << varIndex << "\n";
    glp_smcp parm;
    glp_init_smcp(&parm);
    //parm.meth = GLP_DUAL;
    int result = glp_simplex(lp, &parm);
    glp_print_sol(lp, "solution.txt");
    // Finalizando
    glp_delete_prob(lp);
    glp_free_env();
    return 0;
}
