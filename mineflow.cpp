/* mineflow.cpp
Contact: matthewvdeutsch@gmail.com
You should also have 'mineflow.h', or else you'll have trouble. See that file
for more information.

Citation:
@article{deutsch2022mineflow,
  author={Deutsch, Matthew and Da{\u{g}}delen, Kadri and Johnson, Thys},
  title={An Open-Source Program for Efficiently Computing Ultimate Pit Limits: MineFlow},
  journal={Natural Resources Research},
  year={2022},
  month={Mar},
  day={17},
  issn={1573-8981},
  doi={10.1007/s11053-022-10035-w},
  url={https://doi.org/10.1007/s11053-022-10035-w}
}

This file contains roughly 3 sections:

    - The actual implementation of the library
    - A testing framework and tests: Access with compile definition: MVD_MINEFLOW_TESTS
    - The executable: Access with compile definition: MVD_MINEFLOW_EXE

It is not great practice to recompile the same implementation multiple times
with different #defines to get at what you want, but note that the decision here
was to optimize for integration, not elegance.


LICENSE
-------
Copyright 2022 Matthew Deutsch

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <iomanip>
#include <functional>
#include <numeric>
#include <chrono>
#include <unordered_map>
#include <unordered_set>

#ifdef MVD_USE_GMP
#include <forward_list>
#include <cstring>
#endif

#include "mineflow.h"

using namespace mvd::mineflow;
using namespace mvd::mineflow::impl;

#ifndef MVD_ASSERT
    #ifdef NDEBUG
        #undef NDEBUG
    #endif

    #include <cassert>
    #define MVD_ASSERT(x) assert(x)
#endif

namespace mvd::mineflow::impl {
class NodePool;
struct AntecedentsInfo {
    std::vector<Node*> OutOfTree;
    IndexType NextArc;
    NodePool* Init;
};

struct Node {
    ValueType Excess; // Positive for excess, negative for deficit
    Arc* ToRoot;    // Normalized Tree arcs
    IndexType Label; // The 'distance' label

    Node* FirstChild; 
    Node* NextChild;
    Node* NextScan;
    AntecedentsInfo Antecedents;


    // Operations
    void AddChild(Node* child);
    void RemoveChild(Node* child);
    void IncrementLabel();
    void ForNodeAndChildren(std::function<void(Node*)> cback);
    Node* FindWeakAbove();
    void InitPrecedence(Node* node);
};
struct Arc {
    Node* Tail;     // null for main 'root'
    Node* Head;     // null for main 'root'
    ValueType Flow;
};

class PrecedenceArcPool {
public:
    PrecedenceArcPool();
    ~PrecedenceArcPool();

    Arc* NewArc(Node* from, Node* to);
    void DeleteArc(Arc* arc);

    IndexType NumUsed() const;

private:
#ifndef MVD_USE_GMP
    ObjectPoolBase<8192, 16> m_ObjectPool; // tunable
#else
    template <size_t N>
    struct ArcSet {
        ArcSet() {
            for (size_t i = 0; i < N; i++) {
                mpz_init(Arcs[i].Flow);
            }
        }
        ~ArcSet() {
            for (size_t i = 0; i < N; i++) {
                mpz_clear(Arcs[i].Flow);
            }
        }

        Arc Arcs[N];
    };
    static constexpr inline int N = 1024;
    std::forward_list<ArcSet<N>> m_ArcSets; // tunable
    int m_Remaining;
    ArcSet<N>* m_Set;
#endif
    IndexType m_NumUsed;
};

class NodePool {
public:
    NodePool(std::shared_ptr<const IPrecedenceConstraints> pre);
    ~NodePool();

    void InitializeNodeValue(IndexType nodeIndex, std::function<void(ValueType*)> getValue);

    void GetNodeValue(IndexType nodeIndex, ValueType* value) const;
    Node* GetNode(IndexType nodeIndex);
    IndexType GetNodeIndex(const Node* node) const;

    void ReconnectToRoot(Node* node);

    void PushStrongRoot(Node* node);
    bool NextStrongRoot(Node** nodep);

    void IncrementLabel(Node* node);

    // Output
    IndexType NumNodes() const;
    bool InMinimumCut(IndexType nodeIndex) const;
    void InitPrecedence(Node* node);


private:
    std::shared_ptr<const IPrecedenceConstraints> m_PrecedenceConstraints;
    IndexType m_NumNodes;

    std::vector<IndexType> m_LabelCount;
    std::vector<std::queue<Node*>> m_Buckets;

    std::vector<Node> m_Nodes;
    std::vector<Arc> m_RootArcs;
};
}

IBlockValues::IBlockValues()
{
}

IBlockValues::~IBlockValues() {
}

////////////////////////////////////////////////////////////////////////////////

#ifndef MVD_USE_GMP
VecBlockValues::VecBlockValues(IndexType numBlocks)
{
    m_Values.assign(numBlocks, 0);
}

VecBlockValues::VecBlockValues(std::vector<ValueType>&& values)
    : m_Values(values)
{
}

VecBlockValues::VecBlockValues(std::initializer_list<int> values)
{
    m_Values.reserve(values.size());
    for (auto & v : values) {
        m_Values.push_back(static_cast<ValueType>(v));
    }
}

VecBlockValues::VecBlockValues(IBlockValues* valuesToCopy)
{
    mineflow::IndexType n = valuesToCopy->NumBlocks();

    m_Values.resize(n);
    for (mineflow::IndexType blockIndex = 0; blockIndex < n; blockIndex++) {
        valuesToCopy->BlockValue(blockIndex, &m_Values[blockIndex]);
    }
}

VecBlockValues::~VecBlockValues()
{
}

IndexType VecBlockValues::NumBlocks() const
{
    return m_Values.size();
}

void VecBlockValues::BlockValue(IndexType blockIndex, ValueType* value) const
{
    *value = m_Values.at(blockIndex);
}

VecBlockValues::const_iterator VecBlockValues::begin() const
{
    return m_Values.begin();
}

VecBlockValues::const_iterator VecBlockValues::end() const
{
    return m_Values.end();
}

VecBlockValues::iterator VecBlockValues::begin()
{
    return m_Values.begin();
}

VecBlockValues::iterator VecBlockValues::end()
{
    return m_Values.end();
}

ValueType VecBlockValues::operator[](IndexType blockIndex) const
{
    return m_Values[blockIndex];
}

ValueType& VecBlockValues::operator[](IndexType blockIndex)
{
    return m_Values[blockIndex];
}

void VecBlockValues::SetBlockValueSI(IndexType blockIndex, int64_t si)
{
    m_Values[blockIndex] = static_cast<ValueType>(si);
}

#else

GMPBlockValues::GMPBlockValues(IndexType numBlocks)
    : m_NumBlocks(numBlocks)
    , m_BlockValues(std::make_unique<ValueType[]>(numBlocks))
{
    for (IndexType i = 0; i < m_NumBlocks; i++) {
        mpz_init(m_BlockValues[i]);
    }
}

GMPBlockValues::GMPBlockValues(IndexType numBlocks, int64_t initialValue)
    : m_NumBlocks(numBlocks)
    , m_BlockValues(std::make_unique<ValueType[]>(numBlocks))
{
    for (IndexType i = 0; i < m_NumBlocks; i++) {
        mpz_init_set_si(m_BlockValues[i], initialValue);
    }
}

GMPBlockValues::GMPBlockValues(const std::vector<int64_t>& initialValues)
    : GMPBlockValues(initialValues.size())
{
    for (IndexType i = 0; i < m_NumBlocks; i++) {
        SetBlockValueSI(i, initialValues[i]);
    }
}

GMPBlockValues::~GMPBlockValues()
{
    for (IndexType i = 0; i < m_NumBlocks; i++) {
        mpz_clear(m_BlockValues[i]);
    }
}

IndexType GMPBlockValues::NumBlocks() const
{
    return m_NumBlocks;
}

void GMPBlockValues::BlockValue(IndexType blockIndex, ValueType* value) const
{
    mpz_set(*value, m_BlockValues[blockIndex]);
}

void GMPBlockValues::SetBlockValueSI(IndexType blockIndex, int64_t si)
{
    mpz_set_si(m_BlockValues[blockIndex], si);
}

#endif

////////////////////////////////////////////////////////////////////////////////

BlockDefinition::BlockDefinition()
{
}

BlockDefinition::~BlockDefinition()
{
}

BlockDefinition::BlockDefinition(
        IndexType iNumX, IndexType iNumY, IndexType iNumZ,
        double iMinX, double iMinY, double iMinZ,
        double iSizeX, double iSizeY, double iSizeZ)
    : NumX(iNumX), NumY(iNumY), NumZ(iNumZ)
    , MinX(iMinX), MinY(iMinY), MinZ(iMinZ)
    , SizeX(iSizeX), SizeY(iSizeY), SizeZ(iSizeZ)
{
}

IndexType BlockDefinition::GridIndex(IndexType x, IndexType y, IndexType z) const
{
    return x + y * NumX + z * NumX * NumY;
}

IndexType BlockDefinition::XIndex(IndexType idx) const
{
    return idx % NumX;
}
IndexType BlockDefinition::YIndex(IndexType idx) const
{
    return (idx / NumX) % NumY;
}
IndexType BlockDefinition::ZIndex(IndexType idx) const
{
    return idx / (NumX * NumY);
}
std::tuple<IndexType, IndexType, IndexType> BlockDefinition::XYZIndices(IndexType idx) const
{
    return std::make_tuple(XIndex(idx), YIndex(idx), ZIndex(idx));
}

IndexType BlockDefinition::NumBlocks() const
{
    return NumX * NumY * NumZ;
}

IndexType BlockDefinition::OffsetIndex(IndexType idx, IndexType ox, IndexType oy, IndexType oz) const
{
    return idx + ox + oy * NumX + oz * NumX * NumY;
}

bool BlockDefinition::InDef(IndexType x, IndexType y, IndexType z) const
{
    if (x < 0 || x >= NumX ||
        y < 0 || y >= NumY ||
        z < 0 || z >= NumZ) {
        return false;
    }
    return true;
}

bool BlockDefinition::InDef(IndexType idx) const
{
    if (idx < 0 || idx >= NumBlocks()) {
        return false;
    }
    return true;
}

bool BlockDefinition::OffsetInDef(IndexType x, IndexType y, IndexType z, 
        IndexType ox, IndexType oy, IndexType oz) const
{
    return InDef(x + ox, y + oy, z + oz);
}

bool BlockDefinition::OffsetInDef(IndexType idx, IndexType ox, IndexType oy, IndexType oz) const
{
    return InDef(OffsetIndex(idx, ox, oy, oz));
}

BlockDefinition BlockDefinition::UnitModel(IndexType iNumX, IndexType iNumY, IndexType iNumZ)
{
    BlockDefinition def(iNumX, iNumY, iNumZ, 0, 0, 0, 1, 1, 1);
    return def;
}

const double PI  = 3.141592653589793238462643383280;
const double TAU = 6.283185307179586476925286766559;

////////////////////////////////////////////////////////////////////////////////

AzmSlopePair::AzmSlopePair()
{
}

AzmSlopePair::AzmSlopePair(double iazm, double islope) 
    : Azimuth(iazm), Slope(islope)
{
}

AzmSlopePair::~AzmSlopePair()
{
}

std::ostream& mvd::mineflow::operator<<(std::ostream& os, const AzmSlopePair& a)
{
    os << "{" << a.Azimuth << " (" << ToDegrees(a.Azimuth) << "_deg),  " << a.Slope << " (" << ToDegrees(a.Slope) << "_deg)}";
    return os;
}

bool AzmSlopePair::operator<(const AzmSlopePair& other) const
{
    if (Azimuth == other.Azimuth) {
        return Slope < other.Slope;
    }
    return Azimuth < other.Azimuth;
}

bool AzmSlopePair::operator<(double otherAzimuth) const
{
    return Azimuth < otherAzimuth;
}

////////////////////////////////////////////////////////////////////////////////


typedef std::vector<AzmSlopePair>::const_iterator pairIter;
typedef std::pair<pairIter, pairIter> pairIterPair;
static pairIterPair GetLeftRight(const std::vector<AzmSlopePair>& pairs,
        double azimuth)
{
    while (azimuth >= TAU) azimuth -= TAU;
    while (azimuth < 0) azimuth += TAU;

    auto right = std::lower_bound(pairs.begin(), pairs.end(), azimuth);
    if (right == pairs.end()) right = pairs.begin();
    auto left = (right == pairs.begin()) ? std::prev(pairs.end()) : std::prev(right);

    return std::make_pair(left, right);
}
static double GetXval(const pairIter& left, const pairIter& right, double azimuth)
{
    double to_left, to_right;
    if (left->Azimuth > azimuth) {
        to_left = TAU - left->Azimuth + azimuth;
    } else {
        to_left = azimuth - left->Azimuth;
    }
    if (right->Azimuth < azimuth) {
        to_right = TAU - azimuth + right->Azimuth;
    } else {
        to_right = right->Azimuth - azimuth;
    }

    return to_left / (to_left + to_right);
}

////////////////////////////////////////////////////////////////////////////////

SlopeDefinition::SlopeDefinition(std::initializer_list<std::initializer_list<double>> list)
{
    for (auto & rec : list) {
        double azm = *rec.begin();
        double slope = *std::next(rec.begin());

        while (azm >= TAU) azm -= TAU;
        while (azm < 0) azm += TAU;

        AzmSlopePair as;
        as.Azimuth = azm;
        as.Slope = slope;

        m_Pairs.push_back(as);
    }

    std::sort(m_Pairs.begin(), m_Pairs.end());
}

SlopeDefinition::SlopeDefinition(const std::vector<AzmSlopePair>& pairs) 
    : m_Pairs(pairs)
{
    std::sort(m_Pairs.begin(), m_Pairs.end());
}

std::ostream& mvd::mineflow::operator<<(std::ostream& os, const SlopeDefinition& def)
{
    for (auto & pair : def.Pairs()) {
        os << pair << std::endl;
    }
    return os;
}

SlopeDefinition SlopeDefinition::Constant(double slope)
{
    SlopeDefinition def({{0, slope}});
    return def;
}

double SlopeDefinition::operator()(double azimuth) const
{
    return Get(azimuth);
}

double SlopeDefinition::Get(double azimuth) const
{
    if (m_Pairs.empty()) {
        return 0.0;
    }
    if (m_Pairs.size() == 1) {
        return m_Pairs[0].Slope;
    }

    auto lr = GetLeftRight(m_Pairs, azimuth);
    double xval = GetXval(lr.first, lr.second, azimuth);

    double slope = lr.first->Slope + (lr.second->Slope - lr.first->Slope) * xval;
    return slope;
}

bool SlopeDefinition::Within(const Vector3D& vec) const
{
    return Within(vec.x, vec.y, vec.z);
}

bool SlopeDefinition::Within(double dx, double dy, double dz) const
{
    if (dx == 0 && dy == 0) {
        return true;
    }
    double dt = std::sqrt(dx * dx + dy * dy);
    double theta = std::atan(std::abs(dz) / dt);
    double azm = PI / 2 - std::atan2(dy, dx);
    return theta >= Get(azm);
}

double SlopeDefinition::MinSlope() const
{
    if (m_Pairs.empty()) {
        return 0.0;
    }

    double minslope = m_Pairs[0].Slope;
    for (auto & pair : m_Pairs) {
        if (pair.Slope < minslope) {
            minslope = pair.Slope;
        }
    }
    return minslope;
}

uint64_t SlopeDefinition::NumPairs() const
{
    return static_cast<uint64_t>(m_Pairs.size());
}

bool SlopeDefinition::Empty() const
{
    return m_Pairs.empty();
}

const std::vector<AzmSlopePair>& SlopeDefinition::Pairs() const
{
    return m_Pairs;
}

////////////////////////////////////////////////////////////////////////////////

SlopeDefinition mvd::mineflow::CubicInterpolate(const SlopeDefinition& def, int cnt)
{
    if (def.NumPairs() < 4) {
        throw std::runtime_error("must be at least 4 pairs for cubic");
    }

    auto y0it = std::prev(def.Pairs().end());
    auto y1it = def.Pairs().begin();
    auto y2it = std::next(y1it);
    auto y3it = std::next(y2it);

    std::vector<AzmSlopePair> outPairs;
    outPairs.resize(cnt + 1);
    int i = 0;
    for (auto v : Linspace(0, TAU, cnt + 1)) {
        if (v >= y2it->Azimuth && y2it->Azimuth != 0) {
            y0it = y1it;
            y1it = y2it;
            y2it = y3it;
            y3it = std::next(y3it);
            if (y3it == def.Pairs().end()) y3it = def.Pairs().begin();
        }

        double mu = GetXval(y1it, y2it, v);

        double y0 = y0it->Slope;
        double y1 = y1it->Slope;
        double y2 = y2it->Slope;
        double y3 = y3it->Slope;

        double mu2 = mu * mu;
        double a0 = y3 - y2 - y0 + y1;
        double a1 = y0 - y1 - a0;
        double a2 = y2 - y0;
        double a3 = y1;

        outPairs[i].Azimuth = v;

        double yn = a0 * mu * mu2 + a1 * mu2 + a2 * mu + a3;
        outPairs[i].Slope = yn;
        i++;
    }
    outPairs.pop_back();

    return SlopeDefinition(outPairs);
}

SlopeDefinition mvd::mineflow::CosineInterpolate(const SlopeDefinition& def, int cnt)
{
    if (def.NumPairs() < 2) {
        throw std::runtime_error("must be at least 2 pairs for cosine");
    }

    auto y1it = def.Pairs().begin();
    auto y2it = std::next(y1it);

    std::vector<AzmSlopePair> outPairs(cnt + 1);
    int i = 0;
    for (auto v : Linspace(0, TAU, cnt + 1)) {
        if (v >= y2it->Azimuth && y2it->Azimuth != 0) {
            y1it = y2it;
            y2it = std::next(y2it);
            if (y2it == def.Pairs().end()) y2it = def.Pairs().begin();
        }

        double mu = GetXval(y1it, y2it, v);
        double y1 = y1it->Slope;
        double y2 = y2it->Slope;

        double mu2 = (1 - std::cos(mu * PI)) / 2.0;

        outPairs[i].Azimuth = v;
        double yn = y1 * (1 - mu2) + y2 * mu2;
        outPairs[i].Slope = yn;
        i++;
    }
    outPairs.pop_back();

    return SlopeDefinition(outPairs);
}

////////////////////////////////////////////////////////////////////////////////

IBlockIndexInputIteratorSource::IBlockIndexInputIteratorSource()
{
}

IBlockIndexInputIteratorSource::~IBlockIndexInputIteratorSource()
{
}

////////////////////////////////////////////////////////////////////////////////

BlockIndexInputIteratorBase::BlockIndexInputIteratorBase(
        IBlockIndexInputIteratorSource* source)
    : m_Source(source)
{
}
BlockIndexInputIteratorBase::~BlockIndexInputIteratorBase()
{
    if (m_Source) {
        delete m_Source;
    }
}

BlockIndexInputIteratorBase::iterator BlockIndexInputIteratorBase::begin() const
{
    return BlockIndexInputIteratorBase::iterator(m_Source);
}

BlockIndexInputIteratorBase::iterator BlockIndexInputIteratorBase::end() const
{
    return BlockIndexInputIteratorBase::iterator();
}

////////////////////////////////////////////////////////////////////////////////

BlockIndexInputIteratorBase::iterator::iterator(IBlockIndexInputIteratorSource* source)
    : m_Source(source)
{
    Next();
}

BlockIndexInputIteratorBase::iterator::reference BlockIndexInputIteratorBase::iterator::operator*() const
{
    return m_CurrentIndex;
}

BlockIndexInputIteratorBase::iterator& BlockIndexInputIteratorBase::iterator::operator++()
{
    Next();
    return *this;
}

BlockIndexInputIteratorBase::iterator& BlockIndexInputIteratorBase::iterator::operator++(int)
{
    Next();
    return *this;
}

bool BlockIndexInputIteratorBase::iterator::operator==(BlockIndexInputIteratorBase::iterator rhs) const
{
    return m_Source == rhs.m_Source;
}

bool BlockIndexInputIteratorBase::iterator::operator!=(BlockIndexInputIteratorBase::iterator rhs) const
{
    return m_Source != rhs.m_Source;
}

void BlockIndexInputIteratorBase::iterator::Next()
{
    if (m_Source && m_Source->HasMore()) {
        m_CurrentIndex = m_Source->Next();
    } else {
        m_Source = nullptr;
    }
}

////////////////////////////////////////////////////////////////////////////////

IPrecedenceConstraintInputIteratorSource::IPrecedenceConstraintInputIteratorSource()
{
}

IPrecedenceConstraintInputIteratorSource::~IPrecedenceConstraintInputIteratorSource()
{
}

////////////////////////////////////////////////////////////////////////////////

PrecedenceConstraintInputIteratorBase::PrecedenceConstraintInputIteratorBase(
        IPrecedenceConstraintInputIteratorSource* source)
    : m_Source(source)
{
}

PrecedenceConstraintInputIteratorBase::~PrecedenceConstraintInputIteratorBase()
{
    if (m_Source) {
        delete m_Source;
    }
}

PrecedenceConstraintInputIteratorBase::iterator PrecedenceConstraintInputIteratorBase::begin() const
{
    return PrecedenceConstraintInputIteratorBase::iterator(m_Source);
}

PrecedenceConstraintInputIteratorBase::iterator PrecedenceConstraintInputIteratorBase::end() const
{
    return PrecedenceConstraintInputIteratorBase::iterator();
}

////////////////////////////////////////////////////////////////////////////////

PrecedenceConstraintInputIteratorBase::iterator::iterator(
        IPrecedenceConstraintInputIteratorSource* source)
    : m_Source(source)
{
    Next();
}

PrecedenceConstraintInputIteratorBase::iterator::reference PrecedenceConstraintInputIteratorBase::iterator::operator*() const
{
    return m_CurrentPrecedenceConstraint;
}

PrecedenceConstraintInputIteratorBase::iterator& PrecedenceConstraintInputIteratorBase::iterator::operator++()
{
    Next();
    return *this;
}

PrecedenceConstraintInputIteratorBase::iterator& PrecedenceConstraintInputIteratorBase::iterator::operator++(int)
{
    Next();
    return *this;
}

bool PrecedenceConstraintInputIteratorBase::iterator::operator==(PrecedenceConstraintInputIteratorBase::iterator rhs) const
{
    return m_Source == rhs.m_Source;
}

bool PrecedenceConstraintInputIteratorBase::iterator::operator!=(PrecedenceConstraintInputIteratorBase::iterator rhs) const
{
    return m_Source != rhs.m_Source;
}

void PrecedenceConstraintInputIteratorBase::iterator::Next()
{
    if (m_Source && m_Source->HasMore()) {
        m_CurrentPrecedenceConstraint = m_Source->Next();
    } else {
        m_Source = nullptr;
    }
}

////////////////////////////////////////////////////////////////////////////////

class SimplePrecedenceConstraintInputIteratorSource : public IPrecedenceConstraintInputIteratorSource
{
public:
    SimplePrecedenceConstraintInputIteratorSource(IndexType numBlocks,
            std::function<BlockIndexInputIteratorBase(IndexType)> antecedents) 
        : m_NumBlocks(numBlocks)
        , m_CurrentBlockIndex(0)
        , m_AntecedentsFunc(antecedents)
    {
        
        for (auto & v : m_AntecedentsFunc(m_CurrentBlockIndex)) {
            m_Remaining.push(v);
        }
        PopulateRemaining();
    }
    ~SimplePrecedenceConstraintInputIteratorSource() {}

    PrecedenceConstraint Next()
    {
        PrecedenceConstraint c;
        c.From = m_CurrentBlockIndex;
        c.To = m_Remaining.front();
        m_Remaining.pop();
        PopulateRemaining();
        return c;
    }

    bool HasMore() const
    {
        return !m_Remaining.empty();
    }

    void PopulateRemaining()
    {
        while (m_Remaining.empty() && m_CurrentBlockIndex < (m_NumBlocks - 1)) {
            m_CurrentBlockIndex++;
            for (auto & v : m_AntecedentsFunc(m_CurrentBlockIndex)) {
                m_Remaining.push(v);
            }
        }
    }

private:
    IndexType m_NumBlocks;
    IndexType m_CurrentBlockIndex;
    std::function<BlockIndexInputIteratorBase(IndexType)> m_AntecedentsFunc;
    std::queue<IndexType> m_Remaining;
};


////////////////////////////////////////////////////////////////////////////////

IPrecedenceConstraints::IPrecedenceConstraints()
{
}

IPrecedenceConstraints::~IPrecedenceConstraints()
{
}

BlockIndexInputIteratorBase IPrecedenceConstraints::Successors(IndexType toBlockIndex) const
{
    throw std::runtime_error("not implemented");
}

IndexType IPrecedenceConstraints::NumAntecedents(IndexType fromBlockIndex) const
{
    IndexType cnt = 0;
    for (auto & v : Antecedents(fromBlockIndex)) {
        cnt++;
    }
    return cnt;
}

IndexType IPrecedenceConstraints::ApproxNumAntecedents(IndexType fromBlockIndex) const
{
    return 0;
}

void IPrecedenceConstraints::AntecedentsVector(IndexType fromBlockIndex, std::vector<IndexType>* vec) const
{
    MVD_ASSERT(vec);
    vec->clear();
    for (auto & v : Antecedents(fromBlockIndex)) {
        vec->push_back(v);
    }
}

IndexType IPrecedenceConstraints::NumSuccessors(IndexType toBlockIndex) const
{
    IndexType cnt = 0;
    for (auto & v : Successors(toBlockIndex)) {
        cnt++;
    }
    return cnt;
}

IndexType IPrecedenceConstraints::ApproxNumSuccessors(IndexType toBlockIndex) const
{
    return 0;
}

void IPrecedenceConstraints::SuccessorsVector(IndexType toBlockIndex, std::vector<IndexType>* vec) const
{
    MVD_ASSERT(vec);
    vec->clear();
    for (auto & v : Successors(toBlockIndex)) {
        vec->push_back(v);
    }
}

// Sometimes we want all the precedence constraints
IndexType IPrecedenceConstraints::NumPrecedenceConstraints() const
{
    IndexType cnt = 0;
    for (IndexType blockIndex = 0; blockIndex < NumBlocks(); blockIndex++) {
        cnt += NumAntecedents(blockIndex);
    }
    return cnt;
}

IndexType IPrecedenceConstraints::ApproxNumPrecedenceConstraints() const
{
    return 0;
}

PrecedenceConstraintInputIteratorBase IPrecedenceConstraints::PrecedenceConstraints() const
{
    return PrecedenceConstraintInputIteratorBase(new 
        SimplePrecedenceConstraintInputIteratorSource(
            NumBlocks(),
            std::bind(&IPrecedenceConstraints::Antecedents, this, std::placeholders::_1)
    ));
}

void IPrecedenceConstraints::PrecedenceConstraintsVector(std::vector<PrecedenceConstraint>* vec) const
{
    MVD_ASSERT(vec);
    vec->clear();
    for (auto & v : PrecedenceConstraints()) {
        vec->push_back(v);
    }
}


PrecedenceConstraintsReachableSearchBufferPtr IPrecedenceConstraints::GetNewSearchBuffer() const
{
    return std::make_unique<PrecedenceConstraintsReachableSearchBuffer>(NumBlocks());
}

class ReachableBlockSource : public IBlockIndexInputIteratorSource
{
public:
    ReachableBlockSource(
            IndexType blockIndex,
            std::function<BlockIndexInputIteratorBase(IndexType)> func,
            PrecedenceConstraintsReachableSearchBuffer* buffer)
        : m_Func(func)
        , m_Buffer(buffer)
    {
        MVD_ASSERT(m_Buffer);
        m_Buffer->NewSearch();
        for (auto & v : m_Func(blockIndex)) {
            m_Buffer->Queue(v);
        }
    }
    ~ReachableBlockSource() {}

    IndexType Next() override final
    {
        IndexType v;
        m_Buffer->Search(&v);
        for (auto & t : m_Func(v)) {
            m_Buffer->Queue(t);
        }
        return v;
    }

    bool HasMore() const override final
    {
        return m_Buffer->HasMore();
    }

private:
    std::function<BlockIndexInputIteratorBase(IndexType)> m_Func;
    PrecedenceConstraintsReachableSearchBuffer* m_Buffer;
};

BlockIndexInputIteratorBase IPrecedenceConstraints::ReachableAntecedents(IndexType fromBlockIndex, 
        PrecedenceConstraintsReachableSearchBuffer* buffer) const
{
    return BlockIndexInputIteratorBase(new ReachableBlockSource(
        fromBlockIndex,
        std::bind(&IPrecedenceConstraints::Antecedents, this, std::placeholders::_1),
        buffer)
    );
}

BlockIndexInputIteratorBase IPrecedenceConstraints::ReachableSuccessors(IndexType toBlockIndex, 
        PrecedenceConstraintsReachableSearchBuffer* buffer) const
{
    return BlockIndexInputIteratorBase(new ReachableBlockSource(
        toBlockIndex,
        std::bind(&IPrecedenceConstraints::Successors, this, std::placeholders::_1),
        buffer)
    );
}

static void PartialSearch(IndexType start,
        std::function<bool(IndexType v)> cback,
        std::function<BlockIndexInputIteratorBase(IndexType v)> func,
        PrecedenceConstraintsReachableSearchBuffer* buffer)
{
    buffer->NewSearch();
    for (auto & to : func(start)) {
        buffer->Queue(to);
    }

    IndexType v;
    while (buffer->Search(&v)) {
        if (cback(v)) {
            for (auto & to : func(v)) {
                buffer->Queue(to);
            }
        }
    }
}

void IPrecedenceConstraints::PartialReachableAntecedents(IndexType fromBlockIndex,
        std::function<bool(IndexType toBlockIndex)> cback,
        PrecedenceConstraintsReachableSearchBuffer* buffer) const
{
    PartialSearch(fromBlockIndex, cback,
            std::bind(&IPrecedenceConstraints::Antecedents, this, std::placeholders::_1),
            buffer);
}

void IPrecedenceConstraints::PartialReachableSuccessors(IndexType toBlockIndex,
        std::function<bool(IndexType fromBlockIndex)> cback,
        PrecedenceConstraintsReachableSearchBuffer* buffer) const
{
    PartialSearch(toBlockIndex, cback,
            std::bind(&IPrecedenceConstraints::Successors, this, std::placeholders::_1),
            buffer);
}

////////////////////////////////////////////////////////////////////////////////

PrecedenceConstraintsReachableSearchBuffer::PrecedenceConstraintsReachableSearchBuffer(
        IndexType numBlocks) 
    : m_NumBlocks(numBlocks)
{
    m_Tag = 101;
}

PrecedenceConstraintsReachableSearchBuffer::~PrecedenceConstraintsReachableSearchBuffer() {
}

void PrecedenceConstraintsReachableSearchBuffer::NewSearch()
{
    if (m_Tag >= 100) {
        m_Seen.assign(m_NumBlocks, 101);
        m_Tag = 0;
    } else {
        m_Tag++;
    }
    while (!m_Queue.empty()) {
        m_Queue.pop();
    }
}

void PrecedenceConstraintsReachableSearchBuffer::Queue(IndexType v)
{
    if (m_Seen[v] != m_Tag) {
        m_Seen[v] = m_Tag;
        m_Queue.push(v);
    }
}

bool PrecedenceConstraintsReachableSearchBuffer::Search(IndexType* v)
{
    if (m_Queue.empty()) {
        return false;
    }
    *v = m_Queue.front();
    m_Queue.pop();
    return true;
}

bool PrecedenceConstraintsReachableSearchBuffer::HasMore()
{
    return !m_Queue.empty();
}

////////////////////////////////////////////////////////////////////////////////

bool mvd::mineflow::ConsistentPrecedenceConstraints(IPrecedenceConstraints* pre)
{
    MVD_ASSERT(pre);

    IndexType preNumBlocks = pre->NumBlocks();
    IndexType preNumPrecedenceConstraints = pre->NumPrecedenceConstraints();

    std::unordered_map<IndexType, std::unordered_set<IndexType>> antecedents;
    std::unordered_map<IndexType, std::unordered_set<IndexType>> successors;

    std::unordered_map<IndexType, std::unordered_set<IndexType>> mySuccessors;

    for (IndexType blockIndex = 0; blockIndex < preNumBlocks; blockIndex++) {
        IndexType nAnte = 0;
        for (auto & targetBlockIndex : pre->Antecedents(blockIndex)) {
            antecedents[blockIndex].insert(targetBlockIndex);
            mySuccessors[targetBlockIndex].insert(blockIndex);
            nAnte++;
        }
        if (nAnte != antecedents[blockIndex].size()) {
            return false;
        }
        IndexType nSucc = 0;
        for (auto & targetBlockIndex : pre->Successors(blockIndex)) {
            successors[blockIndex].insert(targetBlockIndex);
            nSucc++;
        }
        if (nSucc != successors[blockIndex].size()) {
            return false;
        }
    }

    for (IndexType blockIndex = 0; blockIndex < preNumBlocks; blockIndex++) {
        if (successors[blockIndex] != mySuccessors[blockIndex]) {
            return false;
        }
    }

    IndexType actualNumber = 0;
    for (IndexType blockIndex = 0; blockIndex < preNumBlocks; blockIndex++) {
        actualNumber += static_cast<IndexType>(antecedents[blockIndex].size());
    }
    if (actualNumber != preNumPrecedenceConstraints) {
        return false;
    }

    // todo could check for cycles..

    return true;
}

////////////////////////////////////////////////////////////////////////////////

class BlockOffsetSource : public IBlockIndexInputIteratorSource
{
public:
    BlockOffsetSource(IndexType blockIndex, const IndexType* data, IndexType n);
    ~BlockOffsetSource();
    IndexType Next() override final;
    bool HasMore() const override final;

private:
    IndexType m_BlockIndex;
    const IndexType* m_Offsets;
    IndexType m_NumOffsets;
};


BlockOffsetSource::BlockOffsetSource(IndexType blockIndex, const IndexType* data, IndexType n)
    : m_BlockIndex(blockIndex)
    , m_Offsets(data)
    , m_NumOffsets(n)
{
}

BlockOffsetSource::~BlockOffsetSource()
{
}

IndexType BlockOffsetSource::Next()
{
    MVD_ASSERT(m_Offsets);
    MVD_ASSERT(m_NumOffsets > 0);
    IndexType ret = m_BlockIndex + *m_Offsets;
    m_Offsets++;
    m_NumOffsets--;
    if (m_NumOffsets == 0) {
        m_Offsets = nullptr;
    }
    return ret;
}

bool BlockOffsetSource::HasMore() const
{
    return m_NumOffsets > 0;
}

////////////////////////////////////////////////////////////////////////////////

Regular2DGrid45DegreePrecedence::Regular2DGrid45DegreePrecedence(IndexType numX, IndexType numZ) 
    : m_NumX(numX)
    , m_NumZ(numZ)
{
    if (m_NumX <= 0 || m_NumZ <= 0) {
        throw std::invalid_argument("invalid grid size");
    }

    m_AntecedentOffsets = {
        static_cast<IndexType>(m_NumX - 1),
        static_cast<IndexType>(m_NumX),
        static_cast<IndexType>(m_NumX + 1),
    };
    m_SuccessorOffsets = {
        static_cast<IndexType>(-m_NumX - 1),
        static_cast<IndexType>(-m_NumX),
        static_cast<IndexType>(-m_NumX + 1),
    };
}

Regular2DGrid45DegreePrecedence::~Regular2DGrid45DegreePrecedence()
{
}

IndexType Regular2DGrid45DegreePrecedence::NumBlocks() const
{
    return m_NumX * m_NumZ;
}

BlockIndexInputIteratorBase Regular2DGrid45DegreePrecedence::XAdjustedSource(
        IndexType blockIndex, const std::vector<IndexType>& offsets) const
{
    const IndexType* data = offsets.data();
    IndexType n = 3;
    IndexType x = blockIndex % m_NumX;

    if (x == 0) {
        data++;
        n--;
    }
    if (x == (m_NumX - 1)) {
        n--;
    }

    return BlockIndexInputIteratorBase(new BlockOffsetSource(blockIndex, data, n));
}

BlockIndexInputIteratorBase Regular2DGrid45DegreePrecedence::Antecedents(IndexType fromBlockIndex) const
{
    IndexType z = fromBlockIndex / m_NumX;
    if (z >= (m_NumZ - 1)) {
        return BlockIndexInputIteratorBase(nullptr);
    }

    return XAdjustedSource(fromBlockIndex, m_AntecedentOffsets);

}

BlockIndexInputIteratorBase Regular2DGrid45DegreePrecedence::Successors(IndexType toBlockIndex) const
{
    IndexType z = toBlockIndex / m_NumX;
    if (z <= 0) {
        return BlockIndexInputIteratorBase(nullptr);
    }

    return XAdjustedSource(toBlockIndex, m_SuccessorOffsets);
}

////////////////////////////////////////////////////////////////////////////////

PrecedencePattern::PrecedencePattern()
{
}

PrecedencePattern::~PrecedencePattern()
{
}

uint64_t PrecedencePattern::size() const
{
    return static_cast<uint64_t>(Offsets.size());
}

PrecedencePattern::iterator PrecedencePattern::begin()
{
    return Offsets.begin();
}

PrecedencePattern::iterator PrecedencePattern::end()
{
    return Offsets.end();
}

PrecedencePattern::const_iterator PrecedencePattern::begin() const
{
    return Offsets.begin();
}

PrecedencePattern::const_iterator PrecedencePattern::end() const
{
    return Offsets.end();
}

PrecedencePattern PrecedencePattern::OneFive()
{
    PrecedencePattern ptrn;
    ptrn.Offsets = std::vector<Vector3IT>({
        {0, -1, 1}, {-1, 0, 1}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}
    });
    return ptrn;
}

PrecedencePattern PrecedencePattern::OneNine()
{
    PrecedencePattern ptrn;
    ptrn.Offsets.resize(9);

    IndexType k = 0;
    for (IndexType j = -1; j <= 1; j++) {
        for (IndexType i = -1; i <= 1; i++) {
            ptrn.Offsets[k].x = i;
            ptrn.Offsets[k].y = j;
            ptrn.Offsets[k].z = 1;
            k++;
        }
    }
    return ptrn;
}

PrecedencePattern PrecedencePattern::KnightsMove()
{
    PrecedencePattern ptrn;
    ptrn.Offsets = std::vector<Vector3IT>({
        { 0, -1, 1}, {-1,  0, 1}, { 0,  0,  1}, { 1,  0,  1}, { 0,  1,  1},
        {-1, -2, 2}, { 1, -2, 2}, {-2, -1,  2}, { 2, -1,  2},
        {-2,  1, 2}, { 2,  1, 2}, {-1,  2,  2}, { 1,  2,  2}
    });
    return ptrn;
}

PrecedencePattern PrecedencePattern::Naive(const BlockDefinition& block_def,
                       const SlopeDefinition& slope_def,
                       IndexType numZ)
{
    PrecedencePattern ptrn;
    NaiveSearch(block_def, slope_def, numZ, [&](Vector3IT off){
        ptrn.Offsets.push_back(off);
    });
    return ptrn;
}

PrecedencePattern PrecedencePattern::LessNaive(
        const BlockDefinition& block_def,
        const SlopeDefinition& slope_def,
        IndexType numZ)
{
    // Keeps track of x / y and doesn't make things at the same x / y location (instead
    // relying on the original x/y
    PrecedencePattern ptrn;

    std::unordered_map<IndexType, std::unordered_set<IndexType>> seen;

    NaiveSearch(block_def, slope_def, numZ, [&](Vector3IT off){
        auto itx = seen.find(off.x);
        bool add = false;
        if (itx == seen.end()) {
            add = true;
        } else {
            auto ity = itx->second.find(off.y);
            if (ity == itx->second.end())
            {
                add = true;
            }
        }
        if (add) {
            ptrn.Offsets.push_back(off);
            seen[off.x].insert(off.y);
        }
    });

    return ptrn;
}

PrecedencePattern PrecedencePattern::MinSearch(
        const BlockDefinition& block_def,
        const SlopeDefinition& slope_def,
        IndexType n)
{
    PrecedencePattern ptrn;
    if (slope_def.Empty()) {
        return ptrn;
    }

    double min_slope = slope_def.MinSlope();
    double max_height = block_def.SizeZ * n;
    double max_throw = max_height / std::tan(min_slope);

    IndexType cx = static_cast<IndexType>(std::ceil(max_throw));
    IndexType nx = cx * 2 + 1;
    IndexType nz = n + 1;

    IndexType total = nx * nx * nz;

    BlockDefinition flag_def = block_def;
    flag_def.NumX = nx;
    flag_def.NumY = nx;
    flag_def.NumZ = nz;

    const char NOT_FLAGGED = 0;
    const char NO_ARCS = 1;
    const char SOME_ARCS = 2;
    std::vector<char> flag(total, NOT_FLAGGED);

    // Now construct the minimum search pattern,
    IndexType origin_index = flag_def.GridIndex(cx, cx, 0);
    flag[origin_index] = NO_ARCS;
    std::vector<IndexType> flagged;
    flagged.push_back(origin_index);
    for (IndexType z = 1; z < nz; z++) {
        double this_height = z * flag_def.SizeZ;
        double this_throw = this_height / std::tan(min_slope);
        IndexType this_max_off = static_cast<IndexType>(std::ceil(this_throw));

        // flag the violating blocks
        std::vector<Vector3IT> new_offsets;
        for (IndexType x = -this_max_off; x <= this_max_off; x++) {
            for (IndexType y = -this_max_off; y <= this_max_off; y++) {
                IndexType fi = flag_def.OffsetIndex(origin_index, x, y, z);

                if (flag[fi] == NOT_FLAGGED && slope_def.Within(x * flag_def.SizeX, 
                                                                y * flag_def.SizeY, 
                                                                z * flag_def.SizeZ)) {
                    flag[fi] = NO_ARCS;

                    flagged.push_back(fi);
                    ptrn.Offsets.emplace_back(x, y, z);
                    new_offsets.emplace_back(x, y, z);
                }
            }
        }

        // Extend flagged blocks
        std::vector<IndexType> extra;
        for (auto & fi : flagged) {
            IndexType fz = flag_def.ZIndex(fi);

            std::vector<Vector3IT>* offsets;
            if (flag[fi] == NO_ARCS) {
                offsets = &ptrn.Offsets;
                flag[fi] = SOME_ARCS;
            } else {
                offsets = &new_offsets;
            }

            for (auto & arc : (*offsets)) {
                if (fz + arc.z >= flag_def.NumZ) {
                    break;
                }

                IndexType idx = flag_def.OffsetIndex(fi, arc.x, arc.y, arc.z);
                if (flag[idx] == NOT_FLAGGED) {
                    flag[idx] = NO_ARCS;
                    extra.push_back(idx);
                }
            }
        }

        if (!extra.empty()) {
            flagged.reserve(flagged.size() + extra.size());
            flagged.insert(flagged.end(), extra.begin(), extra.end());
        }
    }

    return ptrn;
}

PrecedencePattern PrecedencePattern::MinSearch(double slope_rad, IndexType nz)
{
    BlockDefinition block_def = BlockDefinition::UnitModel(1, 1, 1);
    SlopeDefinition slope_def = SlopeDefinition::Constant(slope_rad);
    return MinSearch(block_def, slope_def, nz);
}

////////////////////////////////////////////////////////////////////////////////

void mvd::mineflow::NaiveSearch(const BlockDefinition& block_def,
                        const SlopeDefinition& slope_def,
                        IndexType nz, std::function<void(Vector3IT)> offsetCallback)
{
    if (slope_def.Empty()) {
        return;
    }

    double min_slope = slope_def.MinSlope();
    if (min_slope <= 0) {
        return;
    }

    for (IndexType z = 1; z <= nz; z++) {
        double this_height = block_def.SizeZ * z;
        double this_throw = this_height / std::tan(min_slope);
        IndexType this_max_off = static_cast<IndexType>(std::ceil(this_throw));

        for (IndexType x = -this_max_off; x <= this_max_off; x++) {
            for (IndexType y = -this_max_off; y <= this_max_off; y++) {
                if (slope_def.Within(x * block_def.SizeX, 
                                     y * block_def.SizeY, 
                                     z * block_def.SizeZ)) {
                    offsetCallback({x, y, z});
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

#define TRUE_NEGATIVE 0
#define FALSE_POSITIVE 1
#define TRUE_POSITIVE 2
#define FALSE_NEGATIVE 3

static void GetAccuracyFlag(const BlockDefinition& block_def,
                            const SlopeDefinition& slope_def,
                            const PrecedencePattern& ptrn,
                            std::vector<char>* flagp)
{
    auto& flag = *flagp;
    flag.assign(block_def.NumBlocks(), TRUE_NEGATIVE);

    uint64_t mx = block_def.NumX / 2;
    uint64_t my = block_def.NumY / 2;

    uint64_t start = block_def.GridIndex(mx, my, 0);
    NaiveSearch(block_def, slope_def, block_def.NumZ, [&](Vector3IT off){
        if (block_def.OffsetInDef(mx, my, 0, off.x, off.y, off.z)) {
            uint64_t idx = block_def.OffsetIndex(start, off.x, off.y, off.z);
            flag[idx] = FALSE_POSITIVE;
        }
    });


    // Now apply the arc template, the flag will keep track of duplicates, and
    // avoid being inefficient
    std::vector<uint64_t> stack {start};
    stack.reserve(ptrn.Offsets.size());
    while (!stack.empty()) {
        uint64_t t = stack.back();
        stack.pop_back();

        for (auto & offset : ptrn.Offsets) {
            if (block_def.OffsetInDef(t, offset.x, offset.y, offset.z)) {
                IndexType idx = block_def.OffsetIndex(t, offset.x, offset.y, offset.z);

                if (flag[idx] > 1) {
                    continue;
                }

                if (flag[idx] == FALSE_POSITIVE) {
                    flag[idx] = TRUE_POSITIVE;
                } else if (flag[idx] == TRUE_NEGATIVE) {
                    flag[idx] = FALSE_NEGATIVE;
                }

                stack.push_back(idx);
            }
        }
    }
}

static void ResetAccuracyCount(PatternAccuracy* accuracy)
{
    accuracy->true_positive = 0;
    accuracy->true_negative = 0;
    accuracy->false_positive = 0;
    accuracy->false_negative = 0;
}

static void CalcAccuracyMeasure(PatternAccuracy* accuracy)
{
    double tp = static_cast<double>(accuracy->true_positive);
    double fp = static_cast<double>(accuracy->false_positive);
    double tn = static_cast<double>(accuracy->true_negative);
    double fn = static_cast<double>(accuracy->false_negative);

    accuracy->accuracy = (tp + tn) / (tp + fp + tn + fn);
    accuracy->true_positive_rate = tp / (tp + fp);
    accuracy->false_negative_rate = fn / (tp + fp);

    double numer = tp * tn - fp * fn;
    double denom = std::sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
    if (denom == 0) {
        denom = 1.0;
    }
    accuracy->matthews_correlation = numer / denom;
}

void mvd::mineflow::MeasureAccuracy(const BlockDefinition& block_def,
                            const SlopeDefinition& slope_def,
                            const PrecedencePattern& ptrn, 
                            PatternAccuracy* accuracy)
{
    if (!accuracy) return;
    ResetAccuracyCount(accuracy);

    std::vector<char> flag;
    GetAccuracyFlag(block_def, slope_def, ptrn, &flag);

    for (auto & v : flag)
    {
        switch (v) {
            case TRUE_NEGATIVE: accuracy->true_negative++; break;
            case TRUE_POSITIVE: accuracy->true_positive++; break;
            case FALSE_NEGATIVE: accuracy->false_negative++; break;
            case FALSE_POSITIVE: accuracy->false_positive++; break;
        }
    }

    CalcAccuracyMeasure(accuracy);
}

void mvd::mineflow::MultiMeasureAccuracy(const BlockDefinition& block_def,
        const SlopeDefinition& slope_def,
        const PrecedencePattern& ptrn,
        std::vector<PatternAccuracy>* accuracies)
{
    if (!accuracies) return;

    std::vector<char> flag;
    GetAccuracyFlag(block_def, slope_def, ptrn, &flag);

    IndexType nz = block_def.NumZ;
    IndexType nxy = block_def.NumX * block_def.NumY;

    accuracies->resize(nz);

    accuracies->at(0).true_positive = 1;
    accuracies->at(0).true_negative = nxy - 1;
    accuracies->at(0).false_positive = 0;
    accuracies->at(0).false_negative = 0;

    IndexType k = nxy;
    for (IndexType z = 1; z < block_def.NumZ; z++) {
        accuracies->at(z) = accuracies->at(z - 1);
        for (IndexType yx = 0; yx < block_def.NumX * block_def.NumY; yx++) {
            char v = flag[k++];

            switch (v) {
                case TRUE_NEGATIVE: accuracies->at(z).true_negative++; break;
                case TRUE_POSITIVE: accuracies->at(z).true_positive++; break;
                case FALSE_NEGATIVE: accuracies->at(z).false_negative++; break;
                case FALSE_POSITIVE: accuracies->at(z).false_positive++; break;
            }
        }
    }

    for (auto & acc : *accuracies) {
        CalcAccuracyMeasure(&acc);
    }
}

std::ostream& mvd::mineflow::operator<<(std::ostream& os, const PatternAccuracy& acc)
{
    os << "tp " << acc.true_positive << std::endl;
    os << "tn " << acc.true_negative << std::endl;
    os << "fp " << acc.false_positive << std::endl;
    os << "fn " << acc.false_negative << std::endl;
    os << "ac " << acc.accuracy << std::endl;
    os << "tr " << acc.true_positive_rate << std::endl;
    os << "fr " << acc.false_negative_rate << std::endl;
    os << "mc " << acc.matthews_correlation << std::endl;
    return os;
}

void mvd::mineflow::PrintPattern(const PrecedencePattern& ptrn)
{
    IndexType xlo = std::numeric_limits<IndexType>::max();
    IndexType ylo = std::numeric_limits<IndexType>::max();
    IndexType xhi = std::numeric_limits<IndexType>::min();
    IndexType yhi = std::numeric_limits<IndexType>::min();

    for (auto & off : ptrn.Offsets) {
        if (off.x < xlo) xlo = off.x;
        if (off.y < ylo) ylo = off.y;
        if (off.x > xhi) xhi = off.x;
        if (off.y > yhi) yhi = off.y;
    }
    IndexType pnx = xhi - xlo + 1;
    IndexType pny = yhi - ylo + 1;
    std::vector<IndexType> img(pnx * pny, -1);


    for (auto & off : ptrn.Offsets) {
        IndexType i = (off.y - ylo) * pnx + (off.x - xlo);
        if (img[i] == -1) {
            img[i] = off.z;
        }
        if (off.x == 0 && off.y == 0) {
            img[i] = 0;
        }
    }

    IndexType i = 0;
    for (IndexType y = 0; y < pny; y++) {
        for (IndexType x = 0; x < pnx; x++) {
            if (img[i] == -1) {
                std::cout << "  ";
            } else {
                std::cout << std::setw(2) << img[i];
            }
            i++;
        }
        std::cout << std::endl;
    }
}

////////////////////////////////////////////////////////////////////////////////

class BlockOffsetExtentSource : public IBlockIndexInputIteratorSource
{
public:
    BlockOffsetExtentSource(
            IndexType x, IndexType y, IndexType z, IndexType nx, IndexType ny, const Vector3IT* offsets, IndexType n) 
        : m_X(x)
        , m_Y(y)
        , m_Z(z)
        , m_NX(nx)
        , m_NY(ny)
        , m_Offsets(offsets)
        , m_NumOffsets(n)
    {
        while (m_NumOffsets > 0) {
            IndexType tx = m_X + m_Offsets->x;
            IndexType ty = m_Y + m_Offsets->y;
            IndexType tz = m_Z + m_Offsets->z;
            if (tx < 0 || tx >= m_NX ||
                ty < 0 || ty >= m_NY) {
                m_NumOffsets--;
                m_Offsets++;
            } else {
                break;
            }
        }
    }

    ~BlockOffsetExtentSource() {};

    IndexType Next() override final
    {
        IndexType tx = m_X + m_Offsets->x;
        IndexType ty = m_Y + m_Offsets->y;
        IndexType tz = m_Z + m_Offsets->z;

        m_NumOffsets--;
        m_Offsets++;
        while (m_NumOffsets > 0) {
            IndexType tx = m_X + m_Offsets->x;
            IndexType ty = m_Y + m_Offsets->y;
            IndexType tz = m_Z + m_Offsets->z;
            if (tx < 0 || tx >= m_NX ||
                ty < 0 || ty >= m_NY) {
                m_NumOffsets--;
                m_Offsets++;
            } else {
                break;
            }
        }
        return tx + ty * m_NX + tz * m_NX * m_NY;
    }

    bool HasMore() const override final
    {
        return m_NumOffsets > 0;
    }

private:
    IndexType m_X, m_Y, m_Z;
    IndexType m_NX, m_NY;

    const Vector3IT* m_Offsets;
    IndexType m_NumOffsets;
};

////////////////////////////////////////////////////////////////////////////////

Regular3DBlockModelPatternPrecedence::Regular3DBlockModelPatternPrecedence(
        const BlockDefinition& blockDef,
        const PrecedencePattern& pattern) 
    : Regular3DBlockModelPatternPrecedence(blockDef.NumX, blockDef.NumY, blockDef.NumZ, pattern)
{
}

Regular3DBlockModelPatternPrecedence::Regular3DBlockModelPatternPrecedence(
        IndexType numX, IndexType numY, IndexType numZ,
        const PrecedencePattern& pattern) 
    : m_NumX(numX)
    , m_NumY(numY)
    , m_NumZ(numZ)
    , m_Offsets(pattern.Offsets)
{
    if (m_Offsets.empty()) throw std::runtime_error("invalid pattern");

    std::sort(m_Offsets.begin(), m_Offsets.end(), [&](const Vector3IT& a, const Vector3IT& b){
        if (a.z == b.z) {
            if (a.y == b.y) {
                return a.x < b.x;
            }
            return a.y < b.y;
        }
        return a.z < b.z;
    });
    m_Precomputed1DOffsets.resize(m_Offsets.size());
    for (size_t i = 0; i < m_Offsets.size(); i++) {
        Vector3IT& off = m_Offsets[i];
        m_Precomputed1DOffsets[i] = off.x + off.y * m_NumX + off.z * m_NumX * m_NumY;
    }

    m_InnerRegion.xLo = 0;
    m_InnerRegion.xHi = m_NumX;
    m_InnerRegion.yLo = 0;
    m_InnerRegion.yHi = m_NumY;
    for (auto & off : m_Offsets) {
        if (off.x < 0) {
            IndexType xlo = -off.x;
            if (xlo > m_InnerRegion.xLo) {
                m_InnerRegion.xLo = xlo;
            }
        }
        if (off.x > 0) {
            IndexType xhi = m_NumX - off.x;
            if (xhi < m_InnerRegion.xHi) {
                m_InnerRegion.xHi = xhi;
            }
        }
        if (off.y < 0) {
            IndexType ylo = -off.y;
            if (ylo > m_InnerRegion.yLo) {
                m_InnerRegion.yLo = ylo;
            }
        }
        if (off.y > 0) {
            IndexType yhi = m_NumY - off.y;
            if (yhi < m_InnerRegion.yHi) {
                m_InnerRegion.yHi = yhi;
            }
        }
    }

    m_MaxOffsetZ = m_Offsets.back().z;
    m_NumOffsetsByZMinus.resize(m_MaxOffsetZ + 1);
    for (auto & off : m_Offsets) {
        m_NumOffsetsByZMinus[off.z]++;
    }
    std::partial_sum(m_NumOffsetsByZMinus.begin(), m_NumOffsetsByZMinus.end(), m_NumOffsetsByZMinus.begin());
}

Regular3DBlockModelPatternPrecedence::~Regular3DBlockModelPatternPrecedence()
{
}

IndexType Regular3DBlockModelPatternPrecedence::NumBlocks() const
{
    return m_NumX * m_NumY * m_NumZ;
}

BlockIndexInputIteratorBase Regular3DBlockModelPatternPrecedence::Antecedents(
        IndexType fromBlockIndex) const
{
    auto [x, y, z] = XYZ(fromBlockIndex);

    if (z == m_NumZ - 1) {
        return BlockIndexInputIteratorBase(nullptr);
    }

    IndexType zMinus = m_NumZ - z  - 1;
    IndexType n = m_Offsets.size();
    if (zMinus <= m_MaxOffsetZ) {
        n = m_NumOffsetsByZMinus[zMinus];
    }

    if (x >= m_InnerRegion.xLo && x < m_InnerRegion.xHi &&
        y >= m_InnerRegion.yLo && y < m_InnerRegion.yHi) {
        return BlockIndexInputIteratorBase(new BlockOffsetSource(fromBlockIndex, m_Precomputed1DOffsets.data(), n));
    } else {
        return BlockIndexInputIteratorBase(new BlockOffsetExtentSource(
            x, y, z, m_NumX, m_NumY, m_Offsets.data(), n));
    }
    return BlockIndexInputIteratorBase(nullptr);
}

BlockIndexInputIteratorBase Regular3DBlockModelPatternPrecedence::Successors(
        IndexType toBlockIndex) const
{
    return BlockIndexInputIteratorBase(nullptr);
}

IndexType Regular3DBlockModelPatternPrecedence::ApproxNumAntecedents(IndexType fromBlockIndex) const
{
    return m_Offsets.size();
}


std::tuple<IndexType, IndexType, IndexType> Regular3DBlockModelPatternPrecedence::XYZ(IndexType k) const
{
    return std::make_tuple(
        k % m_NumX, 
        (k / m_NumX) % m_NumY, 
        k / (m_NumX * m_NumY)
    );
}

////////////////////////////////////////////////////////////////////////////////

Regular3DBlockModelKeyedPatternsPrecedence::Regular3DBlockModelKeyedPatternsPrecedence(const BlockDefinition& blockDef,
        const std::vector<PrecedencePattern>& patterns,
        std::shared_ptr<std::vector<IndexType>> patternIndices)
    : m_PatternIndices(patternIndices)
{
    if (patterns.empty()) {
        throw std::invalid_argument("a non zero number of patterns are required");
    }
    if (patternIndices->size() != blockDef.NumBlocks()) {
        throw std::invalid_argument("invalid pattern indices count");
    }

    m_Patterns.reserve(patterns.size());
    for (size_t i = 0; i < patterns.size(); i++) {
        m_Patterns.emplace_back(blockDef, patterns[i]);
    }
}
Regular3DBlockModelKeyedPatternsPrecedence::~Regular3DBlockModelKeyedPatternsPrecedence()
{
}

IndexType Regular3DBlockModelKeyedPatternsPrecedence::NumBlocks() const 
{
    return m_Patterns.front().NumBlocks();
}

BlockIndexInputIteratorBase Regular3DBlockModelKeyedPatternsPrecedence::Antecedents(IndexType fromBlockIndex) const 
{
    return m_Patterns.at(m_PatternIndices->at(fromBlockIndex)).Antecedents(fromBlockIndex);
}

BlockIndexInputIteratorBase Regular3DBlockModelKeyedPatternsPrecedence::Successors(IndexType toBlockIndex) const 
{
    return m_Patterns.at(m_PatternIndices->at(toBlockIndex)).Successors(toBlockIndex);
}

IndexType Regular3DBlockModelKeyedPatternsPrecedence::ApproxNumAntecedents(IndexType fromBlockIndex) const 
{
    return m_Patterns.at(m_PatternIndices->at(fromBlockIndex)).ApproxNumAntecedents(fromBlockIndex);
}

////////////////////////////////////////////////////////////////////////////////

class BlockVectorSource : public IBlockIndexInputIteratorSource
{
public:
    BlockVectorSource(const std::vector<IndexType>& vec);
    ~BlockVectorSource();

    IndexType Next() override final;
    bool HasMore() const override final;
private:
    const IndexType* m_Ptr;
    size_t m_Remaining;
};

BlockVectorSource::BlockVectorSource(const std::vector<IndexType>& vec)
    : m_Ptr(vec.data())
    , m_Remaining(vec.size())
{
}

BlockVectorSource::~BlockVectorSource()
{
}

bool BlockVectorSource::HasMore() const
{
    return m_Remaining > 0;
}

IndexType BlockVectorSource::Next()
{
    IndexType v = *m_Ptr;
    m_Remaining--;
    m_Ptr++;
    return v;
}

ExplicitPrecedence::ExplicitPrecedence(IndexType numBlocks)
    : m_NumBlocks(numBlocks)
{
}

ExplicitPrecedence::ExplicitPrecedence(IndexType numBlocks,
        std::initializer_list<std::initializer_list<int>> list)
    : m_NumBlocks(numBlocks)
{
    for (auto & pair : list) {
        if (pair.size() != 2) {
            throw std::invalid_argument("invalid pair size");
        }
        m_Antecedents[*pair.begin()].push_back(*(pair.begin() + 1));
    }
}

ExplicitPrecedence::ExplicitPrecedence(IndexType numBlocks,
            std::unordered_map<IndexType, std::vector<IndexType>>&& antecedents)
    : m_NumBlocks(numBlocks)
    , m_Antecedents(antecedents)
{
}

ExplicitPrecedence::~ExplicitPrecedence()
{
}

IndexType ExplicitPrecedence::NumBlocks() const
{
    return m_NumBlocks;
}

BlockIndexInputIteratorBase ExplicitPrecedence::Antecedents(IndexType fromBlockIndex) const
{
    auto it = m_Antecedents.find(fromBlockIndex);
    if (it == m_Antecedents.end()) {
        return BlockIndexInputIteratorBase(nullptr);
    }
    return BlockIndexInputIteratorBase(new BlockVectorSource(it->second));
}

BlockIndexInputIteratorBase ExplicitPrecedence::Successors(IndexType fromBlockIndex) const
{
    throw std::logic_error("not supported with explicit precedence");
    return BlockIndexInputIteratorBase(nullptr);
}

void ExplicitPrecedence::AddPrecedenceConstraint(IndexType fromBlockIndex, IndexType toBlockIndex)
{
    m_Antecedents[fromBlockIndex].push_back(toBlockIndex);
}

////////////////////////////////////////////////////////////////////////////////

VecBlockSource::VecBlockSource(std::vector<IndexType>&& blocks)
    : m_Blocks(blocks)
    , m_Index(0)
{
}

VecBlockSource::~VecBlockSource()
{
}

IndexType VecBlockSource::Next()
{
    return m_Blocks[m_Index++];
}

bool VecBlockSource::HasMore() const
{
    return m_Index < m_Blocks.size();
}

////////////////////////////////////////////////////////////////////////////////

PseudoSolverSolveInfo::PseudoSolverSolveInfo()
{
#ifdef MVD_USE_GMP
    mpz_init(ContainedValue);
#endif
}

PseudoSolverSolveInfo::~PseudoSolverSolveInfo()
{
#ifdef MVD_USE_GMP
    mpz_clear(ContainedValue);
#endif
}

std::string mvd::mineflow::PseudoSolverSolveInfoToString(const PseudoSolverSolveInfo& info)
{
    std::ostringstream os;
    os << "PseudoSolverSolveInfo: " << info.NumNodes << " input nodes" << std::endl;
    os << "  Contained : " << info.NumContainedNodes << " / " << info.NumNodes << std::endl;
    os << "  Used : " << info.NumUsedPrecedenceConstraints << " precedence constraints" << std::endl;
#ifdef MVD_USE_GMP
    char* str = mpz_get_str(NULL, 10, info.ContainedValue);
    void (*freefunc)(void *, size_t);
    mp_get_memory_functions(NULL, NULL, &freefunc);

    os << "  Value : " << str << std::endl;
    freefunc(str, std::strlen(str) + 1);
#else
    os << "  Value : " << info.ContainedValue << std::endl;
#endif
    os << "  Elapsed : " << std::setw(8) << std::fixed << std::setprecision(2) << info.ElapsedSeconds << "s";
    return os.str();
}

PseudoSolver::PseudoSolver(
        std::shared_ptr<const IPrecedenceConstraints> pre,
        const IBlockValues* values)
    : m_NodePool(std::make_unique<NodePool>(pre))
    , m_PrecedenceArcs(std::make_unique<PrecedenceArcPool>())
    , m_NodePoolHasBeenInitialized(false)
    , m_MinCutHasBeenSolved(false)
    , m_PrecedenceConstraints(pre)
{
#ifdef MVD_USE_GMP
    mpz_init(m_PrevExcess);
#endif

    if (!pre) {
        throw std::invalid_argument("precedence constraints must be defined");
    }

    if (values) {
        UpdateValues(values);
    }
}

PseudoSolver::PseudoSolver(
        std::shared_ptr<const IPrecedenceConstraints> pre,
        std::shared_ptr<const IBlockValues> values)
    : PseudoSolver(pre, values.get())
{
}

PseudoSolver::~PseudoSolver()
{
#ifdef MVD_USE_GMP
    mpz_clear(m_PrevExcess);
#endif
}

void PseudoSolver::UpdateValues(const IBlockValues* values)
{
    if (!values) {
        throw std::invalid_argument("values must be non null");
    }
    if (values->NumBlocks() != m_NodePool->NumNodes()) {
        throw std::invalid_argument("argument num blocks disagree");
    }

    if (!m_NodePoolHasBeenInitialized) {
        for (IndexType nodeIndex = 0; nodeIndex < m_NodePool->NumNodes(); nodeIndex++) {
            m_NodePool->InitializeNodeValue(nodeIndex, [&](ValueType* v){
                values->BlockValue(nodeIndex, v);
            });
        }
        m_NodePoolHasBeenInitialized = true;
    } else {
        // TODO, only need to renormalize some of the branches.
        // But the questions currently remaining are:
        //  - How do I best identify them?
        //  - How do I update labels correctly? and other flow values and such
        //  - For the precedence arc pool what arcs are lost?
        //  - Probably more
        // So for now we just reset everything..
        m_NodePool = std::make_unique<NodePool>(m_PrecedenceConstraints);
        m_PrecedenceArcs = std::make_unique<PrecedenceArcPool>();

        for (IndexType nodeIndex = 0; nodeIndex < m_NodePool->NumNodes(); nodeIndex++) {
            m_NodePool->InitializeNodeValue(nodeIndex, [&](ValueType* v){
                values->BlockValue(nodeIndex, v);
            });
        }
    }
    m_MinCutHasBeenSolved = false;
}

void PseudoSolver::UpdateValues(std::shared_ptr<const IBlockValues> values)
{
    return UpdateValues(values.get());
}

Node* PseudoSolver::WalkToRoot(Node* strongNode, Node* weakNode, Arc* newArc)
{
    Node* current = strongNode;
    Node* newParent = weakNode;

    Node* oldParent = nullptr;
    Arc* oldArc = nullptr;

    Arc* q = current->ToRoot;
    while (q->Tail && q->Head) {
        oldArc = current->ToRoot;
        current->ToRoot = newArc;
        oldParent = (q->Tail == current) ? q->Head : q->Tail;

        oldParent->RemoveChild(current);
        newParent->AddChild(current);

        newParent = current;
        current = oldParent;
        newArc = oldArc;
        q = current->ToRoot;
    }

    current->ToRoot = newArc;
    newParent->AddChild(current);
    return current;
}

void PseudoSolver::Split(impl::Node* current, impl::Node* parent, impl::Arc* arc)
{
#ifdef MVD_USE_GMP
    mpz_sub(current->Excess, current->Excess, arc->Flow);
    mpz_add(parent->Excess, parent->Excess, arc->Flow);
#else
    current->Excess -= arc->Flow;
    parent->Excess += arc->Flow;
#endif
    m_PrecedenceArcs->DeleteArc(arc);
    parent->Antecedents.OutOfTree.push_back(current);
    parent->RemoveChild(current);
    m_NodePool->ReconnectToRoot(current);
    m_NodePool->PushStrongRoot(current);
}

void PseudoSolver::PushFlow(Node* strongRoot)
{
#ifdef MVD_USE_GMP
    mpz_set_si(m_PrevExcess, 1);
#else
    m_PrevExcess = 1;
#endif
    Node* parent = nullptr;
    Node* current = strongRoot;
    while (true) {
        Arc* tr = current->ToRoot;
        Node* parent = current->ToRoot->Tail;
        if (parent == current) {
            parent = current->ToRoot->Head;
        }
#ifdef MVD_USE_GMP
        bool currentExcessGTZero = mpz_sgn(current->Excess) == 1;
#else
        bool currentExcessGTZero = current->Excess > 0;
#endif
        if (currentExcessGTZero && parent) {
#ifdef MVD_USE_GMP
            mpz_set(m_PrevExcess, parent->Excess);
#else
            m_PrevExcess = parent->Excess;
#endif

            bool up = tr->Tail == current;
            if (up) {
#ifdef MVD_USE_GMP
                mpz_add(parent->Excess, parent->Excess, current->Excess);
                mpz_add(tr->Flow, tr->Flow, current->Excess);
                mpz_set_ui(current->Excess, 0);
#else
                parent->Excess += current->Excess;
                tr->Flow += current->Excess;
                current->Excess = 0;
#endif
            } else {
#ifdef MVD_USE_GMP
                bool trFlowGECurrentExcess = mpz_cmp(tr->Flow, current->Excess) >= 0;
#else
                bool trFlowGECurrentExcess = tr->Flow >= current->Excess;
#endif
                if (trFlowGECurrentExcess) {
#ifdef MVD_USE_GMP
                    mpz_add(parent->Excess, parent->Excess, current->Excess);
                    mpz_sub(tr->Flow, tr->Flow, current->Excess);
                    mpz_set_ui(current->Excess, 0);
#else
                    parent->Excess += current->Excess;
                    tr->Flow -= current->Excess;
                    current->Excess = 0;
#endif
                } else {
                    Split(current, parent, tr);
                }
            }
        } else {
            break;
        }
        current = parent;
    }

#ifdef MVD_USE_GMP
    bool currentExcessGTZero = mpz_sgn(current->Excess) == 1;
#else
    bool currentExcessGTZero = current->Excess > 0;
#endif
    if (currentExcessGTZero) {
#ifdef MVD_USE_GMP
        bool prevExcessLEZero = mpz_sgn(m_PrevExcess) <= 0;
#else
        bool prevExcessLEZero = m_PrevExcess <= 0;
#endif
        if (prevExcessLEZero) {
            m_NodePool->PushStrongRoot(current);
        }
    }
}

void PseudoSolver::Merge(Node* strongNode, Node* weakNode)
{
    Arc* newArc = m_PrecedenceArcs->NewArc(strongNode, weakNode);
    Node* strongRoot = WalkToRoot(strongNode, weakNode, newArc);
    PushFlow(strongRoot);
}

void PseudoSolver::ProcessChildren(Node* node)
{
    MVD_ASSERT(node);

    // Loop over the remaining children (might be all of them!)
    while (node->NextScan) {
        MVD_ASSERT(node->NextScan->Label >= node->Label);
        if (node->NextScan->Label == node->Label) {
            return;
        }

        node->NextScan = node->NextScan->NextChild;
    }

    m_NodePool->IncrementLabel(node);
    node->Antecedents.NextArc = 0;
}

void PseudoSolver::ProcessStrongRoot(Node* strongRoot)
{
    IndexType inLabel = strongRoot->Label;
    strongRoot->NextScan = strongRoot->FirstChild;

    Node* weak = strongRoot->FindWeakAbove();
    if (weak) {
        Merge(strongRoot, weak);
        return;
    }

    Node* strongNode = strongRoot;
    ProcessChildren(strongRoot);

    while (strongNode) {
        while (strongNode->NextScan) {
            Node* temp = strongNode->NextScan;
            strongNode->NextScan = strongNode->NextScan->NextChild;
            strongNode = temp;
            strongNode->NextScan = strongNode->FirstChild;

            weak = strongNode->FindWeakAbove();
            if (weak) {
                Merge(strongNode, weak);
                return;
            }

            ProcessChildren(strongNode);
        }

        Node* temp = strongNode->ToRoot->Head;
        if (temp == strongNode) {
            temp = strongNode->ToRoot->Tail;
        } else {
            MVD_ASSERT(strongNode->ToRoot->Tail == strongNode);
        }
        strongNode = temp;

        if (strongNode) {
            ProcessChildren(strongNode);
        }
    }

    MVD_ASSERT(strongRoot->Label > inLabel);
    m_NodePool->PushStrongRoot(strongRoot);
}

void PseudoSolver::Solve(PseudoSolverSolveInfo* info)
{
    Node* strongRoot;

    auto start = std::chrono::steady_clock::now();
    if (m_NodePoolHasBeenInitialized) {
        while (m_NodePool->NextStrongRoot(&strongRoot)) {
            ProcessStrongRoot(strongRoot);
        }
    }
    auto end = std::chrono::steady_clock::now();

    if (info) {
        info->ElapsedSeconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;
        info->NumNodes = m_NodePool->NumNodes();
        info->NumContainedNodes = 0;
        ValueType temp;

#ifdef MVD_USE_GMP
        mpz_set_si(info->ContainedValue, 0);
        mpz_init(temp);
#else
        info->ContainedValue = 0;
#endif
        info->NumUsedPrecedenceConstraints = m_PrecedenceArcs->NumUsed();
        for (IndexType nodeIndex = 0; nodeIndex < NumNodes(); nodeIndex++) {
            if (m_NodePool->InMinimumCut(nodeIndex)) {
                info->NumContainedNodes++;
                m_NodePool->GetNodeValue(nodeIndex, &temp);
#ifdef MVD_USE_GMP
                mpz_add(info->ContainedValue, info->ContainedValue, temp);
#else
                info->ContainedValue += temp;
#endif
            }
        }
#ifdef MVD_USE_GMP
        mpz_clear(temp);
#endif
    }
    m_MinCutHasBeenSolved = true;
}

void PseudoSolver::SolveLargest(PseudoSolverSolveInfo* info)
{
    auto start = std::chrono::steady_clock::now();
    if (!m_MinCutHasBeenSolved) {
        Solve();
    }
    
    constexpr uint8_t UNKNOWN = 10;
    constexpr uint8_t DEFINITELY_IN = 1;
    constexpr uint8_t DEFINITELY_OUT = 0;
    constexpr uint8_t IN_PROCESS = 2;

    IndexType numNodes = NumNodes();

    m_LargestSolution.assign(numNodes, UNKNOWN);

    std::vector<std::vector<IndexType>> toCheck;
    for (IndexType nodeIndex = 0; nodeIndex < numNodes; nodeIndex++) {
        if (m_LargestSolution[nodeIndex] == UNKNOWN) {
            if (m_NodePool->InMinimumCut(nodeIndex)) {
                m_LargestSolution[nodeIndex] = DEFINITELY_IN;
            } else {
                Node* start = m_NodePool->GetNode(nodeIndex);
                Node* n = start;
                Arc* q = n->ToRoot;
                while (q->Tail && q->Head) {
                    n = (q->Tail == n) ? q->Head : q->Tail;
                    q = n->ToRoot;
                }

#ifdef MVD_USE_GMP
                bool nExcessZero = mpz_sgn(n->Excess) == 0;
#else
                bool nExcessZero = n->Excess == 0;
#endif
                uint8_t setBranchTo = (nExcessZero) ? IN_PROCESS : DEFINITELY_OUT;
                //MVD_ASSERT(n->Excess <= 0);

                std::vector<IndexType> thisBranch;

                n->ForNodeAndChildren([&](Node* v){
                    IndexType vi = m_NodePool->GetNodeIndex(v);
                    m_LargestSolution[vi] = setBranchTo;
                    if (nExcessZero) {
                        thisBranch.push_back(vi);
                    }
                });
                if (nExcessZero) {
                    toCheck.emplace_back(std::move(thisBranch));
                }
            }
        }
    }

    // These roots have an excess of zero
    if (!toCheck.empty()) {
        auto buffer = m_PrecedenceConstraints->GetNewSearchBuffer();
        for (auto & branch : toCheck) {
            uint8_t whatItIs = UNKNOWN;
            for (auto & v : branch) {
                if (m_LargestSolution[v] != IN_PROCESS) {
                    whatItIs = m_LargestSolution[v];
                    break;
                }
            }
            if (whatItIs != UNKNOWN) {
                for (auto & v : branch) {
                    m_LargestSolution[v] = whatItIs;
                }
            } else {
                bool foundDefOut = false;

                std::vector<IndexType> thisSearch;

                for (auto & l : branch) {
                    thisSearch.push_back(l);

                    // do search
                    if (!foundDefOut) {
                        m_PrecedenceConstraints->PartialReachableAntecedents(l,
                        [&](IndexType v){
                            if (m_LargestSolution[v] == DEFINITELY_OUT) {
                                foundDefOut = true;
                                return false;
                            } else if (m_LargestSolution[v] == DEFINITELY_IN) {
                                return false;
                            } else if (m_LargestSolution[v] == IN_PROCESS) {
                                thisSearch.push_back(v);
                                return !foundDefOut;
                            }
                            MVD_ASSERT(false);
                            return false;
                        }, buffer.get());
                    }
                }

                uint8_t setSearchTo = (foundDefOut) ? DEFINITELY_OUT : DEFINITELY_IN;
                for (auto & v : thisSearch) {
                    m_LargestSolution[v] = setSearchTo;
                }
            }
        }
    }
    auto end = std::chrono::steady_clock::now();

    if (info) {
        info->ElapsedSeconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;
        info->NumNodes = m_NodePool->NumNodes();
        info->NumContainedNodes = 0;

        ValueType temp;
#ifdef MVD_USE_GMP
        mpz_set_si(info->ContainedValue, 0);
        mpz_init(temp);
#else
        info->ContainedValue = 0;
#endif
        info->NumUsedPrecedenceConstraints = 0;
        for (IndexType nodeIndex = 0; nodeIndex < numNodes; nodeIndex++) {
            if (m_LargestSolution[nodeIndex] > 0) {
                info->NumContainedNodes++;
                m_NodePool->GetNodeValue(nodeIndex, &temp);
#ifdef MVD_USE_GMP
                mpz_add(info->ContainedValue, info->ContainedValue, temp);
#else
                info->ContainedValue += temp;
#endif
            }
        }
#ifdef MVD_USE_GMP
        mpz_clear(temp);
#endif
    }
}

////////////////////////////////////////////////////////////////////////////////

IndexType PseudoSolver::NumNodes() const
{
    return m_NodePool->NumNodes();
}

bool PseudoSolver::InMinimumCut(IndexType nodeIndex) const
{
    return m_NodePool->InMinimumCut(nodeIndex);
}

bool PseudoSolver::InLargestMinimumCut(IndexType nodeIndex) const
{
    if (m_LargestSolution.size() == m_NodePool->NumNodes()) {
        return m_LargestSolution.at(nodeIndex) > 0;
    }
    throw std::runtime_error("call solve largest");
    return false;
}

////////////////////////////////////////////////////////////////////////////////

void Node::AddChild(Node* child)
{
    MVD_ASSERT(child->NextChild == nullptr);
    child->NextChild = FirstChild;
    FirstChild = child;
}

void Node::RemoveChild(Node* child)
{
    MVD_ASSERT(FirstChild != nullptr);

    if (FirstChild == child) {
        FirstChild = child->NextChild;
        child->NextChild = nullptr;
        return;
    }

    Node* current = FirstChild;
    MVD_ASSERT(current->NextChild != nullptr);
    while (current->NextChild != child) {
        MVD_ASSERT(current->NextChild != nullptr);
        current = current->NextChild;
    }
    MVD_ASSERT(current->NextChild == child);

    current->NextChild = child->NextChild;
    child->NextChild = nullptr;
}

void Node::IncrementLabel()
{
    Label++;
}

void Node::ForNodeAndChildren(std::function<void(Node*)> cback)
{
    cback(this);
    Node* v = FirstChild;
    while (v) {
        v->ForNodeAndChildren(cback);
        v = v->NextChild;
    }
}

Node* Node::FindWeakAbove()
{
    if (Antecedents.Init) {
        Antecedents.Init->InitPrecedence(this);
        Antecedents.Init = nullptr;
    }
    for (IndexType i = Antecedents.NextArc; 
         i < static_cast<IndexType>(Antecedents.OutOfTree.size());
         i++) {

        Node* to = Antecedents.OutOfTree[i];
        if (to->Label == (Label - 1)) {
            Antecedents.NextArc = i;
            Antecedents.OutOfTree[i] = Antecedents.OutOfTree.back();
            Antecedents.OutOfTree.pop_back();
            return to;
        }
    }

    Antecedents.NextArc = Antecedents.OutOfTree.size();
    return nullptr;
}

////////////////////////////////////////////////////////////////////////////////


PrecedenceArcPool::PrecedenceArcPool() 
    : m_NumUsed(0) 
{
#ifdef MVD_USE_GMP
    m_Remaining = 0;
#endif
}

PrecedenceArcPool::~PrecedenceArcPool()
{
}

Arc* PrecedenceArcPool::NewArc(Node* tail, Node* head)
{
#ifdef MVD_USE_GMP
    if (m_Remaining == 0) {
        m_Remaining = N - 1;
        m_ArcSets.emplace_front();
        m_Set = &(m_ArcSets.front());
    }

    Arc* arc = &(m_Set->Arcs[m_Remaining]);
    m_Remaining--;
    mpz_set_si(arc->Flow, 0);
#else
    Arc* arc = m_ObjectPool.Alloc<Arc>();
    arc->Flow = 0;
#endif

    arc->Tail = tail;
    arc->Head = head;
    m_NumUsed++;
    return arc;
}
void PrecedenceArcPool::DeleteArc(Arc* arc)
{
    // TODO may be efficient to reclaim this memory when allowing resolves 
    
    //MVD_ASSERT(arc->Head);
    //MVD_ASSERT(arc->Tail);
    arc->Head = nullptr;
    arc->Tail = nullptr;
#ifdef MVD_USE_GMP
    mpz_set_si(arc->Flow, 0);
#else
    arc->Flow = 0;
#endif
}

IndexType PrecedenceArcPool::NumUsed() const
{
    return m_NumUsed;
}

////////////////////////////////////////////////////////////////////////////////

NodePool::NodePool(std::shared_ptr<const IPrecedenceConstraints> pre)
    : m_PrecedenceConstraints(pre)
    , m_NumNodes(pre->NumBlocks())
{
    m_Nodes.resize(m_NumNodes);
    m_RootArcs.resize(m_NumNodes);
    m_Buckets.resize(2);
    m_LabelCount.resize(2);

    for (IndexType nodeIndex = 0; nodeIndex < m_NumNodes; nodeIndex++) {
        Node* node = &m_Nodes[nodeIndex];
        Arc* arc = &m_RootArcs[nodeIndex];

#ifdef MVD_USE_GMP
        mpz_init_set_si(node->Excess, 0);
        mpz_init(arc->Flow);
#else
        node->Excess = 0;
#endif
        node->ToRoot = arc;
        node->Label = 0;
        node->FirstChild = nullptr;
        node->NextChild = nullptr;
        node->NextScan = nullptr;
        node->Antecedents.Init = this;
        node->Antecedents.NextArc = 0;
    }
}

NodePool::~NodePool()
{
#ifdef MVD_USE_GMP
    for (IndexType nodeIndex = 0; nodeIndex < m_NumNodes; nodeIndex++) {
        Node* node = &m_Nodes[nodeIndex];
        Arc* arc = &m_RootArcs[nodeIndex];

        mpz_clear(node->Excess);
        mpz_clear(arc->Flow);
    }
#endif
}

void NodePool::ReconnectToRoot(Node* node)
{
    IndexType nodeIndex = node - &m_Nodes[0];
    node->ToRoot = &m_RootArcs[nodeIndex];
}

IndexType NodePool::NumNodes() const
{
    return m_NumNodes;
}

bool NodePool::InMinimumCut(IndexType nodeIndex) const
{
    return m_Nodes[nodeIndex].Label == m_NumNodes;
}

void NodePool::InitializeNodeValue(IndexType nodeIndex, std::function<void(ValueType*)> getValue)
{
    Node* node = &m_Nodes[nodeIndex];
    Arc* arc = &m_RootArcs[nodeIndex];

    getValue(&node->Excess);
#ifdef MVD_USE_GMP
    bool nodeExcessGTZero = mpz_sgn(node->Excess) > 0;
#else
    bool nodeExcessGTZero = node->Excess > 0;
#endif

    if (nodeExcessGTZero) {
        node->Label = 1;
        m_LabelCount[1]++;
        PushStrongRoot(node);

        arc->Tail = nullptr;
        arc->Head = node;
#ifdef MVD_USE_GMP
        mpz_set(arc->Flow, node->Excess);
#else
        arc->Flow = node->Excess;
#endif
    } else {
        node->Label = 0;
        m_LabelCount[0]++;

        arc->Tail = node;
        arc->Head = nullptr;
#ifdef MVD_USE_GMP
        mpz_neg(arc->Flow, node->Excess);
#else
        arc->Flow = -node->Excess;
#endif
    }
}

Node* NodePool::GetNode(IndexType nodeIndex)
{
    return &m_Nodes[nodeIndex];
}

IndexType NodePool::GetNodeIndex(const Node* n) const
{
    return static_cast<IndexType>(n - &m_Nodes[0]);
}

void NodePool::GetNodeValue(IndexType nodeIndex, ValueType* value) const
{
    const Arc* arc = &m_RootArcs[nodeIndex];

    if (arc->Tail) {
#ifdef MVD_USE_GMP
        mpz_neg(*value, arc->Flow);
#else
        *value = -arc->Flow;
#endif
    } else {
#ifdef MVD_USE_GMP
        mpz_set(*value, arc->Flow);
#else
        *value = arc->Flow;
#endif
    }
}

void NodePool::PushStrongRoot(Node* node)
{
    if (static_cast<IndexType>(m_Buckets.size()) <= node->Label) {
        m_Buckets.resize(node->Label + 1);
    }
    m_Buckets[node->Label].push(node);
}

bool NodePool::NextStrongRoot(Node** nodep)
{
    for (int64_t i = m_Buckets.size() - 1; i > 0; i--) {
        std::queue<Node*>& queue = m_Buckets[i];
        if (!queue.empty()) {
            if (m_LabelCount[i - 1] > 0) {
                *nodep = queue.front();
                queue.pop();
                return true;
            } else {

                while (!queue.empty()) {
                    queue.front()->ForNodeAndChildren([&](Node* v){
                        m_LabelCount[v->Label]--;
                        v->Label = m_NumNodes;
                    });
                    queue.pop();
                }
            }
        } else {
            m_Buckets.pop_back();
        }
    }
    
    if (m_Buckets[0].empty()) {
        *nodep = nullptr;
        return false;
    }

    std::queue<Node*>& queue = m_Buckets[0];
    while (!queue.empty()) {
        Node* root = queue.front();
        queue.pop();

        MVD_ASSERT(root->Label == 0);
        IncrementLabel(root);
        PushStrongRoot(root);
    }

    MVD_ASSERT(!m_Buckets[1].empty());
    *nodep = m_Buckets[1].front();
    m_Buckets[1].pop();
    return true;
}

void NodePool::IncrementLabel(Node* node)
{
    MVD_ASSERT(m_LabelCount[node->Label] > 0);
    m_LabelCount[node->Label]--;
    node->IncrementLabel();
    if (static_cast<IndexType>(m_LabelCount.size()) <= node->Label) {
        m_LabelCount.resize(node->Label + 1);
    }
    m_LabelCount[node->Label]++;
}

void NodePool::InitPrecedence(Node* node)
{
    IndexType nodeIndex = node - &m_Nodes[0];
    node->Antecedents.OutOfTree.reserve(
            m_PrecedenceConstraints->ApproxNumAntecedents(nodeIndex));
    for (auto & targetIndex : m_PrecedenceConstraints->Antecedents(nodeIndex)) {
        node->Antecedents.OutOfTree.push_back(&m_Nodes[targetIndex]);
    }
}

////////////////////////////////////////////////////////////////////////////////

SolveLargestValuesAdapter::SolveLargestValuesAdapter(
        std::shared_ptr<const IBlockValues> values)
    : m_Values(values)
{
    if (!m_Values) {
        throw std::invalid_argument("must supply values to solve largest adapter");
    }
    IndexType n = m_Values->NumBlocks();

#ifdef MVD_USE_GMP
    mpz_init_set_si(m_NumNonNegativeBlocks, 0);
    ValueType temp;
    mpz_init(temp);
    for (IndexType i = 0; i < n; i++) {
        m_Values->BlockValue(i, &temp);
        if (mpz_sgn(temp) >= 0) {
            mpz_add_ui(m_NumNonNegativeBlocks, m_NumNonNegativeBlocks, 1);
        }
    }
    mpz_clear(temp);
#else
    m_NumNonNegativeBlocks = 0;
    ValueType temp;
    for (IndexType i = 0; i < n; i++) {
        m_Values->BlockValue(i, &temp);
        if (temp >= 0) {
            m_NumNonNegativeBlocks++;
        }
    }
    m_NumNonNegativeBlocks++;
#endif
}

SolveLargestValuesAdapter::~SolveLargestValuesAdapter()
{
#ifdef MVD_USE_GMP
    mpz_clear(m_NumNonNegativeBlocks);
#endif
}

IndexType SolveLargestValuesAdapter::NumBlocks() const
{
    return m_Values->NumBlocks();
}

void SolveLargestValuesAdapter::BlockValue(IndexType blockIndex, ValueType* value) const
{
#ifdef MVD_USE_GMP
    m_Values->BlockValue(blockIndex, value);
    mpz_mul(*value, *value, m_NumNonNegativeBlocks);
    if (mpz_sgn(*value) >= 0) {
        mpz_add_ui(*value, *value, 1);
    }
#else
    ValueType v;
    m_Values->BlockValue(blockIndex, &v);

    if (v >= 0) {
        *value = v * m_NumNonNegativeBlocks + 1;
    } else { 
        *value = v * m_NumNonNegativeBlocks;
    }
#endif
}

////////////////////////////////////////////////////////////////////////////////
#ifdef MVD_MINEFLOW_TESTS
#ifdef MVD_MINEFLOW_EXE
#error "Must not define both 'MVD_MINEFLOW_EXE' and 'MVD_MINEFLOW_TESTS'"
#endif
#include <sstream>
#include <fstream>

// Oh the lengths I go to to avoid a dependency.
// The following just sets up a simple test framework
class ITest
{
public:
    ITest(){};
    virtual ~ITest(){};
    virtual void DoTest() = 0;

    double m_ExtraElapsed = -1; // Hacky, but... yeah
};
struct RegisteredTest 
{
    std::shared_ptr<ITest> Test;
    std::string ClassName;
    std::string InstanceName;
};
class TestException : public std::exception 
{
public:
    virtual const std::string& GetTestFailedMessage() const = 0;
};
class TestSException : public TestException
{
public:
    const std::string& GetTestFailedMessage() const {
        return m_Message;
    }
    std::string m_Message;
};

static std::string FailedLine(const std::string& assertType,
        const std::string& file, int lineno) {
    std::ostringstream os;
    os << assertType << " FAILED: " << file << ":" << lineno << " failure" << std::endl;
    return os.str();
};

class AssertTrueException : public TestSException
{
public:
    AssertTrueException(const std::string& exprs, const std::string& file, 
            int lineno)
    {
        std::ostringstream os;
        os << FailedLine("ASSERT", file, lineno);
        os << "  (" << exprs << ") was false";
        m_Message = os.str();
    }
};
#define ASSERT_TRUE(expr)                                                       \
    if (!(expr)) {                                                              \
        throw AssertTrueException(#expr, __FILE__, __LINE__);                   \
    }

class AssertEqException : public TestSException
{
public:
    AssertEqException(const std::string& lhs, const std::string& rhs,
            const std::string& values,
            const std::string& file, int lineno)
    {
        std::ostringstream os;
        os << FailedLine("ASSERT_EQ", file, lineno);
        os << "  (" << lhs << " == " << rhs << ")" << std::endl;
        os << values;
        m_Message = os.str();
    }
};
#define ASSERT_EQ(lhs, rhs) {                                                       \
    auto TESTING_l = lhs;                                                           \
    auto TESTING_r = rhs;                                                           \
    if (TESTING_l != TESTING_r) {                                                   \
        std::ostringstream TESTING_OS;                                              \
        TESTING_OS << "lhs:  " << TESTING_l << std::endl << "rhs:  " << TESTING_r;  \
        throw AssertEqException(#lhs, #rhs, TESTING_OS.str(), __FILE__, __LINE__);  \
    }}

class AssertNearException : public TestSException
{
public:
    AssertNearException(const std::string& lhs, const std::string& rhs, 
            const std::string& epsilon,
            const std::string& values,
            const std::string& file, int lineno)
    {
        std::ostringstream os;
        os << FailedLine("ASSERT_NEAR", file, lineno);
        os << "  (std::abs(" << lhs << " - " << rhs << ") < " << epsilon << ")" << std::endl;
        os << values;
        m_Message = os.str();
    }
};
#define ASSERT_NEAR(lhs, rhs, epsilon) {                                            \
    auto TESTING_l = lhs;                                                           \
    auto TESTING_r = rhs;                                                           \
    auto diff = lhs - rhs;                                                          \
    if (std::abs(diff) > epsilon) {                                                 \
        std::ostringstream TESTING_OS;                                              \
        TESTING_OS << "lhs:  " << std::setprecision(16) << TESTING_l << std::endl;  \
        TESTING_OS << "rhs:  " << std::setprecision(16) << TESTING_r;               \
        throw AssertNearException(#lhs, #rhs, #epsilon,                             \
                TESTING_OS.str(), __FILE__, __LINE__);                              \
    }}

class AssertVecEqException : public TestSException
{
public:
    AssertVecEqException(const std::string& lhs, const std::string& rhs,
            const std::string& file, int lineno, const std::string& message)
    {
        std::ostringstream os;
        os << FailedLine("ASSERT_VEC_EQ", file, lineno);
        os << "  " << message;
        m_Message = os.str();
    }
};
#define ASSERT_VEC_EQ(lhs, rhs) {                                                           \
    if (lhs.size() != rhs.size()) {                                                         \
        std::ostringstream TESTING_OS;                                                      \
        TESTING_OS << "vector sizes not equal: " << lhs.size() << " " << rhs.size();        \
        throw AssertVecEqException(#lhs, #rhs, __FILE__, __LINE__, TESTING_OS.str());       \
    }                                                                                       \
    for (size_t TESTING_I = 0; TESTING_I < lhs.size(); TESTING_I++) {                       \
        if (lhs[TESTING_I] != rhs[TESTING_I]) {                                             \
            std::ostringstream TESTING_OS;                                                  \
            TESTING_OS << "vector elements not equal: index: " << TESTING_I << std::endl;   \
            TESTING_OS << " lhs[" << TESTING_I << "]: " << lhs[TESTING_I] << std::endl;     \
            TESTING_OS << " rhs[" << TESTING_I << "]: " << rhs[TESTING_I] << std::endl;     \
            throw AssertVecEqException(#lhs, #rhs, __FILE__, __LINE__,                      \
                    TESTING_OS.str());                                                      \
        }                                                                                   \
    }                                                                                       \
}

class AssertValueTypeEqException : public TestSException
{
public:
    AssertValueTypeEqException(const std::string& lhs, const std::string& rhs,
            const std::string& values,
            const std::string& file, int lineno)
    {
        std::ostringstream os;
        os << FailedLine("ASSERT_VALUETYPE_EQ_INT", file, lineno);
        os << "  (" << lhs << " == " << rhs << ")" << std::endl;
        os << values;
        m_Message = os.str();
    }
};
#ifdef MVD_USE_GMP
#define ASSERT_VALUETYPE_EQ_INT(lhs, rhs) {                                                 \
    if (mpz_cmp_si(lhs, rhs) != 0) {                                                        \
        char* lhss = mpz_get_str(NULL, 10, lhs);                                            \
        std::ostringstream TESTING_OS;                                                      \
        TESTING_OS << "lhs:  " << lhss << std::endl << "rhs:  " << rhs;                     \
        void (*freefunc)(void *, size_t);                                                   \
        mp_get_memory_functions(NULL, NULL, &freefunc);                                     \
        freefunc(lhss, std::strlen(lhss) + 1);                                              \
        throw AssertValueTypeEqException(#lhs, #rhs, TESTING_OS.str(), __FILE__, __LINE__); \
    }                                                                                       \
}
#else
#define ASSERT_VALUETYPE_EQ_INT(lhs, rhs) {                                                     \
    auto TESTING_l = lhs;                                                                   \
    auto TESTING_r = rhs;                                                                   \
    if (TESTING_l != TESTING_r) {                                                           \
        std::ostringstream TESTING_OS;                                                      \
        TESTING_OS << "lhs:  " << TESTING_l << std::endl << "rhs:  " << TESTING_r;          \
        throw AssertValueTypeEqException(#lhs, #rhs, TESTING_OS.str(), __FILE__, __LINE__); \
    }                                                                                       \
}
#endif




static std::vector<RegisteredTest>& RegisteredTests() 
{
    static std::vector<RegisteredTest> tests;
    return tests;
}
static std::string RegisterTest(const std::string& className, const std::string& instanceName,
        std::shared_ptr<ITest>&& test) 
{
    auto& ts = RegisteredTests();
    ts.emplace_back();
    auto& t = ts.back();
    t.Test = std::move(test);
    t.ClassName = className;
    t.InstanceName = instanceName;

    std::string testName = className + "_" + instanceName;
    return testName;
}

#define TEST_NAME_(class_name, instance_name) class_name ## instance_name ## _iinstance

#define TEST(class_name, instance_name)                                         \
class TEST_NAME_(class_name, instance_name) : public ITest                      \
{                                                                               \
public:                                                                         \
    TEST_NAME_(class_name, instance_name)(){};                                  \
    virtual ~TEST_NAME_(class_name, instance_name)(){};                         \
    virtual void DoTest();                                                      \
    static std::string s_TestName;                                              \
};                                                                              \
std::string TEST_NAME_(class_name, instance_name)::s_TestName =                 \
    RegisterTest(#class_name, #instance_name,                                   \
            std::make_shared<TEST_NAME_(class_name, instance_name)>());         \
void TEST_NAME_(class_name, instance_name)::DoTest()

// MVD_MINEFLOW_TESTS_BEGIN

////////////////////////////////////////////////////////////////////////////////

TEST(Vector, SizeOf) 
{
    ASSERT_EQ(sizeof(double) * 3, sizeof(VectorBase<double, 3>));
    ASSERT_EQ(sizeof(double) * 2, sizeof(VectorBase<double, 2>));
    ASSERT_EQ(sizeof(float) * 3, sizeof(VectorBase<float, 3>));
    ASSERT_EQ(sizeof(int) * 2, sizeof(VectorBase<int, 2>));
}

TEST(Vector, Properties) 
{
    bool polymorphic = std::is_polymorphic<VectorBase<double, 3>>::value;
    ASSERT_TRUE(!polymorphic);
}

TEST(Vector, BasicConstructor) 
{
    Vector3D a(1.2, -13.4, 5.41);
    ASSERT_EQ(1.2, a.x);
    ASSERT_EQ(-13.4, a.y);
    ASSERT_EQ(5.41, a.z);
}

TEST(Vector, Origin) 
{
    // By default a vector is initialized with garbage (for speed) however there
    // are convenience methods
    Vector3D origin = Vector3D::Origin();
    ASSERT_EQ(0.0, origin.x);
    ASSERT_EQ(0.0, origin.y);
    ASSERT_EQ(0.0, origin.z);
}

TEST(Vector, Axes) 
{
    Vector3D x_axis = Vector3D::XAxis();
    ASSERT_EQ(1.0, x_axis.x);
    ASSERT_EQ(0.0, x_axis.y);
    ASSERT_EQ(0.0, x_axis.z);

    Vector3D y_axis = Vector3D::YAxis();
    ASSERT_EQ(0.0, y_axis.x);
    ASSERT_EQ(1.0, y_axis.y);
    ASSERT_EQ(0.0, y_axis.z);

    Vector3D z_axis = Vector3D::ZAxis();
    ASSERT_EQ(0.0, z_axis.x);
    ASSERT_EQ(0.0, z_axis.y);
    ASSERT_EQ(1.0, z_axis.z);
}

TEST(Vector, VectorAddition)
{
    Vector3D a(1, 2, 3);
    Vector3D b(-1, 7.2, 0);

    Vector3D c = a + b;
    ASSERT_EQ(0, c.x);
    ASSERT_EQ(9.2, c.y);
    ASSERT_EQ(3, c.z);

    a += a;
    ASSERT_EQ(2, a.x);
    ASSERT_EQ(4, a.y);
    ASSERT_EQ(6, a.z);
}

TEST(Vector, VectorSubtraction)
{
    Vector3D a(12, 4, 2);
    Vector3D b(3, 8, 1);

    Vector3D c = a - b;
    ASSERT_EQ(9, c.x);
    ASSERT_EQ(-4, c.y);
    ASSERT_EQ(1, c.z);

    a -= a;
    ASSERT_EQ(0, a.x);
    ASSERT_EQ(0, a.y);
    ASSERT_EQ(0, a.z);
}

TEST(Vector, VectorAddSubtractConstant)
{
    Vector3D a(0.5, 12.5, -12.3);
    a += 0.5;
    ASSERT_EQ(1, a.x);
    ASSERT_EQ(13, a.y);
    ASSERT_EQ(-11.8, a.z);

    a -= 11.9;
    ASSERT_NEAR(-10.9, a.x, 0.0000001);
    ASSERT_NEAR(1.1, a.y, 0.0000001);
    ASSERT_NEAR(-23.7, a.z, 0.0000001);

    Vector3D b = a + 10.9;
    ASSERT_NEAR(0, b.x, 0.0000001);
    ASSERT_NEAR(12, b.y, 0.0000001);
    ASSERT_NEAR(-12.8, b.z, 0.0000001);
}

TEST(Vector, Assignment)
{
    Vector3D a(1.2, 4.2, -1.8);

    Vector3D b;
    b = a;
    ASSERT_EQ(1.2, b.x);
    ASSERT_EQ(4.2, b.y);
    ASSERT_EQ(-1.8, b.z);
}

TEST(Vector, ComparisonEquals)
{
    Vector3D a(1.3, 2.1, 0);
    Vector3D b(1.3, 2.1, 0);
    Vector3D c(1.1, 1.0, 0);

    ASSERT_EQ(a, b);
    ASSERT_TRUE(a == b);
    ASSERT_TRUE(!(a != b));
    ASSERT_TRUE(!(a == c));
    ASSERT_TRUE(a != c);
}

TEST(Vector, ComparisonLessthan)
{
    std::vector<Vector3D> points {
        Vector3D(1, 1, 0),
        Vector3D(1, 2, 0),
        Vector3D(2, 1, 0),
        Vector3D(0, 1, 0),
        Vector3D(3, 1, 0)
    };
    std::sort(points.begin(), points.end());

    ASSERT_EQ(Vector3D(0, 1, 0), points[0]);
    ASSERT_EQ(Vector3D(1, 1, 0), points[1]);
    ASSERT_EQ(Vector3D(1, 2, 0), points[2]);
    ASSERT_EQ(Vector3D(2, 1, 0), points[3]);
    ASSERT_EQ(Vector3D(3, 1, 0), points[4]);
}

TEST(Vector, Multiplication)
{
    Vector3D a(1.2, 2.4, 3.6);

    Vector3D b = a * 4.0;
    ASSERT_EQ(Vector3D(4.8, 9.6, 14.4), b);

    b *= (1.0/4.0);
    ASSERT_EQ(a, b);
}

TEST(Vector, Division)
{
    Vector3D c(12, 15, -3);
    c /= 3;
    ASSERT_EQ(Vector3D(4, 5, -1), c);

    Vector3D a = c / 4.0;
    ASSERT_EQ(Vector3D(1, 1.25, -0.25), a);
}

TEST(Vector, DotProduct)
{
    Vector3D a(1, 3, -5);
    Vector3D b(4, -2, -1);

    ASSERT_EQ(3, Dot(a, b));

    a = Vector3D(1, 0, 0);
    b = Vector3D(0, 1, 0);
    ASSERT_EQ(0, Dot(a, b));
}

TEST(Vector, Magnitude)
{
    Vector3D a(1, 1, 0);

    ASSERT_EQ(sqrt(2), Magnitude(a));

    Vector3D b(3, 4, 0);
    ASSERT_EQ(5, Magnitude(b));
    ASSERT_EQ(25, MagnitudeSquared(b));

    Vector3D c(12, 16, 25);
    ASSERT_NEAR(32.0156, Magnitude(c), 0.001);
    ASSERT_EQ(1025, MagnitudeSquared(c));
}

TEST(Vector, Theta)
{
    Vector3D a(1, 0, 0);
    Vector3D b(0, 1, 0);

    ASSERT_NEAR(1.5708, Theta(a, b), 0.0001);

    Vector3D c(1, 1, 0);
    ASSERT_NEAR(0.7854, Theta(a, c), 0.0001);
    ASSERT_NEAR(0.7854, Theta(c, b), 0.0001);
}

TEST(Vector, Cross)
{
    Vector3D a(1, 0, 0);
    Vector3D b(0, 1, 0);

    Vector3D c = Cross(a, b);
    ASSERT_EQ(Vector3D(0, 0, 1), c);

    Vector3D d(3, -3, 1);
    Vector3D e(4, 9, 2);
    ASSERT_EQ(Vector3D(-15, -2, 39), Cross(d, e));
}

TEST(Vector, Normalize)
{
    Vector3D a(1.5, 0, 0);
    Normalize(a);
    ASSERT_EQ(Vector3D(1, 0, 0), a);

    Vector3D b(22.6, 22.6, 0);
    Normalize(b);
    ASSERT_EQ(Vector3D(sqrt(2.0)/2, sqrt(2.0)/2, 0), b);
}

TEST(Vector, LeftRight)
{
    Vector2D a(0, 0);
    Vector2D b(1, 1);

    ASSERT_TRUE(IsLeft(a, b, Vector2D(0, 1)));
    ASSERT_TRUE(!(IsRight(a, b, Vector2D(0, 1))));
    ASSERT_TRUE(!(IsCollinear(a, b, Vector2D(0, 1))));

    ASSERT_TRUE(!(IsLeft(a, b, Vector2D(2, 2))));
    ASSERT_TRUE(!(IsRight(a, b, Vector2D(2, 2))));
    ASSERT_TRUE(IsCollinear(a, b, Vector2D(2, 2)));

    ASSERT_TRUE(!(IsLeft(a, b, Vector2D(20, 1))));
    ASSERT_TRUE(IsRight(a, b, Vector2D(20, 1)));
    ASSERT_TRUE(!(IsCollinear(a, b, Vector2D(20, 1))));
}

TEST(Vector, TriArea)
{
    Vector2D a(-1, 0);
    Vector2D b(0, 1);
    Vector2D c(1, 0);

    ASSERT_EQ(-1, TriArea(a, b, c));
    ASSERT_EQ(1, TriArea(c, b, a));
}

TEST(Vector, InOut)
{
    Vector3D a(0, 1, 2), b;

    std::ostringstream os;
    os << a;
    ASSERT_EQ("{0, 1, 2}", os.str());

    std::istringstream is(os.str());
    is >> b;

    ASSERT_EQ(0, b.x);
    ASSERT_EQ(1, b.y);
    ASSERT_EQ(2, b.z);
}

TEST(Vector, Input)
{
    Vector3D a;

    std::istringstream is1(" {0.1, 0.5, -12.4}  ");
    is1 >> a;

    ASSERT_EQ(0.1, a.x);
    ASSERT_EQ(0.5, a.y);
    ASSERT_EQ(-12.4, a.z);
}

TEST(Angles, ToDegrees) 
{
    ASSERT_NEAR(ToDegrees(3.14159265), 180.0, 0.00001);
}

TEST(Angles, ToRadians)
{
    ASSERT_NEAR(ToRadians(180.0), 3.14159265, 0.00001);
}

TEST(Linspace, Base)
{
    std::vector<double> a(11);
    InplaceLinspace(a.begin(), a.end(), 0, 100);

    ASSERT_EQ(0, a[0]);
    ASSERT_EQ(10, a[1]);
    ASSERT_EQ(20, a[2]);
    ASSERT_EQ(30, a[3]);
    ASSERT_EQ(100, a[10]);
}

TEST(Linspace, NonZeroStart)
{
    std::vector<double> a(10);
    InplaceLinspace(a.begin(), a.end(), 20.0, 34.4);

    ASSERT_EQ(20.0, a[0]);
    ASSERT_EQ(21.6, a[1]);
    ASSERT_EQ(32.8, a[8]);
    ASSERT_EQ(34.4, a[9]);
}

TEST(Linspace, NegativeRange)
{
    std::vector<double> a(20);
    InplaceLinspace(a.begin(), a.end(), 83.1, -10.0);


    ASSERT_NEAR( 83.1, a[0], 0.0001);
    ASSERT_NEAR( 78.2, a[1], 0.0001);
    ASSERT_NEAR( -5.1, a[18], 0.0001);
    ASSERT_NEAR(-10.0, a[19], 0.0001);
}

TEST(Linspace, Generator)
{
    double i = 0;

    for (auto v : Linspace(0.0, 1.0, 11)) {
        ASSERT_NEAR(i / 10, v, 0.0001);
        i++;
    }
}

TEST(Linspace, GeneratorGauss)
{
    double sum = 0;
    for (auto v : Linspace(0, 100, 101)) {
        sum += v;
    }
    ASSERT_NEAR(5050, sum, 0.00001);
}

TEST(Linspace, GeneratorNegative)
{
    double sum = 0;
    for (auto v : Linspace(50, -50, 10)) {
        sum += v;
    }
    ASSERT_EQ(0, sum);
}

TEST(Block, 1DIndices)
{
    BlockDefinition def = BlockDefinition::UnitModel(10, 8, 5);

    ASSERT_EQ(0, def.XIndex(0));
    ASSERT_EQ(0, def.YIndex(0));
    ASSERT_EQ(0, def.ZIndex(0));

    ASSERT_EQ(1, def.XIndex(1));
    ASSERT_EQ(0, def.YIndex(1));
    ASSERT_EQ(0, def.ZIndex(1));

    ASSERT_EQ(2, def.XIndex(2));
    ASSERT_EQ(0, def.YIndex(2));
    ASSERT_EQ(0, def.ZIndex(2));

    ASSERT_EQ(0, def.XIndex(10));
    ASSERT_EQ(1, def.YIndex(10));
    ASSERT_EQ(0, def.ZIndex(10));

    ASSERT_EQ(9, def.XIndex(79));
    ASSERT_EQ(7, def.YIndex(79));
    ASSERT_EQ(0, def.ZIndex(79));

    ASSERT_EQ(0, def.XIndex(80));
    ASSERT_EQ(0, def.YIndex(80));
    ASSERT_EQ(1, def.ZIndex(80));

    ASSERT_EQ(1, def.XIndex(81));
    ASSERT_EQ(0, def.YIndex(81));
    ASSERT_EQ(1, def.ZIndex(81));

    ASSERT_EQ(3, def.XIndex(173));
    ASSERT_EQ(1, def.YIndex(173));
    ASSERT_EQ(2, def.ZIndex(173));
}

TEST(Block, 3DIndices)
{
    int64_t nx = 10, ny = 8, nz = 5;
    BlockDefinition def = BlockDefinition::UnitModel(nx, ny, nz);

    int64_t k = 0;
    for (int64_t z = 0; z < nz; z++) {
        for (int64_t y = 0; y < ny; y++) {
            for (int64_t x = 0; x < nx; x++) {
                ASSERT_EQ(x, def.XIndex(k));
                ASSERT_EQ(y, def.YIndex(k));
                ASSERT_EQ(z, def.ZIndex(k));
                ASSERT_EQ(k, def.GridIndex(x, y, z));

                k++;
            }
        }
    }
}

TEST(Block, NumBlocks)
{
    BlockDefinition def1 = BlockDefinition::UnitModel(10, 8, 5);
    ASSERT_EQ(400, def1.NumBlocks());

    BlockDefinition def2 = BlockDefinition::UnitModel(10, 10, 1);
    ASSERT_EQ(100, def2.NumBlocks());

    BlockDefinition def3 = BlockDefinition::UnitModel(32, 50, 20);
    ASSERT_EQ(32000, def3.NumBlocks());
}

TEST(Block, OffsetIndex)
{
    BlockDefinition def = BlockDefinition::UnitModel(10, 8, 5);

    ASSERT_EQ(5, def.OffsetIndex(0, 5, 0, 0));
    ASSERT_EQ(5, def.OffsetIndex(2, 3, 0, 0));
    ASSERT_EQ(15, def.OffsetIndex(5, 0, 1, 0));
    ASSERT_EQ(95, def.OffsetIndex(15, 0, 0, 1));
}

TEST(Block, InDef)
{
    int64_t nx = 10, ny = 8, nz = 5;
    BlockDefinition def = BlockDefinition::UnitModel(nx, ny, nz);

    ASSERT_EQ(false, def.InDef(-1, 0, 0));
    ASSERT_EQ(false, def.InDef(0, -1, 0));
    ASSERT_EQ(false, def.InDef(0, 0, -1));

    ASSERT_EQ(false, def.InDef(-1));
    ASSERT_EQ(false, def.InDef(nx * ny * nz));

    int64_t k = 0;
    for (int64_t z = 0; z < nz; z++) {
        for (int64_t y = 0; y < ny; y++) {
            for (int64_t x = 0; x < nx; x++) {
                ASSERT_EQ(true, def.InDef(x, y, z));
                ASSERT_EQ(true, def.InDef(k));

                k++;
            }
        }
    }
}

TEST(Precedence, Regular2DGrid45DegreePrecedenceBase)
{
    Regular2DGrid45DegreePrecedence pre(10, 6);
    ASSERT_EQ(60, pre.NumBlocks());

    uint64_t ne = 8 * 5 * 3 + 1 * 5 * 2 + 1 * 5 * 2;
    ASSERT_EQ(ne, pre.NumPrecedenceConstraints());

    std::vector<IndexType> to;
    pre.AntecedentsVector(5, &to);
    ASSERT_EQ(3, to.size());
    ASSERT_EQ(14, to[0]);
    ASSERT_EQ(15, to[1]);
    ASSERT_EQ(16, to[2]);

    ASSERT_TRUE(ConsistentPrecedenceConstraints(&pre));
}

TEST(Precedence, Regular2DGrid45DegreePrecedenceOneWide)
{
    Regular2DGrid45DegreePrecedence pre(1, 6);
    ASSERT_EQ(6, pre.NumBlocks());

    ASSERT_EQ(5, pre.NumPrecedenceConstraints());

    std::vector<IndexType> to;
    pre.AntecedentsVector(0, &to);
    ASSERT_EQ(1, to.size());
    ASSERT_EQ(1, to[0]);
}

TEST(Precedence, Regular2DGrid45DegreePrecedenceReachableAntecedents)
{
    Regular2DGrid45DegreePrecedence pre(10, 6);
    ASSERT_EQ(60, pre.NumBlocks());

    std::vector<int> expected {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 1, 1, 0, 0, 0,
        0, 0, 0, 1, 1, 1, 1, 1, 0, 0,
        0, 0, 1, 1, 1, 1, 1, 1, 1, 0,
        0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    }; // flipped

    std::vector<int> actual(60, 0);

    auto buffer = pre.GetNewSearchBuffer();
    for (auto & v : pre.ReachableAntecedents(5, buffer.get())) {
        actual[v] = 1;
    }

    ASSERT_VEC_EQ(expected, actual);
}

TEST(Precedence, Regular2DGrid45DegreePrecedenceReachableSuccessors)
{
    Regular2DGrid45DegreePrecedence pre(10, 6);
    ASSERT_EQ(60, pre.NumBlocks());

    std::vector<int> expected {
        0, 1, 1, 1, 2, 2, 1, 0, 0, 0,
        0, 0, 1, 1, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    }; // flipped

    std::vector<int> actual(60, 0);

    auto buffer = pre.GetNewSearchBuffer();
    for (auto & v : pre.ReachableSuccessors(15, buffer.get())) {
        actual[v]++;
    }
    for (auto & v : pre.ReachableSuccessors(23, buffer.get())) {
        actual[v]++;
    }

    ASSERT_VEC_EQ(expected, actual);
}

TEST(Precedence, Regular2DGrid45DegreePrecedenceAllConstraints)
{
    Regular2DGrid45DegreePrecedence pre(10, 6);

    IndexType count = 0;
    for (auto & v : pre.PrecedenceConstraints()) {
        ASSERT_TRUE(v.From < v.To);
        count++;
    }
    ASSERT_EQ(140, count);
}

TEST(Precedence, SlopePairLessthan)
{
    AzmSlopePair a( 90_deg, 45_deg);
    AzmSlopePair b(  0_deg, 50_deg);
    AzmSlopePair c(180_deg, 40_deg);

    ASSERT_TRUE(a < c);
    ASSERT_TRUE(b < c);
    ASSERT_TRUE(b < a);
}

TEST(Precedence, SlopeGetSimple)
{
    SlopeDefinition def({
        {static_cast<double>(  0_deg), static_cast<double>(45_deg)}, 
        {static_cast<double>( 90_deg), static_cast<double>(50_deg)}, 
        {static_cast<double>(180_deg), static_cast<double>(40_deg)}
    });

    ASSERT_NEAR(45.0_deg, def.Get(  0_deg), 0.0000001);
    ASSERT_NEAR(47.5_deg, def.Get( 45_deg), 0.0000001);
    ASSERT_NEAR(50.0_deg, def.Get( 90_deg), 0.0000001);
    ASSERT_NEAR(45.0_deg, def.Get(135_deg), 0.0000001);
    ASSERT_NEAR(40.0_deg, def.Get(180_deg), 0.0000001);
    ASSERT_NEAR(42.5_deg, def.Get(270_deg), 0.0000001);
    ASSERT_NEAR(45.0_deg, def.Get(360_deg), 0.0000001);
}

TEST(Precedence, SlopeGetSingle)
{
    SlopeDefinition def({
        {static_cast<double>(20_deg), static_cast<double>(45_deg)}
    });

    ASSERT_NEAR(45_deg, def.Get(  0_deg), 0.0000001);
    ASSERT_NEAR(45_deg, def.Get( 45_deg), 0.0000001);
    ASSERT_NEAR(45_deg, def.Get(270_deg), 0.0000001);
    ASSERT_NEAR(45_deg, def.Get(730_deg), 0.0000001);
    ASSERT_NEAR(45_deg, def.Get(-30_deg), 0.0000001);
    ASSERT_NEAR(45_deg, def.Get(999_deg), 0.0000001);
}

TEST(Precedence, SlopeGetRound)
{
    SlopeDefinition def({
        {static_cast<double>(  0_deg), static_cast<double>(40_deg)}, 
        {static_cast<double>(180_deg), static_cast<double>(50_deg)}
    });

    ASSERT_NEAR(40.0_deg, def.Get(  0_deg), 0.0000001);
    ASSERT_NEAR(42.5_deg, def.Get( 45_deg), 0.0000001);
    ASSERT_NEAR(45.0_deg, def.Get( 90_deg), 0.0000001);
    ASSERT_NEAR(47.5_deg, def.Get(135_deg), 0.0000001);
    ASSERT_NEAR(50.0_deg, def.Get(180_deg), 0.0000001);
    ASSERT_NEAR(45.0_deg, def.Get(270_deg), 0.0000001);
}

TEST(Precedence, SlopeCubic)
{
    SlopeDefinition def({
        {static_cast<double>(  0_deg), static_cast<double>(45_deg)}, 
        {static_cast<double>( 45_deg), static_cast<double>(45_deg)}, 
        {static_cast<double>( 90_deg), static_cast<double>(30_deg)}, 
        {static_cast<double>(135_deg), static_cast<double>(40_deg)}, 
        {static_cast<double>(180_deg), static_cast<double>(45_deg)}, 
        {static_cast<double>(270_deg), static_cast<double>(45_deg)}
    });

    SlopeDefinition def2 = CubicInterpolate(def, 512);
    ASSERT_EQ(512, def2.NumPairs());

    ASSERT_NEAR(45.0000_deg, def2.Get(  0_deg), 0.000001);
    ASSERT_NEAR(43.1476_deg, def2.Get(150_deg), 0.000001);
}

TEST(Precedence, SlopeCosine)
{
    SlopeDefinition def({
        {static_cast<double>(  0_deg), static_cast<double>(45_deg)}, 
        {static_cast<double>( 45_deg), static_cast<double>(45_deg)}, 
        {static_cast<double>( 90_deg), static_cast<double>(30_deg)}, 
        {static_cast<double>(135_deg), static_cast<double>(40_deg)}, 
        {static_cast<double>(180_deg), static_cast<double>(45_deg)}, 
        {static_cast<double>(270_deg), static_cast<double>(45_deg)}
    });

    SlopeDefinition def2 = CosineInterpolate(def, 512);
    ASSERT_EQ(512, def2.NumPairs());

    ASSERT_NEAR(45.0000_deg, def2.Get(  0_deg), 0.000001);
    ASSERT_NEAR(41.2503_deg, def2.Get(150_deg), 0.000001);
}

TEST(Precedence, SlopeViolateBase)
{
    SlopeDefinition def = SlopeDefinition::Constant(45_deg);

    ASSERT_EQ(true, def.Within( 1, 0, 1));
    ASSERT_EQ(true, def.Within(-1, 0, 1));
    ASSERT_EQ(true, def.Within( 0,  1, 1));
    ASSERT_EQ(true, def.Within( 0, -1, 1));
    ASSERT_EQ(false, def.Within( 1, 1, 1));
    ASSERT_EQ(true, def.Within( 2, 0, 2));
    ASSERT_EQ(true, def.Within( 4, 0, 4));
    ASSERT_EQ(true, def.Within( 2, 2, 4));
}

TEST(Precedence, SlopeViolateDual)
{
    SlopeDefinition def({
        {static_cast<double>(  0_deg), static_cast<double>(40_deg)}, 
        {static_cast<double>(180_deg), static_cast<double>(60_deg)}
    });

    ASSERT_EQ( true, def.Within(0,  1, 1));
    ASSERT_EQ(false, def.Within(0,  2, 1));
    ASSERT_EQ(false, def.Within(0, -1, 1));
}

TEST(Precedence, PatternOneFive)
{
    PrecedencePattern ptrn = PrecedencePattern::OneFive();
    ASSERT_EQ(5, ptrn.Offsets.size());

    ASSERT_EQ( 0, ptrn.Offsets[0].x);
    ASSERT_EQ(-1, ptrn.Offsets[0].y);
    ASSERT_EQ( 1, ptrn.Offsets[0].z);

    ASSERT_EQ(-1, ptrn.Offsets[1].x);
    ASSERT_EQ( 0, ptrn.Offsets[1].y);
    ASSERT_EQ( 1, ptrn.Offsets[1].z);

    ASSERT_EQ( 0, ptrn.Offsets[2].x);
    ASSERT_EQ( 0, ptrn.Offsets[2].y);
    ASSERT_EQ( 1, ptrn.Offsets[2].z);

    ASSERT_EQ( 1, ptrn.Offsets[3].x);
    ASSERT_EQ( 0, ptrn.Offsets[3].y);
    ASSERT_EQ( 1, ptrn.Offsets[3].z);

    ASSERT_EQ( 0, ptrn.Offsets[4].x);
    ASSERT_EQ( 1, ptrn.Offsets[4].y);
    ASSERT_EQ( 1, ptrn.Offsets[4].z);
}

TEST(Precedence, PatternOneNine)
{
    PrecedencePattern ptrn = PrecedencePattern::OneNine();
    ASSERT_EQ(9, ptrn.Offsets.size());

    int k = 0;
    for (int j = -1; j <= 1; j++) {
        for (int i = -1; i <= 1; i++) {
            ASSERT_EQ(i, ptrn.Offsets[k].x);
            ASSERT_EQ(j, ptrn.Offsets[k].y);
            ASSERT_EQ(1, ptrn.Offsets[k].z);
            k++;
        }
    }
}

TEST(Precedence, PatternMinSearch)
{
    PrecedencePattern ptrn = PrecedencePattern::MinSearch(45_deg, 10);
    ASSERT_EQ(25, ptrn.size());
}

static IBlockValuesSPtr ValuesFromVec(const std::vector<int64_t>& vs)
{
#ifdef MVD_USE_GMP
    auto values = std::make_shared<GMPBlockValues>(vs.size());
#else
    auto values = std::make_shared<VecBlockValues>(vs.size());
#endif

    for (IndexType i = 0; i < vs.size(); i++) {
        values->SetBlockValueSI(i, vs[i]);
    }
    return values;
}

TEST(MFlow, LargestMinCutTiny)
{
    std::vector<int64_t> v {7, 2, -2, -2, -2};
    auto values = ValuesFromVec(v);
    auto values2 = std::make_shared<SolveLargestValuesAdapter>(values);
    auto pre = std::make_shared<ExplicitPrecedence>(
            values->NumBlocks(), 
            std::initializer_list<std::initializer_list<int>>
                {{0, 2}, {0, 3}, {1, 3}, {1, 4}});

    PseudoSolverSolveInfo info, info2;
    PseudoSolver solver(pre, values);

    solver.Solve(&info);
    ASSERT_EQ(3, info.NumContainedNodes);

    PseudoSolver solver2(pre, values2);
    solver2.Solve(&info2);

    ASSERT_EQ(5, info2.NumContainedNodes);
}

TEST(MFlow, LargestMinCutMMW)
{
    std::vector<int64_t> v {3, 8, 1, -2, -2, -2, -2, 0, 0, 0, 0, 0};
    auto values = ValuesFromVec(v);
    auto pre = std::make_shared<ExplicitPrecedence>(
            values->NumBlocks(), 
            std::initializer_list<std::initializer_list<int>>
                {{0, 3}, {0, 4}, {1, 4}, {1, 5}, {2, 5}, {2, 6},
                {7, 0}, {7, 1}, {8, 1}, {8, 2}, {9, 3}, {9, 4},
                {10, 4}, {10, 5}, {11, 5}, {11, 6}});

    auto values2 = std::make_shared<SolveLargestValuesAdapter>(values);
    PseudoSolverSolveInfo info, info2;
    PseudoSolver solver(pre, values);

    solver.Solve(&info);
    ASSERT_EQ(5, info.NumContainedNodes);

    PseudoSolver solver2(pre, values2);
    solver2.Solve(&info2);
    ASSERT_EQ(8, info2.NumContainedNodes);
}

IBlockValuesSPtr ReadTestDataValues(const std::string& stem)
{
    std::ostringstream os;
    os << "../data/" << stem << ".dat";

    std::ifstream input(os.str());
    std::string line;

    std::vector<int64_t> values;
    while (std::getline(input, line)) {
        int64_t v = static_cast<int64_t>(std::stoi(line));
        values.push_back(v);
    }
    return ValuesFromVec(values);
}

TEST(MFlow, Sim2D76)
{
    BlockDefinition bdef = BlockDefinition::UnitModel(75, 1, 40);
    auto values = ReadTestDataValues("sim2d76");
    auto pre = std::make_shared<Regular2DGrid45DegreePrecedence>(bdef.NumX, bdef.NumZ);

    ASSERT_TRUE(values);
    ASSERT_EQ(values->NumBlocks(), bdef.NumBlocks());
    ASSERT_EQ(values->NumBlocks(), pre->NumBlocks());

    PseudoSolverSolveInfo info;

    PseudoSolver solver(pre, values.get());
    solver.Solve(&info);

    ASSERT_EQ(info.NumNodes, values->NumBlocks());
    ASSERT_EQ(info.NumContainedNodes, 945);
    ASSERT_VALUETYPE_EQ_INT(info.ContainedValue, 295932);

    m_ExtraElapsed = info.ElapsedSeconds;
}

TEST(MFlow, Sim2D76Largest)
{
    BlockDefinition bdef = BlockDefinition::UnitModel(75, 1, 40);
    auto values = ReadTestDataValues("sim2d76");
    auto pre = std::make_shared<Regular2DGrid45DegreePrecedence>(bdef.NumX, bdef.NumZ);

    ASSERT_TRUE(values);
    ASSERT_EQ(values->NumBlocks(), bdef.NumBlocks());
    ASSERT_EQ(values->NumBlocks(), pre->NumBlocks());

    PseudoSolverSolveInfo info, info2;

    PseudoSolver solver(pre, values.get());
    solver.Solve(&info);

    solver.UpdateValues(std::make_shared<SolveLargestValuesAdapter>(values));

    solver.Solve(&info2);

    ASSERT_EQ(info.NumNodes, values->NumBlocks());
    ASSERT_EQ(info.NumContainedNodes, 945);
    ASSERT_VALUETYPE_EQ_INT(info.ContainedValue, 295932);

    ASSERT_EQ(info2.NumNodes, values->NumBlocks());
    ASSERT_EQ(info2.NumContainedNodes, 946);
//    ASSERT_VALUETYPE_EQ_INT(info2.ContainedValue, 295932); This is of course not
//    true because we use the value based adapter!
}

TEST(MFlow, BauxiteMed)
{
    BlockDefinition bdef = BlockDefinition::UnitModel(120, 120, 26);
    auto values = ReadTestDataValues("bauxitemed");
    auto pattern = PrecedencePattern::MinSearch(45_deg, 8);
    auto pre = std::make_shared<Regular3DBlockModelPatternPrecedence>(bdef, pattern);

    ASSERT_TRUE(values);
    ASSERT_EQ(values->NumBlocks(), bdef.NumBlocks());
    ASSERT_EQ(values->NumBlocks(), pre->NumBlocks());

    PseudoSolverSolveInfo info;

    PseudoSolver solver(pre, values.get());
    solver.Solve(&info);

    ASSERT_EQ(info.NumNodes, values->NumBlocks());
    ASSERT_EQ(info.NumContainedNodes, 74412);
    ASSERT_VALUETYPE_EQ_INT(info.ContainedValue, 28416592);
    m_ExtraElapsed = info.ElapsedSeconds;
}

TEST(MFlow, BauxiteMedLargest)
{
    BlockDefinition bdef = BlockDefinition::UnitModel(120, 120, 26);
    auto values = ReadTestDataValues("bauxitemed");
    auto pattern = PrecedencePattern::MinSearch(45_deg, 8);
    auto pre = std::make_shared<Regular3DBlockModelPatternPrecedence>(bdef, pattern);

    ASSERT_TRUE(values);
    ASSERT_EQ(values->NumBlocks(), bdef.NumBlocks());
    ASSERT_EQ(values->NumBlocks(), pre->NumBlocks());

    PseudoSolverSolveInfo info;

    PseudoSolver solver(pre, values.get());
    solver.SolveLargest(&info);

    ASSERT_EQ(info.NumNodes, values->NumBlocks());
    ASSERT_EQ(info.NumContainedNodes, 76813);
    ASSERT_VALUETYPE_EQ_INT(info.ContainedValue, 28416592);
    m_ExtraElapsed = info.ElapsedSeconds;
}

TEST(MFlow, CuCase)
{
    BlockDefinition bdef = BlockDefinition::UnitModel(170, 215, 50);
    auto values = ReadTestDataValues("cucase");
    auto pattern = PrecedencePattern::MinSearch(45_deg, 8);
    auto pre = std::make_shared<Regular3DBlockModelPatternPrecedence>(bdef, pattern);

    ASSERT_TRUE(values);
    ASSERT_EQ(values->NumBlocks(), bdef.NumBlocks());
    ASSERT_EQ(values->NumBlocks(), pre->NumBlocks());

    PseudoSolverSolveInfo info;

    PseudoSolver solver(pre, values.get());
    solver.Solve(&info);

    ASSERT_EQ(info.NumNodes, values->NumBlocks());
    ASSERT_EQ(info.NumContainedNodes, 357304);
    ASSERT_VALUETYPE_EQ_INT(info.ContainedValue, 19175685);
    m_ExtraElapsed = info.ElapsedSeconds;
}

TEST(MFlow, CuPipe)
{
    BlockDefinition bdef = BlockDefinition::UnitModel(180, 180, 85);
    auto values = ReadTestDataValues("cupipe");
    auto pattern = PrecedencePattern::MinSearch(45_deg, 8);
    auto pre = std::make_shared<Regular3DBlockModelPatternPrecedence>(bdef, pattern);

    ASSERT_TRUE(values);
    ASSERT_EQ(values->NumBlocks(), bdef.NumBlocks());
    ASSERT_EQ(values->NumBlocks(), pre->NumBlocks());

    PseudoSolverSolveInfo info;

    PseudoSolver solver(pre, values.get());
    solver.Solve(&info);

    ASSERT_EQ(info.NumNodes, values->NumBlocks());
    ASSERT_EQ(info.NumContainedNodes, 198078);
    ASSERT_VALUETYPE_EQ_INT(info.ContainedValue, 102306787);
    m_ExtraElapsed = info.ElapsedSeconds;
}

TEST(MFlow, McLaughlinGeo)
{
    BlockDefinition bdef = BlockDefinition::UnitModel(140, 296, 68);
    auto values = ReadTestDataValues("mclaughlingeo");
    auto pattern = PrecedencePattern::MinSearch(45_deg, 8);
    auto pre = std::make_shared<Regular3DBlockModelPatternPrecedence>(bdef, pattern);

    ASSERT_TRUE(values);
    ASSERT_EQ(values->NumBlocks(), bdef.NumBlocks());
    ASSERT_EQ(values->NumBlocks(), pre->NumBlocks());

    PseudoSolverSolveInfo info;

    PseudoSolver solver(pre, values.get());
    solver.Solve(&info);

    ASSERT_EQ(info.NumNodes, values->NumBlocks());
    ASSERT_EQ(info.NumContainedNodes, 345936);
    ASSERT_VALUETYPE_EQ_INT(info.ContainedValue, 1145395060);
    m_ExtraElapsed = info.ElapsedSeconds;
}

//TEST(MFlow, BigGold)
//{
//    BlockDefinition bdef = BlockDefinition::UnitModel(483, 333, 101);
//    auto values = ReadTestDataValues("biggold");
//    auto pattern = PrecedencePattern::MinSearch(45_deg, 8);
//    auto pre = std::make_shared<Regular3DBlockModelPatternPrecedence>(bdef, pattern);
//
//    ASSERT_TRUE(values);
//    ASSERT_EQ(values->NumBlocks(), bdef.NumBlocks());
//    ASSERT_EQ(values->NumBlocks(), pre->NumBlocks());
//
//    PseudoSolverSolveInfo info;
//
//    PseudoSolver solver(pre, values.get());
//    solver.Solve(&info);
//
//    ASSERT_EQ(info.NumNodes, values->NumBlocks());
//    ASSERT_EQ(info.NumContainedNodes, 602150);
//    ASSERT_VALUETYPE_EQ_INT(info.ContainedValue, 23734996);
//    m_ExtraElapsed = info.ElapsedSeconds;
//}

#ifndef MVD_USE_GMP
TEST(README, Explicit)
{
    int64_t numBlocks = 5;
    auto values = std::make_shared<VecBlockValues>(numBlocks);
    values->SetBlockValueSI(0,  7); // 0 being the index, 7 the economic block value
    values->SetBlockValueSI(1,  2);
    values->SetBlockValueSI(2, -2);
    values->SetBlockValueSI(3, -2);
    values->SetBlockValueSI(4, -2);

    auto precedence = std::make_shared<ExplicitPrecedence>(numBlocks); // Again 5 being the number of blocks
    precedence->AddPrecedenceConstraint(0, 2);
    precedence->AddPrecedenceConstraint(0, 3);
    precedence->AddPrecedenceConstraint(1, 3);
    precedence->AddPrecedenceConstraint(1, 4);

    PseudoSolver solver(precedence, values);
    solver.Solve();

    ASSERT_EQ(solver.InMinimumCut(0), true);
    ASSERT_EQ(solver.InMinimumCut(1), false);
    ASSERT_EQ(solver.InMinimumCut(2), true);
    ASSERT_EQ(solver.InMinimumCut(3), true);
    ASSERT_EQ(solver.InMinimumCut(4), false);
}
#endif

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    auto& tests = RegisteredTests();

    // TODO: Filter tests based on argc / argv?
    // TODO: Handle alternative path to test data?
    
    // Running tests
    std::cout << "\033[0;32m" << "[=== 0_0 ===]" << "\033[0m" << " ";
    std::cout << tests.size() << " tests" << std::endl;
    std::vector<RegisteredTest> failedTests;
    for (auto & test : tests) {
        std::cout << "\033[0;32m" << "[ RUN       ]" << "\033[0m" << " " << test.ClassName << "." << test.InstanceName << std::endl;
        bool failed = false;
        auto start = std::chrono::steady_clock::now();
        try {
            test.Test->DoTest();
        } catch (TestException& e) {
            failed = true;
            std::cout << e.GetTestFailedMessage() << std::endl;
        }
        auto end = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        if (!failed) {
            std::cout << "\033[0;32m" << "[        OK ]" << "\033[0m" << " " << test.ClassName << "." << test.InstanceName << " (" << elapsed << " ms)";
            if (test.Test->m_ExtraElapsed >= 0) {
                std::cout << " [" << std::fixed << std::setprecision(2) << test.Test->m_ExtraElapsed << "s]";
            }
            std::cout << std::endl;
        } else {
            std::cout << "\033[0;31m" << "[   FAILED  ]" << "\033[0m" << " " << test.ClassName << "." << test.InstanceName << std::endl;
        }
    }
    std::cout << "\033[0;32m" << "[===========]" << "\033[0m" << std::endl;

    return 0;
}
#endif // MVD_MINEFLOW_TESTS

////////////////////////////////////////////////////////////////////////////////

#if (defined(MVD_MINEFLOW_EXE) && !defined(MVD_MINEFLOW_TESTS))
#include <fstream>

using mclock = std::chrono::steady_clock;

bool ArgPeekString(int* argc, char*** argv, std::string* str)
{
    if (*argc <= 0) return false;
    *str = (*argv)[0];
    return true;
}

void ArgNext(int* argc, char*** argv)
{
    (*argc)--;
    (*argv)++;
}

bool ArgReadString(int* argc, char*** argv, std::string* str)
{
    if (!ArgPeekString(argc, argv, str)) return false;
    ArgNext(argc, argv);
    return true;
}

bool ArgReadInt64(int* argc, char*** argv, int64_t* v)
{
    std::string arg;
    if (!ArgReadString(argc, argv, &arg)) return false;
    std::istringstream is(arg);
    int64_t q;
    is >> q;
    if (is.fail()) {
        return false;
    }
    *v = q;
    return true;
}

bool ArgReadFloat(int* argc, char*** argv, float* v)
{
    std::string arg;
    if (!ArgReadString(argc, argv, &arg)) return false;

    std::istringstream is(arg);
    float q;
    is >> q;
    if (is.fail()) {
        return false;
    }
    *v = q;
    return true;
}

std::string Elapsed(mclock::time_point start, mclock::time_point end)
{
    int milliseconds = 
        static_cast<int>(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    int seconds = milliseconds / 1000;
    int minutes = seconds / 60;
    int hours = minutes / 60;

    milliseconds %= 1000;
    seconds %= 60;;
    minutes %= 60;

    std::ostringstream os;

    os << std::setw(2) << std::setfill('0') << hours << ":"
       << std::setw(2) << std::setfill('0') << minutes << ":"
       << std::setw(2) << std::setfill('0') << seconds << "."
       << std::setw(3) << std::setfill('0') << milliseconds;
    return os.str();
}

////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<IPrecedenceConstraints> InitRegular(int* argc, char*** argv, IndexType* numBlocks)
{
    BlockDefinition blockDef = BlockDefinition::UnitModel();
    float slope = 0.0f;
    IndexType nBenches = 9; // could be input?
    if (!ArgReadInt64(argc, argv, &blockDef.NumX) ||
        !ArgReadInt64(argc, argv, &blockDef.NumY) ||
        !ArgReadInt64(argc, argv, &blockDef.NumZ) ||
        !ArgReadFloat(argc, argv, &slope)) {
        std::cerr << "Failed reading numx numy numz slope" << std::endl;
        return nullptr;
    }

    *numBlocks = blockDef.NumBlocks();

    SlopeDefinition slopeDef = SlopeDefinition::Constant(ToRadians(slope));
    PrecedencePattern pattern = PrecedencePattern::MinSearch(blockDef, slopeDef, nBenches);
    return std::make_shared<Regular3DBlockModelPatternPrecedence>(blockDef, pattern);
}

std::shared_ptr<IPrecedenceConstraints> InitMinSearch(int* argc, char*** argv, IndexType* numBlocks)
{
    const long double TAU = 6.283185307179586476925286766559;
    BlockDefinition blockDef = BlockDefinition::UnitModel();

    std::string minSearchFile;
    if (!ArgReadString(argc, argv, &minSearchFile)) {
        std::cerr << "Failed reading min search file argument" << std::endl;
        return nullptr;
    }

    std::ifstream in(minSearchFile);
    if (!in.good()) {
        std::cerr << "Failed opening min search file" << std::endl;
        return nullptr;
    }

    std::string line;
    {
        if (!std::getline(in, line)) {
            std::cerr << "Failed reading NumX NumY NumZ" << std::endl;
            return nullptr;
        }
        std::istringstream is(line);

        is >> blockDef.NumX >> blockDef.NumY >> blockDef.NumZ;
        if (is.fail()) {
            std::cerr << "Failed reading NumX NumY NumZ" << std::endl;
            return nullptr;
        }
        if (blockDef.NumX <= 0 ||
            blockDef.NumY <= 0 ||
            blockDef.NumZ <= 0) {
            std::cerr << "Invalid NumX NumY NumZ" << std::endl;
            return nullptr;
        }
    }

    {
        if (!std::getline(in, line)) {
            std::cerr << "Failed reading SizeX SizeY SizeZ" << std::endl;
            return nullptr;
        }
        std::istringstream is(line);

        is >> blockDef.SizeX >> blockDef.SizeY >> blockDef.SizeZ;
        if (is.fail()) {
            std::cerr << "Failed reading SizeX SizeY SizeZ" << std::endl;
            return nullptr;
        }
        if (blockDef.SizeX <= 0 ||
            blockDef.SizeY <= 0 ||
            blockDef.SizeZ <= 0) {
            std::cerr << "Invalid SizeX SizeY SizeZ" << std::endl;
            return nullptr;
        }
    }

    IndexType numBenches = 0;
    {
        if (!std::getline(in, line)) {
            std::cerr << "Failed reading numBenches" << std::endl;
            return nullptr;
        }
        std::istringstream is(line);

        int64_t nb;
        is >> nb;
        numBenches = static_cast<IndexType>(nb);
        if (is.fail()) {
            std::cerr << "Failed reading numBenches" << std::endl;
            return nullptr;
        }
        if (numBenches <= 1) {
            std::cerr << "Invalid num benches" << std::endl;
            return nullptr;
        }
    }

    std::vector<AzmSlopePair> pairs;
    while (std::getline(in, line)) {
        std::istringstream is(line);
        double azimuth, slope;
        is >> azimuth >> slope;
        if (is.fail()) {
            std::cerr << "Failed reading azimuth slope" << std::endl;
            return nullptr;
        }
        azimuth = ToRadians(azimuth);
        slope = ToRadians(slope);

        if (slope <= 0) {
            std::cerr << "Invalid slope" << std::endl;
            return nullptr;
        }

        while (azimuth >= TAU) azimuth -= TAU;
        while (azimuth < 0) azimuth += TAU;
        pairs.emplace_back(azimuth, slope);
    }

    if (pairs.empty()) {
        std::cerr << "Failed reading slope definition" << std::endl;
        return nullptr;
    }

    SlopeDefinition slopeDef(pairs);

    *numBlocks = blockDef.NumBlocks();

    PrecedencePattern pattern = PrecedencePattern::MinSearch(blockDef, slopeDef, numBenches);
    return std::make_shared<Regular3DBlockModelPatternPrecedence>(blockDef, pattern);
}

std::shared_ptr<IPrecedenceConstraints> InitExplicit(int* argc, char*** argv, IndexType* numBlocks)
{
    std::string explicitFile;
    if (!ArgReadString(argc, argv, &explicitFile)) {
        std::cerr << "Failed reading explicit precedence file argument" << std::endl;
        return nullptr;
    }

    std::ifstream in(explicitFile);
    if (!in.good()) {
        std::cerr << "Failed opening explicit precedence file" << std::endl;
        return nullptr;
    }

    std::string line;
    if (!std::getline(in, line)) {
        std::cerr << "Failed reading num blocks line" << std::endl;
        return nullptr;
    }
    {
        std::istringstream is(line);
        is >> *numBlocks;
        if (is.fail()) {
            std::cerr << "Failed reading num blocks" << std::endl;
            return nullptr;
        }
        if (*numBlocks <= 0) {
            std::cerr << "Invalid num blocks" << std::endl;
            return nullptr;
        }
    }


    std::unordered_map<IndexType, std::vector<IndexType>> antecedents;
    while (std::getline(in, line)) {
        IndexType fromIndex;
        std::istringstream is(line);
        is >> fromIndex;
        if (is.fail()) {
            std::cerr << "Failed reading from index" << std::endl;
        }

        if (fromIndex < 0 || fromIndex >= *numBlocks) {
            std::cerr << "Invalid block index in precedence file: " << fromIndex << std::endl;
            return nullptr;
        }

        std::vector<IndexType>& toBlocks = antecedents[fromIndex];
        IndexType toIndex;
        while (is >> toIndex) {
            toBlocks.push_back(toIndex);

            if (toIndex < 0 || toIndex >= *numBlocks) {
                std::cerr << "Invalid block index in precedence file: " << toIndex << std::endl;
                return nullptr;
            }
        }
    }

    return std::make_shared<ExplicitPrecedence>(*numBlocks, std::move(antecedents));
}

bool InitValues(const std::string& valuesPath, IndexType numBlocks, std::shared_ptr<IBlockValues>* values)
{
    std::ifstream in(valuesPath);
    if (!in.good()) {
        std::cerr << "Failed opening values file" << std::endl;
        return false;
    }

    std::vector<ValueType> vs(numBlocks);
    for (size_t i = 0; i < static_cast<size_t>(numBlocks); i++) {
        int64_t v;
        in >> v;
        vs[i] = v;
        if (in.fail()) {
            std::cerr << "Failed reading values line: " << i + 1 << std::endl;
            return false;
        }
    }

    *values = std::make_shared<VecBlockValues>(std::move(vs));

    return true;
}

int main(int argc, char** argv)
{
    ArgNext(&argc, &argv);
    if (argc == 0) {
        std::cerr << "Usage: mineflow [options] data.dat output.dat" << std::endl;
        std::cerr << "options:" << std::endl;
        std::cerr << " --regular <nx> <ny> <nz> <slope> Use a single constant slope angle (deg)" << std::endl;
        std::cerr << " --minsearch <file>               Use a single minimum search pattern" << std::endl;
        std::cerr << " --explicit <file>                Use explicit precedence constants (slow!)" << std::endl;
        std::cerr << " --to_dimacs                      Outputs a dimacs file to stdout then exits" << std::endl;
        std::cerr << std::endl;
        std::cerr << "minsearch format:" << std::endl;
        std::cerr << "<NumX> <NumY> <NumZ>      # Number of blocks in x, y, and z" << std::endl;
        std::cerr << "<SizeX> <SizeY> <SizeZ>   # Size of blocks in x, y, and z" << std::endl;
        std::cerr << "<NumBenches>              # Number of benches to extent pattern" << std::endl;
        std::cerr << "<Azimuth> <Slope>" << std::endl;
        std::cerr << "<Azimuth> <Slope>" << std::endl;
        std::cerr << "..." << std::endl;
        std::cerr << std::endl;
        std::cerr << "explicit format:" << std::endl;
        std::cerr << "<num blocks>" << std::endl;
        std::cerr << "<from_block_id> <to_block_id_0> <to_block_id_1> ..." << std::endl;
        std::cerr << "<from_block_id> <to_block_id_0> <to_block_id_1> ..." << std::endl;
        std::cerr << "..." << std::endl;
        std::cerr << std::endl;
        std::cerr << "'data.dat' format:" << std::endl;
        std::cerr << "<value_block_0>" << std::endl;
        std::cerr << "<value_block_1>" << std::endl;
        std::cerr << "..." << std::endl;
        std::cerr << std::endl;
        std::cerr << "'output.dat' format:" << std::endl;
        std::cerr << "<mine_block_0>" << std::endl;
        std::cerr << "<mine_block_1>" << std::endl;
        std::cerr << "..." << std::endl;
        return 1;
    }

    mclock::time_point programStart = mclock::now();

    
    std::shared_ptr<IPrecedenceConstraints> pre = nullptr;
    std::shared_ptr<IBlockValues> values = nullptr;
    IndexType numBlocks = 0;
    bool outputtingDimacs = false;
    float multiplier = 100.0f;

    std::string argument;
    while (ArgPeekString(&argc, &argv, &argument)) {
        if (argument.rfind("--", 0) == 0) {
            ArgNext(&argc, &argv);
            if (argument == "--to_dimacs") {
                outputtingDimacs = true;
            } else if (argument == "--regular") {
                pre = InitRegular(&argc, &argv, &numBlocks);
            } else if (argument == "--minsearch") {
                pre = InitMinSearch(&argc, &argv, &numBlocks);
            } else if (argument == "--explicit") {
                pre = InitExplicit(&argc, &argv, &numBlocks);
            } else {
                std::cerr << "Unknown argument: " << argument << std::endl;
                return 1;
            }
        } else {
            break;
        }
    }

    if (!pre || numBlocks <= 0) {
        std::cerr << "No precedence specified, or no blocks in input" << std::endl;
        return 1;
    }


    if (!ArgReadString(&argc, &argv, &argument)) {
        std::cerr << "No data file argument" << std::endl;
        return 1;
    } else if (!InitValues(argument, numBlocks, &values)) {
        std::cerr << "failure initializing values" << argument << std::endl;
        return 1;
    }

    if (outputtingDimacs) {
// TODO
//    uint64_t numBlocks = values->NumBlocks();
//    uint64_t numPrecedence = pre->NumPrecedenceConstraints();
//
//    uint64_t numNodes = numBlocks + 2;
//    uint64_t numArcs = numPrecedence + numBlocks;
//
//    uint64_t sourceIndex = numBlocks;
//    uint64_t sinkIndex = numBlocks + 1;
//
//    std::ofstream of(name);
//    of << "p max " << numNodes << " " << numArcs << std::endl;
//    of << "n " << sourceIndex + 1 << " s" << std::endl;
//    of << "n " << sinkIndex + 1 << " t" << std::endl;
//    double maxval = 0;
//    for (uint64_t blockIndex = 0; blockIndex < numBlocks; blockIndex++) {
//        double v = values->Value(blockIndex);
//        if (v > 0) {
//            of << "a " << sourceIndex + 1 << " " << blockIndex + 1 << " " << 
//                static_cast<int>(v) << std::endl;
//            maxval += v;
//        } else {
//            of << "a " << blockIndex + 1 << " " << sinkIndex + 1 << " " << 
//                static_cast<int>(-v) << std::endl;
//        }
//    }
//
//    for (uint64_t blockIndex = 0; blockIndex < numBlocks; blockIndex++) {
//        pre->Antecedents(blockIndex, [&](uint64_t targetIndex){
//            of << "a " << blockIndex + 1 << " " << targetIndex + 1 << " " << 
//                static_cast<int>(maxval) << std::endl;
//            return true;
//        });
//    }
    } else {
        if (!ArgReadString(&argc, &argv, &argument)) {
            std::cerr << "No output file argument" << std::endl;
            return 1;
        }

        mclock::time_point readInput = mclock::now();
        std::cout << "  MineFlow - Version 1.0" << std::endl;
        std::cout << "--------------------------" << std::endl;
        std::cout << "Num blocks  : " << numBlocks << std::endl;

        PseudoSolver solver(pre, values.get());
        mclock::time_point initialized = mclock::now();

        PseudoSolverSolveInfo info;
        solver.Solve(&info);
        mclock::time_point solved = mclock::now();

        std::ofstream of(argument);
        for (IndexType i = 0; i < info.NumNodes; i++) {
            if (solver.InMinimumCut(i)) {
                of << i << "\n";
            }
        }
        mclock::time_point output = mclock::now();

        std::cout << "Num mined   : " << info.NumContainedNodes << std::endl;
        std::cout << "Value       : " << info.ContainedValue << std::endl;
        std::cout << "--------------------------" << std::endl;
        std::cout << "Read data   : " << Elapsed(programStart, readInput) << std::endl;
        std::cout << "Init solver : " << Elapsed(readInput, initialized) << std::endl;
        std::cout << "Solved      : " << Elapsed(initialized, solved) << std::endl;
        std::cout << "Saved       : " << Elapsed(solved, output) << std::endl;
        std::cout << "--------------------------" << std::endl;
        std::cout << "Total       : " << Elapsed(programStart, output) << std::endl;
    }
    return 0;
}

#endif // MVD_MINEFLOW_EXE

