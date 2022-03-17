/* mineflow.h
Contact: matthewvdeutsch@gmail.com
You should also have 'mineflow.cpp', or else you'll have trouble.
 
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

================================================================================


    +-------+                                               +-----+
    |        \                                             /      |
    |         \   _                                       /       |
    |          \_/ \                                 ____/        |
    |               \                               /             |
    |                \                             /              |
    |                 \            ____           /               |
    |                  \      ____/    \_________/                |
    |                   \____/                                    |
    |                                                             |
    +-------------------------------------------------------------+

These two files help answer:

    "Given an economic block model, which blocks should I mine?"

Conventionally called "The Ultimate Pit Problem", as introduced by Lerchs and
Grossmann in 1965, and solved in a wide variety of ways. These files implement
the "Pseudoflow" algoritihm from Hochbaum, modified to contend only with the
ultimate pit problem.

Also contains the 'minimum search patterns' from Caccetta and Giannini to
generate efficient (both geometrically, and memory wise) precedence schemes.
These files also includes a few other fancy novel things.

================================================================================

Build instructions
------------------
These two files should compile readily with any c++17 compliant compiler.
Contributions are welcome to fix that sort of thing!

Example commands:

To compile and run the tests:
    g++ -std=c++1z -O3 -DMVD_MINEFLOW_TESTS mineflow.cpp -o mineflow_tests
    ./mineflow_tests

To compile the standalone executable:
    g++ -std=c++1z -O3 -DMVD_MINEFLOW_EXE mineflow.cpp -o mineflow
    ./mineflow

To compile the static library:
    g++ -std=c++1z -O3 mineflow.cpp -c
    ar rvs mineflow.a mineflow.o

Also see the CMakeLists.txt


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

#ifndef MVD_INCLUDE_MINEFLOW_H
#define MVD_INCLUDE_MINEFLOW_H

#include <algorithm>
#include <iterator>
#include <memory>
#include <cmath>
#include <queue>
#include <vector>
#include <iostream>
#include <cstdint>

#ifdef MVD_USE_GMP
#include "gmp.h"
#endif

namespace mvd::mineflow {

////////////////////////////////////////////////////////////////////////////////

#ifdef MVD_USE_GMP
typedef mpz_t ValueType; // careful
#else
typedef int64_t ValueType;
#endif

typedef int64_t IndexType; // Must be signed.
struct BlockDefinition; // 3d regular block model

// General Usage:
//    - Define a class that inherits from IBlockValues and implements:
//          IBlockValues::NumBlocks, and
//          IBlockValues::BlockValue
//              (look at VecBlockValues as an example)
//
//    - Define a class that inherits from IPrecedenceConstraints and implements:
//          IPrecedenceConstraints::NumBlocks, and
//          IPrecedenceConstraints::Antecedents
//              (look at Regular2DGrid45DegreePrecedence as an example)
//
//    - Then use: mineflow::PseudoSolver directly to solve.
//
// You can search in mineflow.cpp for 'MVD_MINEFLOW_TESTS_BEGIN' for examples!


////////////////////////////////////////////////////////////////////////////////
// The fundamental interfaces / operations
////////////////////////////////////////////////////////////////////////////////
class IBlockValues 
{
public:
    IBlockValues();
    virtual ~IBlockValues();

    virtual IndexType NumBlocks() const = 0;
    virtual void BlockValue(IndexType blockIndex, ValueType* value) const = 0;
};
typedef std::shared_ptr<IBlockValues> IBlockValuesSPtr;

class BlockIndexInputIteratorBase;
class PrecedenceConstraintInputIteratorBase;
class PrecedenceConstraintsReachableSearchBuffer;
typedef std::unique_ptr<PrecedenceConstraintsReachableSearchBuffer> 
    PrecedenceConstraintsReachableSearchBufferPtr;

// The block at 'From', requires the block at 'To' to be mined first
struct PrecedenceConstraint
{
    IndexType From, To;
};

// The main precedence constraints interface
class IPrecedenceConstraints
{
public:
    IPrecedenceConstraints();
    virtual ~IPrecedenceConstraints();

    virtual IndexType NumBlocks() const = 0;

    // Input iterators are lightweight and /single/ use.
    // Example usage is:
    // for (auto to : pre->Antecedents(from)) {
    //     // do something
    // }
    virtual BlockIndexInputIteratorBase Antecedents(IndexType fromBlockIndex) const = 0;
    virtual BlockIndexInputIteratorBase Successors(IndexType toBlockIndex) const;

    // If you need to store the information use the following info to help: 
    // Counting (generally) requires iterating, so it should be avoided
    virtual IndexType NumAntecedents(IndexType fromBlockIndex) const; // May be expensive
    virtual IndexType ApproxNumAntecedents(IndexType fromBlockIndex) const; // May return 0
    virtual void AntecedentsVector(IndexType fromBlockIndex, std::vector<IndexType>* vec) const;

    virtual IndexType NumSuccessors(IndexType toBlockIndex) const; // May be expensive
    virtual IndexType ApproxNumSuccessors(IndexType toBlockIndex) const; // May return 0
    virtual void SuccessorsVector(IndexType toBlockIndex, std::vector<IndexType>* vec) const;

    // Sometimes we want all the precedence constraints.
    // Example usage is:
    // for (auto [from, to] : pre->PrecedenceConstraints()) {
    //     // do something
    // }
    virtual IndexType NumPrecedenceConstraints() const; // May be expensive!!
    virtual IndexType ApproxNumPrecedenceConstraints() const; // May return 0
    virtual PrecedenceConstraintInputIteratorBase PrecedenceConstraints() const;
    virtual void PrecedenceConstraintsVector(std::vector<PrecedenceConstraint>* vec) const;

    // When doing reachable checks provide a re-useable search buffer
    // Example usage is:
    // auto buffer = pre->GetNewSearchBuffer();
    // for (auto reachable : pre->ReachableAntecedents(from, buffer.get())) {
    //     // do something
    // }
    virtual PrecedenceConstraintsReachableSearchBufferPtr GetNewSearchBuffer() const;
    virtual BlockIndexInputIteratorBase ReachableAntecedents(IndexType fromBlockIndex, 
            PrecedenceConstraintsReachableSearchBuffer* buffer) const;
    virtual BlockIndexInputIteratorBase ReachableSuccessors(IndexType toBlockIndex, 
            PrecedenceConstraintsReachableSearchBuffer* buffer) const;

    // Does an check allowing for 'partial' searching 
    // return true on continue from that block and false to not.
    virtual void PartialReachableAntecedents(IndexType fromBlockIndex,
            std::function<bool(IndexType toBlockIndex)> cback,
            PrecedenceConstraintsReachableSearchBuffer* buffer) const;
    virtual void PartialReachableSuccessors(IndexType toBlockIndex,
            std::function<bool(IndexType fromBlockIndex)> cback,
            PrecedenceConstraintsReachableSearchBuffer* buffer) const;
};
typedef std::shared_ptr<IPrecedenceConstraints> IPrecedenceConstraintsSPtr;


// And the solver
namespace impl {
struct Arc;
struct Node;
class NodePool;
class PrecedenceArcPool;
}

struct PseudoSolverSolveInfo
{
    PseudoSolverSolveInfo();
    ~PseudoSolverSolveInfo();

    double ElapsedSeconds;
    IndexType NumNodes;
    IndexType NumContainedNodes;
    IndexType NumUsedPrecedenceConstraints;

    ValueType ContainedValue;
};
std::string PseudoSolverSolveInfoToString(const PseudoSolverSolveInfo& info);
inline std::ostream& operator<<(std::ostream& os, const PseudoSolverSolveInfo& info)
{
    os << PseudoSolverSolveInfoToString(info);
    return os;
}

class PseudoSolver
{
public:
    PseudoSolver(
        std::shared_ptr<const IPrecedenceConstraints> pre, // retained
        const IBlockValues* values = nullptr // read once to init structure
    );
    PseudoSolver(std::shared_ptr<const IPrecedenceConstraints> pre,
        std::shared_ptr<const IBlockValues> values);
    ~PseudoSolver();

    IndexType NumNodes() const;

    void Solve(PseudoSolverSolveInfo* info = nullptr);
    bool InMinimumCut(IndexType nodeIndex) const;

    // I discourage these two in favor of the SolveLargestValuesAdapter
    // but that overflows if you're not using gmp, so welp:
    void SolveLargest(PseudoSolverSolveInfo* info = nullptr);
    bool InLargestMinimumCut(IndexType nodeIndex) const;

    void UpdateValues(const IBlockValues* values); // you must then call 'Solve' again.
    void UpdateValues(std::shared_ptr<const IBlockValues> values); // for convenience

private:
    void ProcessStrongRoot(impl::Node* strongRoot);
    void ProcessChildren(impl::Node* node);

    void Merge(impl::Node* strongNode, impl::Node* weakNode);
    impl::Node* WalkToRoot(impl::Node* strongNode, impl::Node* weakNode, 
            impl::Arc* newArc);
    void Split(impl::Node* current, impl::Node* parent, impl::Arc* arc);
    void PushFlow(impl::Node* strongRoot);

private:
    bool m_NodePoolHasBeenInitialized;
    bool m_MinCutHasBeenSolved;
    std::unique_ptr<impl::NodePool> m_NodePool;
    std::unique_ptr<impl::PrecedenceArcPool> m_PrecedenceArcs;
    std::shared_ptr<const IPrecedenceConstraints> m_PrecedenceConstraints;
    std::vector<uint8_t> m_LargestSolution;

    ValueType m_PrevExcess;
};


////////////////////////////////////////////////////////////////////////////////
// UNDERSTAND ABOVE BEFORE PROCEEDING BELOW
//  What follows are implementations, helpers, and auxiliary functions / classes
////////////////////////////////////////////////////////////////////////////////

#ifndef MVD_USE_GMP
class VecBlockValues : public IBlockValues
{
public:
    VecBlockValues(IndexType numBlocks);
    VecBlockValues(std::vector<ValueType>&& values);
    VecBlockValues(std::initializer_list<int> values);
    VecBlockValues(IBlockValues* valuesToCopy);
    ~VecBlockValues();

    virtual IndexType NumBlocks() const override;
    virtual void BlockValue(IndexType blockIndex, ValueType* value) const override;

    // To behave like a vector
    typedef std::vector<ValueType>::const_iterator const_iterator;
    typedef std::vector<ValueType>::iterator iterator;

    const_iterator begin() const;
    const_iterator end() const;
    iterator begin();
    iterator end();

    ValueType operator[](IndexType blockIndex) const;
    ValueType& operator[](IndexType blockIndex);

    void SetBlockValueSI(IndexType blockIndex, int64_t v);

private:
    std::vector<ValueType> m_Values;
};
#else 
class GMPBlockValues : public IBlockValues {
public:
    GMPBlockValues(IndexType numBlocks);
    GMPBlockValues(IndexType numBlocks, int64_t initialValue);
    GMPBlockValues(const std::vector<int64_t>& initialValues);
    ~GMPBlockValues();

    virtual IndexType NumBlocks() const override;
    virtual void BlockValue(IndexType blockIndex, ValueType* value) const override;

    void SetBlockValueSI(IndexType blockIndex, int64_t v);

private:
    IndexType m_NumBlocks;
    std::unique_ptr<ValueType[]> m_BlockValues;
};
#endif 

// To solve for the largest minimum cut you can modify the values
// You pretty much MUST use gmp with this, unless your problem is teensy weensy
class SolveLargestValuesAdapter : public IBlockValues
{
public:
    SolveLargestValuesAdapter(std::shared_ptr<const IBlockValues> values);
    ~SolveLargestValuesAdapter();

    virtual IndexType NumBlocks() const override final;
    virtual void BlockValue(IndexType blockIndex, ValueType* value) const override final;

private:
    std::shared_ptr<const IBlockValues> m_Values;
    ValueType m_NumNonNegativeBlocks;
};


////////////////////////////////////////////////////////////////////////////////

// A simple regular block model definition. The number of blocks, it's origin,
// and block size/spacing

// Regular block models organized such that first block (1d index == 0) is the
// left-most (lowest x), front-most (lowest y), bottom-most (lowest z) block. 1D
// index increases by x fastest, then y, then z.
struct BlockDefinition 
{
    BlockDefinition();
    BlockDefinition(
            IndexType iNumX, IndexType iNumY, IndexType iNumZ,
            double iMinX, double iMinY, double iMinZ,
            double iSizeX, double iSizeY, double iSizeZ);
    ~BlockDefinition();

    IndexType NumBlocks() const;

    // Computes the 1d grid index from the 3d x,y,z indices
    // x The 3d grid index, x >= 0, x < NumX
    // y The 3d grid index, y >= 0, y < NumY
    // z The 3d grid index, z >= 0, z < NumZ
    IndexType GridIndex(IndexType x, IndexType y, IndexType z) const;

    // Computes the 3d x, y, or z index from the 1d grid index
    // idx The 1d grid index, idx >= 0, idx < NumBlocks()
    IndexType XIndex(IndexType idx) const;
    IndexType YIndex(IndexType idx) const;
    IndexType ZIndex(IndexType idx) const;
    std::tuple<IndexType, IndexType, IndexType> XYZIndices(IndexType idx) const;


    // Computes an offset 1d grid index
    // idx The 1d grid index, idx >= 0, idx < NumBlocks
    // ox, oy, oz The offsets
    IndexType OffsetIndex(IndexType idx, IndexType ox, IndexType oy, IndexType oz) const;

    // Returns if the block at the 3d indices would be inside this def
    // x, y, z The 3d grid indices as signed integers
    bool InDef(IndexType x, IndexType y, IndexType z) const;

    // Returns if the block at the 1d index would be inside this def
    // idx The 1d grid index as a signed integer
    bool InDef(IndexType idx) const;

    // Returns if the block at the offset would be inside this def
    // x, y, z The 3d grid indices
    bool OffsetInDef(IndexType x, IndexType y, IndexType z, IndexType ox, IndexType oy, IndexType oz) const;

    // Returns if the block at the offset would be inside this def
    // idx The 1d grid index
    // ox, oy, oz The 3d offsets
    bool OffsetInDef(IndexType idx, IndexType ox, IndexType oy, IndexType oz) const;

    // Static initializer for unit size block model
    // iNumX, iNumY, iNumZ The number of blocks in x, y, z
    static BlockDefinition UnitModel(IndexType iNumX = 1, IndexType iNumY = 1, IndexType iNumZ = 1);

    IndexType NumX; // Number of blocks in x direction
    IndexType NumY; // Number of blocks in y direction
    IndexType NumZ; // Number of blocks in z direction
    double MinX;    // Origin of blocks x
    double MinY;    // Origin of blocks y
    double MinZ;    // Origin of blocks z
    double SizeX;   // Size/spacing of blocks x
    double SizeY;   // Size/spacing of blocks x
    double SizeZ;   // Size/spacing of blocks x
};

constexpr inline long double ToDegrees(long double radians)
{
    constexpr long double TAU = 6.283185307179586476925286766559;
    return radians * 360.0 / TAU;
}
constexpr inline long double ToRadians(long double degrees)
{
    constexpr long double TAU = 6.283185307179586476925286766559;
    return degrees * TAU / 360.0;
}

constexpr long double operator"" _deg(long double degrees)
{
    return ToRadians(degrees);
}
constexpr long double operator"" _deg(unsigned long long degrees)
{
    return ToRadians(static_cast<long double>(degrees));
}

////////////////////////////////////////////////////////////////////////////////
// A vector class of marginal use

template <typename T, size_t n>
struct VectorBase
{
    VectorBase() {};
    ~VectorBase() {};

    T& operator[](size_t idx) { return data[idx]; };
    const T& operator[](size_t idx) const { return data[idx]; };

    VectorBase& operator+=(const VectorBase& rhs) {
        for (size_t i = 0; i < n; i++)
            data[i] += rhs[i];
        return *this;
    }
    VectorBase& operator-=(const VectorBase& rhs) {
        for (size_t i = 0; i < n; i++)
            data[i] -= rhs[i];
        return *this;
    }
    VectorBase& operator*=(const VectorBase& rhs) {
        for (size_t i = 0; i < n; i++)
            data[i] *= rhs[i];
        return *this;
    }
    VectorBase& operator/=(const VectorBase& rhs) {
        for (size_t i = 0; i < n; i++)
            data[i] /= rhs[i];
        return *this;
    }
    VectorBase& operator+=(T rhs) {
        for (size_t i = 0; i < n; i++)
            data[i] += rhs;
        return *this;
    }
    VectorBase& operator-=(T rhs) {
        for (size_t i = 0; i < n; i++)
            data[i] -= rhs;
        return *this;
    }
    VectorBase& operator*=(T rhs) {
        for (size_t i = 0; i < n; i++)
            data[i] *= rhs;
        return *this;
    }
    VectorBase& operator/=(T rhs) {
        for (size_t i = 0; i < n; i++)
            data[i] /= rhs;
        return *this;
    }
    VectorBase operator-() const {
        auto copy = *this;
        for (size_t i = 0; i < n; i++)
            copy[i] = -copy[i];
        return copy;
    }

    static VectorBase Origin() {
        VectorBase o;
        std::fill(o.data.begin(), o.data.end(), T(0));
        return o;
    }

    std::array<T, n> data;
};

// Binary operations with another vector
template <typename T, size_t n>
inline VectorBase<T, n> operator+(VectorBase<T, n> lhs, const VectorBase<T, n>& rhs)
{
    lhs += rhs; return lhs;
}
template <typename T, size_t n>
inline VectorBase<T, n> operator-(VectorBase<T, n> lhs, const VectorBase<T, n>& rhs)
{
    lhs -= rhs; return lhs;
}
template <typename T, size_t n>
inline VectorBase<T, n> operator*(VectorBase<T, n> lhs, const VectorBase<T, n>& rhs)
{
    lhs *= rhs; return lhs;
}
template <typename T, size_t n>
inline VectorBase<T, n> operator/(VectorBase<T, n> lhs, const VectorBase<T, n>& rhs)
{
    lhs /= rhs; return lhs;
}

// Binary operations with a scalar
template <typename T, size_t n, typename U>
inline VectorBase<T, n> operator+(VectorBase<T, n> lhs, U rhs)
{
    static_assert(std::is_convertible<T, U>::value, "Error: Source type not convertible to destination type.");
    lhs += rhs; return lhs;
}
template <typename T, size_t n, typename U>
inline VectorBase<T, n> operator-(VectorBase<T, n> lhs, U rhs)
{
    static_assert(std::is_convertible<T, U>::value, "Error: Source type not convertible to destination type.");
    lhs -= rhs; return lhs;
}
template <typename T, size_t n, typename U>
inline VectorBase<T, n> operator*(VectorBase<T, n> lhs, U rhs)
{
    static_assert(std::is_convertible<T, U>::value, "Error: Source type not convertible to destination type.");
    lhs *= rhs; return lhs;
}
template <typename T, size_t n, typename U>
inline VectorBase<T, n> operator/(VectorBase<T, n> lhs, U rhs)
{
    static_assert(std::is_convertible<T, U>::value, "Error: Source type not convertible to destination type.");
    lhs /= rhs; return lhs;
}
template <typename T, size_t n, typename U>
inline VectorBase<T, n> operator*(U lhs, VectorBase<T, n> rhs)
{
    static_assert(std::is_convertible<T, U>::value, "Error: Source type not convertible to destination type.");
    rhs *= lhs; return rhs;
}
template <typename T, size_t n, typename U>
inline VectorBase<T, n> operator/(U lhs, VectorBase<T, n> rhs)
{
    static_assert(std::is_convertible<T, U>::value, "Error: Source type not convertible to destination type.");
    rhs /= lhs; return rhs;
}

// Equality operators
template <typename T, size_t n>
inline bool operator==(const VectorBase<T, n>& lhs, const VectorBase<T, n>& rhs)
{
    return std::equal(lhs.data.begin(), lhs.data.end(), rhs.data.begin());
}
template <typename T, size_t n>
inline bool operator!=(const VectorBase<T, n>& lhs, const VectorBase<T, n>& rhs)
{
    return !operator==(lhs, rhs);
}

// Lexicographic comparisons
template <typename T, size_t n>
inline bool operator< (const VectorBase<T, n>& lhs, const VectorBase<T, n>& rhs)
{
    for (size_t i = 0; i < n - 1; i++)
        if (lhs[i] != rhs[i])
            return lhs[i] < rhs[i];
    return lhs[n - 1] < rhs[n - 1];
}
template <typename T, size_t n>
inline bool operator> (const VectorBase<T, n>& lhs, const VectorBase<T, n>& rhs)
{
    return operator<(rhs, lhs);
}
template <typename T, size_t n>
inline bool operator<=(const VectorBase<T, n>& lhs, const VectorBase<T, n>& rhs)
{
    return !operator>(lhs, rhs);
}
template <typename T, size_t n>
inline bool operator>=(const VectorBase<T, n>& lhs, const VectorBase<T, n>& rhs)
{
    return !operator<(lhs, rhs);
}

// General Operations
template <typename T, size_t n>
inline T Distance(const VectorBase<T, n>& a, const VectorBase<T, n>& b)
{
    return Magnitude(b - a);
}
template <typename T, size_t n>
inline T DistanceSquared(const VectorBase<T, n>& a, const VectorBase<T, n>& b)
{
    return MagnitudeSquared(b - a);
}
template <typename T, size_t n>
inline T Dot(const VectorBase<T, n>& lhs, const VectorBase<T, n>& rhs)
{
    auto f1 = lhs.data.begin(); auto l1 = lhs.data.end();
    auto f2 = rhs.data.begin();
    T v = T(0);
    while (f1 != l1) {
        v += *f1 * *f2;
        ++f1; ++f2;
    }
    return v;
}
template <typename T, size_t n>
inline T MagnitudeSquared(const VectorBase<T, n>& vec)
{
    return Dot(vec, vec);
}
template <typename T, size_t n>
inline T Magnitude(const VectorBase<T, n>& vec)
{
    return std::sqrt(MagnitudeSquared(vec));
}
template <typename T, size_t n>
inline T Theta(const VectorBase<T, n>& lhs, const VectorBase<T, n>& rhs)
{
    return std::acos(Dot(lhs, rhs) / (Magnitude(lhs) * Magnitude(rhs)));
}
template <typename T, size_t n>
inline VectorBase<T, n>& Normalize(VectorBase<T, n>& vec)
{
    vec /= Magnitude(vec);
    return vec;
}
template <typename T, size_t n>
inline VectorBase<T, n>* Normalize(VectorBase<T, n>* vec)
{
   (*vec) /= Magnitude(*vec);
   return vec;
}
template <typename T, size_t n>
inline VectorBase<T, n> Normalized(VectorBase<T, n> vec)
{
    return Normalize(vec);
}

// 2D operations, higher dimensions are ignored
template <typename T, size_t n>
inline double TriArea2(const VectorBase<T, n>& a, const VectorBase<T, n>& b, const VectorBase<T, n>& c)
{
    return ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]));
}
template <typename T, size_t n>
inline double TriArea(const VectorBase<T, n>& a, const VectorBase<T, n>& b, const VectorBase<T, n>& c)
{
    return TriArea2(a, b, c) / 2.0;
}
template <typename T, size_t n>
inline bool IsLeft(const VectorBase<T, n>& a, const VectorBase<T, n>& b, const VectorBase<T, n>& c)
{
    return TriArea2(a, b, c) > 0.0;
}
template <typename T, size_t n>
inline bool IsRight(const VectorBase<T, n>& a, const VectorBase<T, n>& b, const VectorBase<T, n>& c)
{
    return TriArea2(a, b, c) < 0.0;
}
template <typename T, size_t n>
inline bool IsCollinear(const VectorBase<T, n>& a, const VectorBase<T, n>& b, const VectorBase<T, n>& c)
{
    return TriArea2(a, b, c) == 0.0;
}

// Output
template <typename T, size_t n>
inline std::ostream& operator<<(std::ostream& os, const VectorBase<T, n>& vec)
{
    os << "{";
    for (size_t i = 0; i < n - 1; i++)
        os << vec[i] << ", ";
    os << vec[n - 1] << "}";
    return os;
}
template <typename T, size_t n>
inline std::istream& operator>>(std::istream& is, VectorBase<T, n>& vec)
{
    std::string str;
    std::getline(is, str, '{');
    for (size_t i = 0; i < n - 1; i++) {
        std::getline(is, str, ',');
        vec[i] = std::stod(str);
    }
    std::getline(is, str, '}');
    vec[n - 1] = std::stod(str);
    return is;
}

// Specializations for the most common dimensions
template <typename T>
struct VectorBase<T, 3>
{
    VectorBase() {};
    VectorBase(T x, T y, T z) : x(x), y(y), z(z) {};

    T& operator[](size_t idx) { return data[idx]; };
    const T& operator[](size_t idx) const { return data[idx]; };

    VectorBase& operator+=(const VectorBase& rhs) {
        x += rhs.x; y += rhs.y; z += rhs.z;
        return *this;
    }
    VectorBase& operator-=(const VectorBase& rhs) {
        x -= rhs.x; y -= rhs.y; z -= rhs.z;
        return *this;
    }
    VectorBase& operator*=(const VectorBase& rhs) {
        x *= rhs.x; y *= rhs.y; z *= rhs.z;
        return *this;
    }
    VectorBase& operator/=(const VectorBase& rhs) {
        x /= rhs.x; y /= rhs.y; z /= rhs.z;
        return *this;
    }
    VectorBase& operator+=(T rhs) {
        x += rhs; y += rhs; z += rhs;
        return *this;
    }
    VectorBase& operator-=(T rhs) {
        x -= rhs; y -= rhs; z -= rhs;
        return *this;
    }
    VectorBase& operator*=(T rhs) {
        x *= rhs; y *= rhs; z *= rhs;
        return *this;
    }
    VectorBase& operator/=(T rhs) {
        x /= rhs; y /= rhs; z /= rhs;
        return *this;
    }
    VectorBase operator-() const {
        return VectorBase(-x, -y, -z);
    }

    static VectorBase Origin() {
        return VectorBase(0.0, 0.0, 0.0);
    }
    static VectorBase XAxis() {
        return VectorBase(1.0, 0.0, 0.0);
    }
    static VectorBase YAxis() {
        return VectorBase(0.0, 1.0, 0.0);
    }
    static VectorBase ZAxis() {
        return VectorBase(0.0, 0.0, 1.0);
    }

    union {
        std::array<T, 3> data;
        struct { T x, y, z;};
    };
};
template <typename T>
VectorBase<T, 3> inline Cross(const VectorBase<T, 3>& lhs, const VectorBase<T, 3>& rhs)
{
    return VectorBase<T, 3>(lhs.y*rhs.z - lhs.z*rhs.y, lhs.z*rhs.x - lhs.x*rhs.z, lhs.x*rhs.y - lhs.y*rhs.x);
}

template <typename T>
struct VectorBase<T, 2>
{
    VectorBase() {};
    VectorBase(VectorBase<T, 3> vec) : x(vec.x), y(vec.y) {};
    VectorBase(T x, T y) : x(x), y(y) {};

    T& operator[](size_t idx) { return data[idx]; };
    const T& operator[](size_t idx) const { return data[idx]; };

    VectorBase& operator+=(const VectorBase& rhs) {
        x += rhs.x; y += rhs.y;
        return *this;
    }
    VectorBase& operator-=(const VectorBase& rhs) {
        x -= rhs.x; y -= rhs.y;
        return *this;
    }
    VectorBase& operator*=(const VectorBase& rhs) {
        x *= rhs.x; y *= rhs.y;
        return *this;
    }
    VectorBase& operator/=(const VectorBase& rhs) {
        x /= rhs.x; y /= rhs.y;
        return *this;
    }
    VectorBase& operator+=(T rhs) {
        x += rhs; y += rhs;
        return *this;
    }
    VectorBase& operator-=(T rhs) {
        x -= rhs; y -= rhs;
        return *this;
    }
    VectorBase& operator*=(T rhs) {
        x *= rhs; y *= rhs;
        return *this;
    }
    VectorBase& operator/=(T rhs) {
        x /= rhs; y /= rhs;
        return *this;
    }
    VectorBase operator-() const {
        return VectorBase(-x, -y);
    }

    static VectorBase Origin() {
        return VectorBase(0.0, 0.0);
    }
    static VectorBase XAxis() {
        return VectorBase(1.0, 0.0);
    }
    static VectorBase YAxis() {
        return VectorBase(0.0, 1.0);
    }

    union {
        std::array<T, 2> data;
        struct { T x, y; };
    };
};

typedef VectorBase<double, 3> Vector3D;
typedef VectorBase<double, 2> Vector2D;
typedef VectorBase<float, 3> Vector3F;
typedef VectorBase<float, 2> Vector2F;
typedef VectorBase<int64_t, 3> Vector3L;
typedef VectorBase<int64_t, 2> Vector2L;
typedef VectorBase<int, 3> Vector3I;
typedef VectorBase<int, 2> Vector2I;
typedef VectorBase<IndexType, 3> Vector3IT;

////////////////////////////////////////////////////////////////////////////////

// This InplaceLinspace will set up a container with a linspace between start
// and stop.
//
// std::vector<double> arr(100);
// InplaceLinspace(arr.begin(), arr.end(), 0, 14);
template <class Iter, typename T>
void InplaceLinspaceBase(Iter begin, Iter end, T start, T stop)
{
    typedef typename std::iterator_traits<Iter>::difference_type diff_t;
    typedef typename std::make_unsigned<diff_t>::type udiff_t;

    if (begin == end) return;
    udiff_t n = end - begin;

    T delta = stop - start;
    T step = delta / static_cast<T>(n - 1);

    udiff_t i = 0;
    for (auto it = begin; it != end; ++it) {
        *it = start + i * step;
        i++;
    }
}

template <class Iter>
void InplaceLinspace(Iter begin, Iter end, double start, double stop)
{
    InplaceLinspaceBase<Iter, double>(begin, end, start, stop);
}

// The linspace generator / iterator allows for the nice:
//
// for (auto v : linspace(0, 1, 100)) {
//     std::cout << v << std::endl;
// }
//
// without storing the entire thing in memory
template <typename T>
struct LinspaceIterator
{
    LinspaceIterator(T start, T step, int i) 
        : m_Start(start)
        , m_Step(step)
        , m_Index(i)
    {}
    ~LinspaceIterator() {}

    bool operator!= (const LinspaceIterator& other) const {
        return m_Index != other.m_Index;
    }

    T operator* () const {
        return m_Start + m_Index * m_Step;
    }

    const LinspaceIterator& operator++ () {
        ++m_Index;
        return *this;
    }

    T m_Start, m_Step;
    int m_Index;
};

// Generator for linspace
template <typename T>
struct LinspaceGeneratorBase
{
    LinspaceGeneratorBase(T start, T stop, int n)
        : m_Start(start)
        , m_N(n) {
        T delta = stop - start;
        m_Step = delta / static_cast<T>(n - 1);
    }
    ~LinspaceGeneratorBase() {};

    LinspaceIterator<T> begin() const {
        return LinspaceIterator<T>(m_Start, m_Step, 0);
    }
    LinspaceIterator<T> end() const {
        return LinspaceIterator<T>(m_Start, m_Step, m_N);
    }

    T m_Start, m_Step;
    int m_N;
};

// Convenience typedefs
typedef LinspaceGeneratorBase<double> Linspace;

////////////////////////////////////////////////////////////////////////////////

// Does not call constructor / destructor preferrable for plain old data types
template <int BlockSize, int WordSize>
class ObjectPoolBase
{
public:
   ObjectPoolBase()
      : m_Remaining(0)
      , m_CurrentBlock(nullptr)
      , m_CurrentLocation(nullptr) {
   }
   ~ObjectPoolBase() {
      while (m_CurrentBlock) {
         void* prev = *(static_cast<void**>(m_CurrentBlock));
         ::free(m_CurrentBlock);
         m_CurrentBlock = prev;
      }
   }

   void* InternalAlloc(size_t requiredSize) {
      size_t size = (requiredSize + WordSize - 1) & ~(WordSize - 1);

      void* location = m_CurrentLocation;
      if (size > m_Remaining) {
         void* newBlock = ::malloc(BlockSize);
         if (!newBlock) {
            return nullptr;
         }

         static_cast<void**>(newBlock)[0] = m_CurrentBlock;
         m_CurrentBlock = newBlock;

         m_Remaining = BlockSize - sizeof(void*);
         location = static_cast<char*>(newBlock) + sizeof(void*);
      }

      m_Remaining -= static_cast<int>(size);
      m_CurrentLocation = static_cast<char*>(location) + size;
      return location;
   }

   template<typename T>
   T* Alloc() {
      void* location = InternalAlloc(sizeof(T));
      if (location == nullptr) {
         return nullptr;
      }
      T* t = new (location) T();
      return t;
   }

private:
   int m_Remaining;
   void* m_CurrentBlock;
   void* m_CurrentLocation;
};

typedef ObjectPoolBase<8192, 64> ObjectPool;

////////////////////////////////////////////////////////////////////////////////

// An azimuth slope pair is a single component of a full slope definition
// Both the Azimuth and Slope are in radians
struct AzmSlopePair
{
    AzmSlopePair();
    AzmSlopePair(double iazm, double islope);
    ~AzmSlopePair();

    double Azimuth; // Azimuth in radians
    double Slope;   // Slope in radians

    bool operator<(const AzmSlopePair& other) const;
    bool operator<(double otherAzimuth) const;

};
std::ostream& operator<<(std::ostream& os, const AzmSlopePair& a);

// A slope definition is a sorted list of azimuth slope pairs.
//
// This will linearly interpolate for any requested azimuth, other interpolation
// techniques are supported by creating a very 'full' slope definition, say 512
// pairs and then linearly interpolating that.
//
// Again always radians
class SlopeDefinition
{
public:
    SlopeDefinition() {};
    SlopeDefinition(std::initializer_list<std::initializer_list<double>> list);
    SlopeDefinition(const std::vector<AzmSlopePair>& pairs);
    ~SlopeDefinition() {};

    // Compute the slope at the given azimuth
    // azimuth in RADIANS
    double Get(double azimuth) const;
    double operator()(double azimuth) const;

    // Computes if the given /vector/ is within the slope definition
    bool Within(double dx, double dy, double dz) const;
    bool Within(const Vector3D& vec) const;

    double MinSlope() const;

    size_t NumPairs() const;
    const std::vector<AzmSlopePair>& Pairs() const;
    bool Empty() const;

    static SlopeDefinition Constant(double slope);

private:
    std::vector<AzmSlopePair> m_Pairs; // The pairs; SORTED
};
std::ostream& operator<<(std::ostream& os, const SlopeDefinition& def);

// Cubic interpolation of the slope definition
SlopeDefinition CubicInterpolate(const SlopeDefinition& def, int cnt = 512);
// Cosine interpolation of the slope definition
SlopeDefinition CosineInterpolate(const SlopeDefinition& def, int cnt = 512);

// A pattern is a set of offsets from a base block
struct PrecedencePattern
{
    PrecedencePattern();
    ~PrecedencePattern();

    std::vector<Vector3IT> Offsets;
    using iterator = std::vector<Vector3IT>::iterator;
    using const_iterator = std::vector<Vector3IT>::const_iterator;

    size_t size() const;
    iterator begin();
    iterator end();
    const_iterator begin() const;
    const_iterator end() const;


    static PrecedencePattern OneFive();
    static PrecedencePattern OneNine();
    static PrecedencePattern KnightsMove();
    static PrecedencePattern Naive(const BlockDefinition& block_def,
                         const SlopeDefinition& slope_def,
                         IndexType numZ);
    static PrecedencePattern LessNaive(const BlockDefinition& block_def,
                             const SlopeDefinition& slope_def,
                             IndexType numZ);

    // The minimum search patterns are the best, they are the 'optimal' pattern
    // for a specific definition
    static PrecedencePattern MinSearch(const BlockDefinition& block_def,
                             const SlopeDefinition& slope_def,
                             IndexType numZ);
    static PrecedencePattern MinSearch(double slope_rad, IndexType numZ);
};
void PrintPattern(const PrecedencePattern& ptrn);

void NaiveSearch(const BlockDefinition& block_def,
                 const SlopeDefinition& slope_def,
                 IndexType numZ, 
                 std::function<void(Vector3IT)> offset_callback);

struct PatternAccuracy
{
    IndexType true_positive;
    IndexType true_negative;
    IndexType false_positive;
    IndexType false_negative;

    double accuracy;
    double true_positive_rate;
    double false_negative_rate;
    double matthews_correlation;
};
std::ostream& operator<<(std::ostream& os, const PatternAccuracy& acc);

void MeasureAccuracy(const BlockDefinition& block_def,
                     const SlopeDefinition& slope_def,
                     const PrecedencePattern& ptrn, PatternAccuracy* accuracy);
void MultiMeasureAccuracy(const BlockDefinition& block_def,
                          const SlopeDefinition& slope_def,
                          const PrecedencePattern& ptrn,
                          std::vector<PatternAccuracy>* accuracies);



class Regular2DGrid45DegreePrecedence : public IPrecedenceConstraints
{
public:
    Regular2DGrid45DegreePrecedence(IndexType numX, IndexType numZ);
    ~Regular2DGrid45DegreePrecedence();

    IndexType NumBlocks() const override final;
    BlockIndexInputIteratorBase Antecedents(IndexType fromBlockIndex) const override final;
    BlockIndexInputIteratorBase Successors(IndexType toBlockIndex) const override final;

private:
    BlockIndexInputIteratorBase XAdjustedSource(
            IndexType blockIndex, const std::vector<IndexType>& offsets) const;
    IndexType m_NumX;
    IndexType m_NumZ;

    std::vector<IndexType> m_AntecedentOffsets;
    std::vector<IndexType> m_SuccessorOffsets;
};

class Regular3DBlockModelPatternPrecedence : public IPrecedenceConstraints
{
public:
    Regular3DBlockModelPatternPrecedence(IndexType numX, IndexType numY, IndexType numZ,
            const PrecedencePattern& pattern);
    Regular3DBlockModelPatternPrecedence(const BlockDefinition& blockDef,
            const PrecedencePattern& pattern);
    ~Regular3DBlockModelPatternPrecedence();

    IndexType NumBlocks() const override final;
    BlockIndexInputIteratorBase Antecedents(IndexType fromBlockIndex) const override final;
    BlockIndexInputIteratorBase Successors(IndexType toBlockIndex) const override final;
    IndexType ApproxNumAntecedents(IndexType fromBlockIndex) const override final;

private:
    IndexType m_NumX, m_NumY, m_NumZ;

    std::vector<Vector3IT> m_Offsets;
    std::vector<IndexType> m_Precomputed1DOffsets;
    IndexType m_MaxOffsetZ;
    std::vector<IndexType> m_NumOffsetsByZMinus;
    struct {
        IndexType xLo, xHi;
        IndexType yLo, yHi;
    } m_InnerRegion;

    std::tuple<IndexType, IndexType, IndexType> XYZ(IndexType k) const;
};

//
class ExplicitPrecedence : public IPrecedenceConstraints
{
public:
    ExplicitPrecedence(IndexType numBlocks);
    ExplicitPrecedence(IndexType numBlocks,
            std::initializer_list<std::initializer_list<int>> l);
    ExplicitPrecedence(IndexType numBlocks, 
            std::unordered_map<IndexType, std::vector<IndexType>>&& antecedents);
    virtual ~ExplicitPrecedence();

    IndexType NumBlocks() const override final;
    BlockIndexInputIteratorBase Antecedents(IndexType fromBlockIndex) const override final;
    BlockIndexInputIteratorBase Successors(IndexType toBlockIndex) const override final;

    void AddPrecedenceConstraint(IndexType fromBlockIndex, IndexType toBlockIndex);

private:
    IndexType m_NumBlocks;
    std::unordered_map<IndexType, std::vector<IndexType>> m_Antecedents;
};

// Checks (primarily for testing) that precedence constraints are consistent:
// - Returns the correct counts
// - Successors and Antecedents are correctly related
// - All precedence constraints are correct
bool ConsistentPrecedenceConstraints(IPrecedenceConstraints* pre);


// 
class PrecedenceConstraintsReachableSearchBuffer
{
public:
    PrecedenceConstraintsReachableSearchBuffer(IndexType numBlocks);
    ~PrecedenceConstraintsReachableSearchBuffer();

    void NewSearch();
    void Queue(IndexType v);
    bool Search(IndexType* v);
    bool HasMore();

private:
    IndexType m_NumBlocks;
    uint8_t m_Tag;
    std::queue<IndexType> m_Queue;
    std::vector<uint8_t> m_Seen;
};

////////////////////////////////////////////////////////////////////////////////

// Input iterators over precedence constraints
class IPrecedenceConstraintInputIteratorSource
{
public:
    IPrecedenceConstraintInputIteratorSource();
    virtual ~IPrecedenceConstraintInputIteratorSource();

    virtual PrecedenceConstraint Next() = 0;
    virtual bool HasMore() const = 0;
};

class PrecedenceConstraintInputIteratorBase
{
public:
    class iterator {
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = PrecedenceConstraint;
        using reference = const PrecedenceConstraint&;
        using pointer = const PrecedenceConstraint*;
        using difference_type = std::ptrdiff_t;

        iterator(IPrecedenceConstraintInputIteratorSource* source = nullptr);

        reference operator*() const;
        iterator& operator++();
        iterator& operator++(int);
        bool operator==(iterator rhs) const;
        bool operator!=(iterator rhs) const;

    private:
        void Next();

        PrecedenceConstraint m_CurrentPrecedenceConstraint;
        IPrecedenceConstraintInputIteratorSource* m_Source;
    };

    PrecedenceConstraintInputIteratorBase(IPrecedenceConstraintInputIteratorSource* source);
    ~PrecedenceConstraintInputIteratorBase();

    iterator begin() const;
    iterator end() const;

private:
    IPrecedenceConstraintInputIteratorSource* m_Source;
};

// Input iterators over blocks
class IBlockIndexInputIteratorSource
{
public:
    IBlockIndexInputIteratorSource();
    virtual ~IBlockIndexInputIteratorSource();

    virtual IndexType Next() = 0;
    virtual bool HasMore() const = 0;
};

class BlockIndexInputIteratorBase
{
public:
    class iterator {
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = IndexType;
        using reference = const IndexType&;
        using pointer = const IndexType*;
        using difference_type = std::ptrdiff_t;

        iterator(IBlockIndexInputIteratorSource* source = nullptr);

        reference operator*() const;
        iterator& operator++();
        iterator& operator++(int);
        bool operator==(iterator rhs) const;
        bool operator!=(iterator rhs) const;

    private:
        void Next();

        IndexType m_CurrentIndex;
        IBlockIndexInputIteratorSource* m_Source;
    };

    BlockIndexInputIteratorBase(IBlockIndexInputIteratorSource* source);
    ~BlockIndexInputIteratorBase();

    iterator begin() const;
    iterator end() const;

private:
    IBlockIndexInputIteratorSource* m_Source;
};

class VecBlockSource : public mineflow::IBlockIndexInputIteratorSource
{
public:
    VecBlockSource(std::vector<IndexType>&& blocks);
    virtual ~VecBlockSource();

    virtual IndexType Next() override;
    virtual bool HasMore() const override;

private:
    size_t m_Index;
    std::vector<IndexType> m_Blocks;
};

}

#endif // MVD_INCLUDE_MINEFLOW_H

