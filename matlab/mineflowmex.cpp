/*
* mineflowmex.cpp
*
* Compute a very simple ultimate pit from a 3d matrix 
*/

#include "mex.hpp"
#include "mexAdapter.hpp"

#include "mineflow.cpp"

typedef std::shared_ptr<const mvd::mineflow::IPrecedenceConstraints> PrecedencePtr;

class MexFunction : public matlab::mex::Function
{
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs);

    static void error(const std::string& message);

private:
    void validateArgumentsOrError(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs);
    PrecedencePtr initializePrecedenceConstraints(const mvd::mineflow::BlockDefinition& bdef, matlab::data::Array input);

    PrecedencePtr initializeScalarPrecedenceConstraints(const mvd::mineflow::BlockDefinition& bdef, double slopeRad);

    void validateSlopeInputOrThrow(const mvd::mineflow::BlockDefinition& bdef, matlab::data::Array input);
    PrecedencePtr initializeLocallyVaryingPrecedenceConstraints(const mvd::mineflow::BlockDefinition& bdef, matlab::data::Array input);

};


PrecedencePtr initializeScalarPrecedenceConstraints(const mvd::mineflow::BlockDefinition& bdef, double slopeRad)
{
    mvd::mineflow::SlopeDefinition sdef = mvd::mineflow::SlopeDefinition::Constant(slopeRad);

    // 12 is arbitrary, but probably fine
    int numz = std::min(static_cast<int>(bdef.NumZ), 12);
    mvd::mineflow::PrecedencePattern pattern = mvd::mineflow::PrecedencePattern::MinSearch(bdef, sdef, numz);

    return std::make_shared<mvd::mineflow::Regular3DBlockModelPatternPrecedence>(bdef, pattern);
}


void MexFunction::validateSlopeInputOrThrow(const mvd::mineflow::BlockDefinition& bdef, matlab::data::Array input)
{
    auto& dim = input.getDimensions();
    int ndim = dim.size();
    bool twod = ndim == 2;
    if (ndim == 2) {
        if (bdef.NumX != dim[1]) {
            throw std::runtime_error("dimension mismatch between EBV and SLOPE");
        }
        if (bdef.NumZ != dim[0]) {
            throw std::runtime_error("dimension mismatch between EBV and SLOPE");
        }
    } else if (ndim == 3) {
        if (bdef.NumX != dim[2]) {
            throw std::runtime_error("dimension mismatch between EBV and SLOPE");
        }
        if (bdef.NumY != dim[1]) {
            throw std::runtime_error("dimension mismatch between EBV and SLOPE");
        }
        if (bdef.NumZ != dim[0]) {
            throw std::runtime_error("dimension mismatch between EBV and SLOPE");
        }
    } else {
        throw std::runtime_error("invalid dimensions for slope input");
    }
}

PrecedencePtr MexFunction::initializeLocallyVaryingPrecedenceConstraints(const mvd::mineflow::BlockDefinition& bdef, matlab::data::Array input)
{
    validateSlopeInputOrThrow(bdef, input);

    int numz = std::min(static_cast<int>(bdef.NumZ), 12);

    std::unordered_map<int, mvd::mineflow::IndexType> inDegreesToPatternIndex;

    std::vector<mvd::mineflow::PrecedencePattern> patterns;
    std::shared_ptr<std::vector<mvd::mineflow::IndexType>> patternIndices(bdef.NumBlocks());

    bool twod = input.getDimensions().size() == 2;
    for (int64_t z = 0; z < bdef.NumZ; z++) {
        for (int64_t y = 0; y < bdef.NumY; y++) {
            for (int64_t x = 0; x < bdef.NumX; x++) {

                double slopeDeg = 0.0;
                if (twod) {
                    slopeDeg = input[z][x];
                } else {
                    slopeDeg = input[z][y][x];
                }

                if (slopeDeg < 5) {
                    throw std::runtime_error("All slopes must be given in degrees and be greater than 5 degrees");
                }
                if (slopeDeg > 90) {
                    throw std::runtime_error("All slopes must be given in degrees and be less than 90 degrees");
                }

                int inDegrees = static_cast<int>(std::round(slopeDeg * 10));
                mvd::IndexType patternIndex = -1;

                auto it = inDegreesToPatternIndex.find(inDegrees);
                if (it == inDegreesToPatternIndex.end()) {
                    // create a new pattern
                    mvd::mineflow::SlopeDefinition sdef = mvd::mineflow::SlopeDefinition::Constant(slopeRad);

                    patternIndex = static_cast<mvd::mineflow::IndexType>(patterns.size());
                    patterns.push_back(mvd::mineflow::PrecedencePattern::MinSearch(bdef, sdef, numz));

                    inDegreesToPatternIndex[inDegrees] = patternIndex;
                } else {
                    patternIndex = *it;
                }
                assert(patternIndex >= 0 && patternIndex < patterns.size());

                int64_t i = bdef.GridIndex(x, y, z);
                patternIndices[i] = patternIndex;
            }
        }
    }

    return std::make_shared<mvd::mineflow::Regular3DBlockModelKeyedPatternsPrecedence>(bdef, patterns, patternIndices);
}

PrecedencePtr MexFunction::initializePrecedenceConstraints(const mvd::mineflow::BlockDefinition& bdef, matlab::data::Array input)
{
    if (input.getNumberOfElements() == 1) {
        return initializeScalarPrecedenceConstraints(bdef, mvd::mineflow::ToRadians(input[0]));
    } else {
        error("not implemented yet");
    }
    return nullptr;
}


void MexFunction::operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
{
    validateArgumentsOrError(outputs, inputs);
    try {
        mvd::mineflow::BlockDefinition bdef = mvd::mineflow::BlockDefinition::UnitModel();
        int ndim = inputs[0].getDimensions().size();
        bool twod = ndim == 2;
        if (ndim == 2) {
            bdef.NumX = inputs[0].getDimensions()[1];
            bdef.NumY = 1;
            bdef.NumZ = inputs[0].getDimensions()[0];
        } else if (ndim == 3) {
            bdef.NumX = inputs[0].getDimensions()[2];
            bdef.NumY = inputs[0].getDimensions()[1];
            bdef.NumZ = inputs[0].getDimensions()[0];
        } else {
            throw std::runtime_error("invalid dimension");
        }

        auto pre = initializePrecedenceConstraints(bdef, inputs[1]);
        if (!pre) {
            throw std::runtime_error("error initializing precedence constraints");
        }

        // TODO could probably avoid copying the data..
        // But MATLAB is column major order, and well. who knows
        mvd::mineflow::VecBlockValues values(inputs[0].getNumberOfElements());

        auto& ebv = inputs[0];
        for (int64_t z = 0; z < bdef.NumZ; z++) {
            for (int64_t y = 0; y < bdef.NumY; y++) {
                for (int64_t x = 0; x < bdef.NumX; x++) {
                    int64_t i = bdef.GridIndex(x, y, z);
                    if (twod) {
                        values[i] = ebv[z][x];
                    } else {
                        values[i] = ebv[z][y][x];
                    }
                }
            }
        }

        mvd::mineflow::PseudoSolver solver(pre, &values);
        mvd::mineflow::PseudoSolverSolveInfo info;
        solver.Solve(&info);

        // TODO could print info somehow?

        matlab::data::ArrayFactory factory;
        auto soln = factory.createArray<bool>(inputs[0].getDimensions());
        for (int64_t x = 0; x < bdef.NumX; x++) {
            for (int64_t y = 0; y < bdef.NumY; y++) {
                for (int64_t z = 0; z < bdef.NumZ; z++) {
                    int64_t i = bdef.GridIndex(x, y, z);

                    if (twod) {
                        soln[z][x] = solver.InMinimumCut(i);
                    } else {
                        soln[z][y][x] = solver.InMinimumCut(i);
                    }
                }
            }
        }

        outputs[0] = std::move(soln);
    } catch (std::exception& e) {
        error(e.what());
    }

}

void MexFunction::error(const std::string& message)
{
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
    matlab::data::ArrayFactory factory;

    matlabPtr->feval(u"error", 
        0, std::vector<matlab::data::Array>({ 
            factory.createScalar(message.c_str()) }));
}

void MexFunction::validateArgumentsOrError(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
{
    if (inputs.size() != 2) {
        error("Two inputs are required. The economic block values, and the slope information (in degrees)");
    }
    if (inputs[1].getNumberOfElements() != 1) {
        error("The slope argument must be a scalar.");
    }
    if (inputs[1].getType() != matlab::data::ArrayType::DOUBLE ||
        inputs[1].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE) {
        error("The slope argument must be a noncomplex scalar double.");
    }

    int ndim = inputs[0].getDimensions().size();
    if (ndim < 2 || ndim > 3) {
        error("The economic block values must be a 2d or 3d matrix.");
    }
    if (inputs[0].getType() != matlab::data::ArrayType::INT64) {
        error("The economic block values must be provided as an int64 matrix.");
    }

    double slope = inputs[1][0];
    if (slope < 5) {
        error("The input slope must be given in degrees and be greater than 5 degrees");
    }
    if (slope > 90) {
        error("The input slope must be given in degrees and be less than 90 degrees");
    }
}

