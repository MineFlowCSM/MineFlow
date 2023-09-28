/*
* mineflowmex.cpp
*
* Compute a very simple ultimate pit from a 3d matrix 
*/

#include "mex.hpp"
#include "mexAdapter.hpp"

#include "../mineflow.cpp"


class MexFunction : public matlab::mex::Function
{
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs);

    void error(const std::string& message);
    void validateArgumentsOrError(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs);
};

void MexFunction::operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
{
    validateArgumentsOrError(outputs, inputs);
    try {
        double slopeRad = mvd::mineflow::ToRadians(inputs[1][0]);
        mvd::mineflow::SlopeDefinition sdef = mvd::mineflow::SlopeDefinition::Constant(slopeRad);
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

        // 12 is arbitrary, but probably fine
        int numz = std::min(static_cast<int>(bdef.NumZ), 12);
        mvd::mineflow::PrecedencePattern pattern = mvd::mineflow::PrecedencePattern::MinSearch(bdef, sdef, numz);

        auto pre = std::make_shared<mvd::mineflow::Regular3DBlockModelPatternPrecedence>(bdef, pattern);

        // TODO could probably avoid copying the data..
        // But MATLAB is column major order, and well. who knows
        mvd::mineflow::VecBlockValues values(inputs[0].getNumberOfElements());

        size_t mexIndex = 0;
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

        // TODO could print info somehow

        matlab::data::ArrayFactory factory;
        auto soln = factory.createArray<bool>(inputs[0].getDimensions());
        mexIndex = 0;
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
        error("Two inputs are required. The economic block values, and the slope (in degrees)");
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

