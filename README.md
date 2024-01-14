MineFlow 
========

An open source library and command line utility for calculating **ultimate pit limits** for open pit mining operations.


Introduction
------------

MineFlow is developed by Matthew Deutsch, with theoretical assistance from Dr. Kadri Daǧdelen and Dr. Thys Johnson. 
MineFlow uses a custom implementation of Dr. Dorit Hochbaum's Pseudoflow algorithm which has been modified herein specifically for mining applications.
With this engine MineFlow is able to solve for ultimate pits quickly and efficiently.
MineFlow is well suited for calculating many different ultimate pits in the presence of uncertainty, or as a part of sensitivity studies.

MineFlow contains techniques for building efficient precedence graphs to model pit slopes that vary by both location and direction.
MineFlow is not limited to regular 3D block models, it supports sub blocks and other irregular models and can be extended for custom applications.


Building Instructions
---------------------

On Linux without cmake:
```
g++ -std=c++1z -O3 -DMVD_MINEFLOW_EXE mineflow.cpp -o mineflow
./mineflow
```

On Linux with cmake:
```
mkdir build && cd build
cmake ..
make
```

On Windows with cmake.
```
# Something like this:
"C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvarsall.bat" amd64

mkdir build
cd build
cmake -G "Visual Studio 16 2019" -A x64 ..
msbuild ALL_BUILD.vcxproj /p:Configuration=Release
```

However, the primary usecase is to include `mineflow.h` and `mineflow.cpp` into your own build system.
Contributions to ease this process are welcome.

Usage
-----

As an executable user follow the above instructions to build the executable for your platform.
Then run it with no arguments and read the printed usage information.
Finally review the examples in the examples directory.
Note: a useful way to open a command prompt on windows from windows explorer is to simply type 'cmd' in the address bar.


To implement in your own code first decide if you want to use GMP or not.
In my non-exhaustive comparison GMP increases runtime by about 80 percent, but you will never have to worry about overflowing which might be appropriate for your application.

In the below examples I will not use GMP just to keep things easy.

```
#include "mineflow.h"
using namespace mvd::mineflow;

int main(int argc, char** argv)
{
    int64_t numBlocks = 5;
    auto values = std::make_shared<VecBlockValues>(numBlocks);
    values->SetBlockValueSI(0,  7); // 0 being the index, 7 the economic block value
    values->SetBlockValueSI(1,  2);
    values->SetBlockValueSI(2, -2);
    values->SetBlockValueSI(3, -2);
    values->SetBlockValueSI(4, -2);

    auto precedence = std::make_shared<ExplicitPrecedence>(numBlocks);
    precedence->AddPrecedenceConstraint(0, 2);
    precedence->AddPrecedenceConstraint(0, 3);
    precedence->AddPrecedenceConstraint(1, 3);
    precedence->AddPrecedenceConstraint(1, 4);

    PseudoSolver solver(precedence, values);
    solver.Solve();

    assert( solver.InMinimumCut(0));
    assert(!solver.InMinimumCut(1));
    assert( solver.InMinimumCut(2));
    assert( solver.InMinimumCut(3));
    assert(!solver.InMinimumCut(4));
    
    return 0;
}
```

```
// Solving a large 3D regular block model using the c++ library
#include "mineflow.h"
using namespace mvd::mineflow;

int main(int argc, char** argv)
{
    BlockDefinition bdef = BlockDefinition::UnitModel(120, 120, 26); // numx, numy, numz
    auto values = ReadTestDataValues(...); // you'll have to define this!
    auto pattern = PrecedencePattern::MinSearch(45_deg, 8); // A minimum search pattern of 45 degrees with 8 benches
    auto precedence = std::make_shared<Regular3DBlockModelPatternPrecedence>(bdef, pattern);

    PseudoSolver solver(precedence, values);

    PseudoSolverSolveInfo info;
    solver.Solve(&info);

    std::cout << info.ContainedValue << " [" << info.NumContainedNodes << "/" << info.NumNodes << "]" << std::endl;
    std::cout << info.ElapsedSeconds << "s" << std::endl;

    // use solver.InMinimumCut(blockIndex) for something!

    return 0;
}
```

Further information is in `mineflow.h`, and there are examples in `mineflow.cpp` in the test section.
Contributions to improve documentation / comments are welcome.


License
-------

MineFlow is licensed under the open source and permissive MIT license.
It can be used with extremely minimal restrictions even for commercial software.
See the LICENSE file for details.

We chose this license so that in our own way we can advance the mining industry as a whole.
MineFlow is available for *everyone* to use, learn from, and improve, without restriction or expense.

If you use this code commercially carefully note the requirements in the license.


Citation and Paper
------------------

If you use this code for your research, please cite our paper:

Deutsch, M., Dağdelen, K. & Johnson, T. "An Open-Source Program for Efficiently Computing Ultimate Pit Limits: MineFlow." Natural Resources Research (2022). https://doi.org/10.1007/s11053-022-10035-w

```
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
```

The paper is included in this repository: [deutsch2022mineflow.pdf](https://github.com/MineFlowCSM/MineFlow/blob/df0f30aabea494371704a926ba47f6166631774d/deutsch2022mineflow.pdf).
However note that this is the Author's "Submitted Manuscript."
This preprint has not undergone any post-submission improvements or corrections.
The **Version of Record** of this article is published in **Natural Resources Research** and is available online at [https://doi.org/10.1007/s11053-022-10035-w](https://doi.org/10.1007/s11053-022-10035-w)


Contact
-------

If MineFlow is useful to you we would love to hear about it.
Contact Matthew Deutsch by email:

* matthewvdeutsch@gmail.com

Acknowledgements
----------------

Hochbaum's pseudoflow algorithm forms the mathematical basis for MineFlow's underlying engine.
More information on Dr. Dorit Hochbaum and the underlying pseudoflow algorithm is [available here](https://riot.ieor.berkeley.edu/Applications/Pseudoflow/maxflow.html).
The minimum search patterns of Caccetta and Giannini 1988 form the basis for the minimum search pattern routines.

