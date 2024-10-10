# Dijkstra Shape Match Model with Coupling

This repository creates the integer linear program introduced in the CVPR2024 best paper award candidate SpiderMatch. For more information take a look at [our project page](https://paulroetzer.github.io/publications/2024-06-19-spidermatch.html) as well as the [main repository](https://github.com/paul0noah/spider-match) to which this repository serves as a helper.

## ⚙️ Installation

This repo comes with c++ code and python bindings 💡

Note: Code only tested on unix systems. Expect issues with Windows 🪟.

🐍 Python

Simply download the repo and `python setup.py install` (you need working c++ compiler + cmake installed on your machine).

🚀 Cpp

Simply include the singlemost header `shapeMatchModel/shapeMatchModel.hpp` file of this project into your project.

## ✨ Usage
🐍 Python
```python
from smm_dijkstra import *
from gurobipy import GRB


smm = ShapeMatchModelDijkstra(vy, ey, vx, ex, absfeatdiffmat, True, True)
E = smm.getCostVector()
RHS = smm.getRHS()
I, J, V = smm.getAVectors()

# Create a new model
m = gp.Model("matrix1")

# Create variables
x = m.addMVar(shape=E.shape[0], vtype=GRB.CONTINUOUS, name="x")

# Set objective
m.setObjective(E.transpose() @ x, GRB.MINIMIZE)

# Build (sparse) constraint matrix
A = sp.csr_matrix((V.flatten(), (I.flatten(), J.flatten())), shape=(RHS.shape[0], E.shape[0]))

# Add constraints
m.addConstr(A @ x == RHS.flatten(), name="c")
m.addConstr(x >= 0, name="x0")
m.addConstr(x <= 1, name="x1")
#m.addConstr(x[0] == 1, name='init')

# Optimize model
m.optimize()
```

🚀 Cpp
Please create a github issue if you want to have instructions on how to use the vanilla cpp implemention.

## Attribution 🎓
When using this code for your own projects please cite the followig:

```bibtex
@inproceedings{roetzer2024spidermatch,
    author     = {Paul Roetzer and Florian Bernard},
    title     = { SpiderMatch: 3D Shape Matching with Global Optimality and Geometric Consistency },
    booktitle = {IEEE Conference on Computer Vision and Pattern Recognition (CVPR)},
    year     = 2024
}
```

## License 🚀
This repo is licensed under MIT licence.
