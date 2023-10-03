# Shape Match Model Dijkstra 3D Coupling

Description: tbd

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


coupling = True
lineIntegral = True
smm = ShapeMatchModelDijkstra(vy, ey, vx, ex, absfeatdiffmat, coupling, lineIntegral)
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
```c++
tbd
```

## Attribution 🎓
When using this code for your own projects please cite the followig:

```bibtex
tbd
```

## License 🚀
This repo is licensed under MIT licence.
