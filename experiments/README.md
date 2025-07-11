# CaWoSched Experiments

This directory contains comprehensive experimental configurations and data for evaluating the CaWoSched framework.

## Dependencies

### Required Software
- **CaWoSched**: Built and deployed 
- **Python 3.7+**: For the ILP solver
- **Gurobi Optimizer**: For optimal ILP solutions 
- **simexpal**: For experiment management

## Running Experiments

### Heuristic Algorithm Experiments

The `simex_CaWoSched` directory contains configurations for testing all heuristic variants. We use [simexpal](https://github.com/hu-macsy/simexpal) to manage our experiments. Use 

```bash
cd simex_CaWoSched
simex e list
```
to list all experiments. Then you can launch the experiments directly via simexpal (see this [link](https://simexpal.readthedocs.io/en/latest/launcher.html) for more information about how to launch with simexpal).

### ILP Solver Experiments

For small instances that can be solved optimally:

```bash
cd simex_ILP
simex e list
```
Again, then you can launch the experiments directly via simexpal (see this [link](https://simexpal.readthedocs.io/en/latest/launcher.html) for more information about how to launch with simexpal).

