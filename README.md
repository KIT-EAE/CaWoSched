# CaWoSched

**CaWoSched: Carbon-aware Workflow Scheduling with Fixed Mapping and Deadline Constraint**

A scheduling framework for optimizing workflow execution on heterogeneous multi-processor systems while minimizing carbon emissions and meeting deadline constraints.

If you use CaWoSched in your research, please cite:

```bibtex
@misc{schweisgut2025carbonawareworkflowschedulingfixed,
      title={Carbon-Aware Workflow Scheduling with Fixed Mapping and Deadline Constraint}, 
      author={Dominik Schweisgut and Anne Benoit and Yves Robert and Henning Meyerhenke},
      year={2025},
      eprint={2507.08725},
      archivePrefix={arXiv},
      primaryClass={cs.DC},
      url={https://arxiv.org/abs/2507.08725}, 
}
```

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Input File Formats](#input-file-formats)
- [Scheduling Algorithms](#scheduling-algorithms)
- [License](#license)

## Overview

CaWoSched is a carbon-aware workflow scheduling framework designed to minimize carbon emissions during workflow execution while ensuring that the workflow execution meets a specified deadline. The system supports multiple heuristic scheduling algorithms and an ILP formulation of the problem to solve small instances optimally.


## Requirements

### System Requirements

- **Compiler**: C++17 compatible compiler (GCC 7+, Clang 5+, ...)
- **Build System**: CMake 3.10 or higher

### Dependencies

#### Core Dependencies (C++)
- CMake 3.10+
- C++17 compatible compiler
- Standard library support for filesystem operations

#### Optional Dependencies (Python ILP solver)
- Python 3.7+
- Gurobi Optimizer (requires license)
- NetworkX library
- Matplotlib (for visualization of experiment results)

### Installing Python Dependencies

```bash
pip install gurobipy networkx matplotlib
```

**Note**: Gurobi requires a valid license. 

## Installation

### Quick Installation

1. **Clone the repository**:
   ```bash
   git clone <repository-url>
   cd CaWoSched
   ```

2. **Build the project**:
   ```bash
   chmod +x build.sh
   ./build.sh
   ```

3. **You should have the following executable now**:
   ```bash
   ./deploy/multi_machine_scheduler
   ```

### Manual Build

If you prefer manual building:

```bash
cd CaWoSched
mkdir build && cd build
cmake ..
cmake --build . --parallel $(nproc)
```

## Usage

### Basic Usage (Heuristics)

```bash
./deploy/multi_machine_scheduler <DAG_file> <Mapping_file> <Setup_file>
```

### Basic Usage (ILP Solver)

```bash
cd ILP
python multiple_machines_no_freezing_hetero.py \
    --graph workflow.dag \
    --mapping mapping.txt \
    --setup setup.txt
```

### Parameters

- **DAG_file**: Workflow description file (Directed Acyclic Graph) in a .dot format
- **Mapping_file**: Stores the given mapping
- **Setup_file**: System setup and green power profile

An example is provided in the ``examples/`` folder.

### Output

The scheduler outputs detailed results for all implemented algorithms:
- **Variant**: Algorithm variant used (slackR, pressureR, etc.)
- **Workflow**: Name of the processed workflow
- **PEs**: Number of processing elements
- **Deadline**: Specified deadline constraint
- **Execution time**: Algorithm runtime
- **Validity**: Whether the solution meets constraints
- **Cost**: Carbon cost

## Input File Formats

### DAG File Format
The DAG file must be a .dot file with a structure following this example:
```
strict digraph "G" {
0 [weight=5];
1 [weight=4];
0 -> 1 [size=3];
}
```

Hence, to specify a vertex use 
```
vertex_id [weight=vertex_weight];
```
and to specify an edge use 
```
source -> target [size=communication_volume];
```

### Mapping File Format
The given mapping should be provided as a .txt file with the following structure:
```
numProcessors
tasks_of_processor_0
tasks_of_processor_1
.
.
.
tasks_of_processor_k-1
```
where the first line specifies the number of processors, and line x specifies the tasks assigned to processor x-1 (given in a order that can be respected and separated by a space).

### Setup File Format
The setup file contains information about the cluster configuration and the green power profile. It should be given as a .txt file. The general structure is as follows:
```
deadline 
numIntervals
idle_power_0 idle_power_1 ... idle_power_P^2-1
work_power_0 work_power_1 ... work_power_P^2-1
brown_carbon_cost
green_power_budget_0 green_power_budget_1 ... green_power_budget_numIntervals-1
b_0=0 e_0=b_1 e_1=b_2 ... e_numIntervals-1=deadline 
```
The first line specifies the given deadline as an integer. The second line specifies the number of given intervals during which the green power budget is constant. The second line contains the idle power consumption of the processors (including communication processors) separated by a space. The third line contains the corresponding work power values.
The fourth line contains the brown carbon cost we have to pay for every energy unit consumed more than the green power budget. The fifth line contains the green power budget for each of the intervals. The last line specifies the start and endpoints of each interval, i.e., $I_j = [b_j, e_j[$.

### Example

Examples for workflow files, mappings and setups can be found in the ``experiments/`` folder.

## Scheduling Algorithms

CaWoSched implements multiple variants of the heuristic framework:

### Primary Algorithms

1. **`slack(-LS)`**: Base score is `slack`. `-LS` if local search is applied.
2. **`press(-LS)`**: Base score is `pressure`. `-LS` if local search is applied. 
3. **`slackR(-LS)`**: Base score is `slack`. Interval refinement is used. `-LS` if local search is applied.
4. **`pressR(-LS)`**: Base score is `pressure`. Interval refinement is used. `-LS` if local search is applied..
5. **`slackW(-LS)`**: Base score is `slack` together with a weight factor to account for the power heterogeneity of the cluster. `-LS` if local search is applied.
6. **`pressW(-LS)`**: Base score is `pressure` together with a weight factor to account for the power heterogeneity of the cluster. `-LS` if local search is applied.
7. **`slackWR(-LS)`**: Base score is `slack` together with a weight factor to account for the power heterogeneity of the cluster and interval refinement is used. `-LS` if local search is applied.
8. **`pressWR(-LS)`**: Base score is `slack` together with a weight factor to account for the power heterogeneity of the cluster and interval refinement is used. `-LS` if local search is applied.
---

## License

License: [MIT](LISENCE)





