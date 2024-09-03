# CoolCore - Poisson's Equation Solver

## Introduction
This project is a Poisson's equation solver designed to facilitate the analysis and visualization of steady-state temperature distributions in various systems, with a particular emphasis on thermal management within microprocessors.

## Getting Started

### Prerequisites
Before you begin, ensure you have the following requirements installed:
- Python 3.x
- Relevant Python libraries: NumPy, Matplotlib

### Installation
1. Unzip the zip file into your desired directory.
2. Verify that all extracted files are located within the same folder.

### Running the Solver
To execute the solver and generate the temperature distribution plots:
1. The source files are conveniently named to correspond with the figures detailed in my report.
2. To run a specific simulation, use the terminal or command line interface to execute the Python script associated with that plot:
  ```
  python fig_5.py
  ```
4. Note that in some of the scripts I have included precomputed values, as certain simulations can be computationally intensive.

### Saving plots
To locally save temperature distribution plots:
1. Locate the `create_plot` function within the 'solver_all.py', 'colver_CM.py', or 'validator.py' script.
2. Modify the function to include your desired file path where the plot image will be saved:
```python
plt.savefig('/your/file/path/here', dpi=500)
```

