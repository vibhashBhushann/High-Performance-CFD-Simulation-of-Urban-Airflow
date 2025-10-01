# High-Performance CFD Simulation of Urban Airflow

This project is a C++ based application developed to model and analyze 3D fluid dynamics around complex architectural geometries using the Lattice Boltzmann Method (LBM). The primary case study investigates the aerodynamic impact of the "Lal Minar" structure on wind flow within the IIT Gandhinagar campus.

The simulation engine is designed for high performance, leveraging parallel computing techniques to handle computationally intensive tasks.

## 🚀 Technical Features

- **C++ Simulation Engine**: The core logic is built in C++ for performance and low-level control over the simulation process.

- **High-Performance Parallel Computing**: Utilizes the Palabos library and MPI (Message Passing Interface) to distribute computations across multiple processor cores, enabling efficient, large-scale 3D simulations.

- **Complex 3D Geometry**: Integrates complex 3D models in STL format (created using Gmsh), allowing for the accurate representation of real-world obstacles.

- **Configurable & Extensible**: Simulation parameters (like viscosity, inlet velocity, and grid resolution) are managed through an external XML file, decoupling configuration from the source code and allowing for easy experimentation.

- **Data Processing & Visualization**: The engine generates output in the VTK file format, creating a data pipeline for analysis and visualization with tools like ParaView.

## 🛠️ Tech Stack

- **Core Language**: C++
- **CFD Library**: Palabos
- **Parallel Computing**: MPI
- **3D Modeling**: Gmsh (for STL generation)
- **Configuration**: XML
- **Data Output**: VTK (Visualization Toolkit)

## 🧠 Engineering Challenges & Solutions

A key part of this project involved overcoming significant performance challenges to enable large-scale 3D simulations.

**Challenge**: The initial implementation using the OpenLB library was effective for 2D models but faced severe performance bottlenecks in 3D. The library's support for MPI was insufficient for the required scale, making complex simulations computationally infeasible.

**Solution**: We strategically migrated the entire codebase from OpenLB to the Palabos library. This decision was driven by Palabos's superior architecture for distributed memory systems and robust MPI support. This migration was critical to the project's success, resulting in a scalable and high-performance application capable of handling complex 3D environments.

## ⚙️ How to Run the Simulation

### Prerequisites

- A C++ compiler (like `g++`)
- An MPI implementation (like Open MPI)
- The Palabos library installed and configured on your system

### 1. Compile the Application

Navigate to the project directory and compile the C++ source file. You will need to link against your Palabos installation.

```bash
g++ -o externalFlow externalFlowAroundObstacle.cpp -I/path/to/palabos/src -L/path/to/palabos/lib -lpalabos
```

Replace `/path/to/palabos/` with the actual path to your Palabos library installation.

### 2. Run the Simulation

Execute the compiled program using `mpirun` to enable parallel processing. The number of processes (`-np`) can be adjusted based on your system's capabilities. The program takes the XML configuration file as an argument.

```bash
mpirun -np 4 ./externalFlow externalFlowAroundObstacle.xml
```

The simulation will start, and output files (including `.vtk` files for visualization) will be generated in the specified output directory.

## 📁 Project Files

- `externalFlowAroundObstacle.cpp`: The main C++ source code for the simulation.
- `externalFlowAroundObstacle.xml`: The XML configuration file that defines all simulation parameters.
- `Unnamed3-final4.stl`: The 3D geometry file representing the physical environment.

## 📊 Visualization

### Simulation Results

<img width="800" alt="Simulation visualization 1" src="https://github.com/user-attachments/assets/f7f77b0c-788e-45a1-ae41-608105df4f7c" />

<img width="800" alt="Simulation visualization 2" src="https://github.com/user-attachments/assets/569ca25f-e2f8-4d02-b3e0-f7bb6bd68f95" />

<img width="800" alt="Simulation visualization 3" src="https://github.com/user-attachments/assets/52844fd3-3ee7-4b88-bdc0-904088dd6247" />

## 👥 Authors & Acknowledgements

**Authors**: Keshav Sobania, Pratyaksh Bhayre, Vibhash Bhushan

**Project Guide**: Prof. Sushobhan Sen

This project was completed as part of the CE-399 course at the Indian Institute of Technology Gandhinagar.
