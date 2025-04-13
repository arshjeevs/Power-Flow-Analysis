Here’s the updated README without the `.gitignore` section:

---

# Power Systems Assignment

This repository contains the implementation of various numerical methods for solving power systems problems. The main methods included are **Gauss-Seidel** and **Newton-Raphson**, which are commonly used to solve power flow problems in electrical power systems. The repository also includes the relevant data files for a 14-bus power system and other auxiliary functions.

## Files in this Repository

- **`CalculateJacobian.m`**: Code for calculating the Jacobian matrix used in the power flow analysis.
- **`Data14Bus.m`**: Contains the data for a 14-bus power system, used as input for power flow calculations.
- **`GaussSiedel.m`**: Implementation of the Gauss-Seidel method for solving the power flow equations.
- **`NewtonRaphson.m`**: Implementation of the Newton-Raphson method for solving nonlinear power flow equations.
- **`YAdmittance.m`**: Function for calculating the admittance matrix of the system.
- **`PowerSystemsAssignment_GuassSeidel.asv`**: A saved session or results file from running the Gauss-Seidel method in MATLAB.
- **`PowerSystemsAssignment_GuassSeidel.asv:Zone.Identifier`**: Metadata file automatically generated when the file was transferred or downloaded.

## Getting Started

To get started with the code in this repository, you need to have MATLAB or Octave installed on your system, as the code is written in MATLAB's `.m` file format.

### Prerequisites

- **MATLAB** or **GNU Octave** (for running the `.m` files)

### Installing MATLAB/Octave

1. **MATLAB**: [Official MATLAB download](https://www.mathworks.com/products/matlab.html) (Paid License)
2. **Octave**: [Download Octave](https://www.gnu.org/software/octave/) (Free Open-Source Alternative)

## Usage

After installing MATLAB or Octave, you can run the scripts directly from the command line in the MATLAB/Octave environment.

### Example Command to Run Gauss-Seidel Method:

```bash
>> GaussSiedel
```

This will execute the Gauss-Seidel method on the 14-bus data system provided in `Data14Bus.m`.

### Example Command to Run Newton-Raphson Method:

```bash
>> NewtonRaphson
```

This will run the Newton-Raphson power flow analysis based on the system data.

## Contributing

If you have any improvements, bug fixes, or additional features you’d like to contribute, feel free to create a pull request. For significant changes, please open an issue first to discuss what you would like to change.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- The methods and data used in this assignment are commonly found in power systems textbooks and research articles.
- Thanks to the open-source community for their contributions to MATLAB and Octave.

---

This version should work well for your repository. Let me know if you need further modifications!
