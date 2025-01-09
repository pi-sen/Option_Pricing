# Options Pricing Project

## Overview
This project implements various numerical methods for pricing European, Asian, and American options. The implementation is done in R and includes several sophisticated pricing techniques and variance reduction methods.

## Team Members
- Piyush Sen
- Aastha Shah
- Srikari Rallabandi

## Project Structure

### 1. European Option Pricing
- Black-Scholes exact solution implementation for benchmarking
- Monte Carlo simulation with:
  - Standard implementation for put options
  - Variance reduction using antithetic variates
  - Comparative analysis with error metrics
  - Visualization of convergence and standard error reduction

Key Features:
- Exact Black-Scholes pricing for European put options
- Monte Carlo path generation with GBM
- Antithetic variates for variance reduction
- Standard error and confidence interval calculations
- Performance benchmarking across sample sizes

### 2. Asian Option Pricing
Three distinct approaches implemented:

1. Standard Monte Carlo:
   - Direct simulation of arithmetic average prices
   - Basic GBM path generation
   - Standard error calculation

2. Control Variate Method:
   - Uses geometric Asian option as control
   - Analytical solution for geometric Asian options
   - Optimal control coefficient through regression
   - Significant variance reduction

3. Moment Matching:
   - Theoretical moment calculations
   - Path adjustment to match moments
   - Preservation of distributional properties

Features:
- Comprehensive comparison across methods
- Variance reduction effectiveness analysis
- Computational efficiency metrics
- Convergence visualization
- Standard error reduction tracking

### 3. American Option Pricing
Multiple pricing approaches implemented:

1. Least Squares Monte Carlo (LSM):
   - Pre-optimized implementation
   - Post-optimized version
   - Basis function selection
   - Early exercise boundary determination

2. Binomial Tree:
   - Standard implementation
   - Parameter optimization
   - Early exercise handling

3. BBSR Method:
   - Richardson extrapolation
   - Multiple time step analysis
   - Convergence improvement

4. Finite Difference Method (FDM):
   - Implicit scheme implementation
   - Grid resolution analysis
   - Boundary condition handling


## Implementation Details

### European Options
- Black-Scholes exact solution for benchmarking
- Monte Carlo with antithetic variates
- Sample sizes: 100 to 100,000
- Standard error and confidence interval calculations

### Asian Options
- Arithmetic average price options
- Three methods compared: Standard MC, Control Variate, Moment Matching
- Sample sizes: 1,000 to 1,024,000
- Monitoring points: 50

### American Options
- Multiple methods implemented
- LSM: 10,000 paths, 50 time steps
- Binomial Tree: 500 steps
- BBSR: 1,000 steps for extrapolation
- FDM: 100 grid points

## Key Results

### European Options
- Antithetic variates reduced variance by approximately 80%
- Convergence to Black-Scholes price achieved
- Computational efficiency improved with vectorization

### Asian Options
1. Control Variate Method:
   - Best performance in variance reduction
   - Consistent convergence across sample sizes
   - Standard error reduction by factor of ~10

2. Moment Matching:
   - Stable convergence behavior
   - Computational efficiency comparable to control variate
   - Theoretical moment preservation

### American Options
- LSM post-optimization showed best balance of accuracy and speed
- BBSR provided improved convergence over standard binomial tree
- FDM showed stable convergence with increasing grid points


## Visualization
The project includes multiple visualizations:
- Convergence plots for all methods
- Standard error comparison charts
- Computational time analysis
- Method comparison graphs

## Performance Metrics
- Pricing accuracy
- Computational efficiency
- Standard error reduction
- Convergence rates
- Memory usage

## Future Improvements
1. Parallel processing implementation
2. Additional variance reduction techniques
3. GPU acceleration for Monte Carlo simulations
4. Interactive visualization dashboard
5. Real-time market data integration

## Contributing
This project was completed as part of an academic course. While it's not open for direct contributions, feedback and suggestions are welcome.

## License
This project is for academic purposes only. Please cite appropriately if used in academic work.

## Acknowledgments
Special thanks to our team members:
- Aastha Shah
- Srikari Rallabandi
for their valuable contributions to this project.
