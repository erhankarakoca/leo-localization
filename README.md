# Single-LEO Snapshots Based Self-UE Localization with TDOA

## Overview

This project implements a Time Difference of Arrival (TDOA) based localization algorithm to estimate the position of a User Equipment (UE) using a constellation of Low Earth Orbit (LEO) satellites. The position of the UE is estimated based on TDOA measurements between single LEO satellite pairs at different times. The problem is formulated as an optimization problem, which is solved using nonlinear least squares optimization techniques to minimize the residuals between the measured and computed TDOA values.

### Key Components:

1. **Satellite Scenario Creation**: The simulation begins by creating a satellite scenario using a provided Two-Line Element (TLE) file. The scenario is defined over a given time window, and the satellites' positions and velocities are propagated based on the TLE data.

2. **TDOA Computation**: The Time of Arrival (TOA) for each satellite to the ground station (UE) is calculated. The TDOA values are computed for all pairs of satellites.

3. **Optimization Problem**: The optimization problem is formulated to minimize the difference between the measured TDOA values and those computed from the distances between the UE and the satellites. This is achieved through the use of a nonlinear least squares approach.

4. **Visualization**: The positions of the satellites, the estimated UE position, and the actual UE position are visualized in 3D space along with the hyperboloids representing the TDOA constraints.

## Mathematical Formulation

### Problem Setup:

Let the UE's position be represented as a vector:

$$
\mathbf{p} = \begin{bmatrix} p_1 \\ p_2 \\ p_3 \end{bmatrix}
$$

where $\( p_1, p_2, p_3 \)$ are the coordinates of the UE in 3D space.

Let the positions of the satellites be stored in a matrix $\( \mathbf{S} \)$:

$$
\mathbf{S} = \begin{bmatrix} 
x_1 & y_1 & z_1 \\
x_2 & y_2 & z_2 \\
\vdots & \vdots & \vdots \\
x_N & y_N & z_N
\end{bmatrix}
$$

where $\( \mathbf{S}_i = [x_i, y_i, z_i]^T \)$ represents the position of the $\( i \)-th$ satellite.

Let the TDOA measurements between the satellite pairs $\( i \) and \( j \)$ be stored in a vector $\( \mathbf{TDOA} \)$:

$$
\mathbf{TDOA} = \begin{bmatrix}
TDOA_{12} \\
TDOA_{13} \\
\vdots \\
TDOA_{ij}
\end{bmatrix}
$$

The distances from the UE to the satellites are:

$$
d_i(\mathbf{p}) = \| \mathbf{p} - \mathbf{S}_i \| = \sqrt{(p_1 - x_i)^2 + (p_2 - y_i)^2 + (p_3 - z_i)^2}
$$

The TDOA values are computed as the difference in the distances between pairs of satellites:

$$
TDOA_{ij} = \frac{1}{c} \left( \| \mathbf{p} - \mathbf{S}_i \| - \| \mathbf{p} - \mathbf{S}_j \| \right)
$$

where $\( c \)$ is the speed of light.

### Objective Function:

The objective function for the optimization problem is to minimize the residual between the measured TDOA and the computed TDOA:

```math
f(\mathbf{p})=\sum_{(i, j)}\left|\left(\frac{1}{c}\left(\left\|\mathbf{p}-\mathbf{S}_i\right\|-\left\|\mathbf{p}-\mathbf{S}_j\right\|\right)\right)-T D O A_{i j}\right|^2
```

or equivalently:

$$
f(\mathbf{p}) = \|\mathbf{r}(\mathbf{p})\|^2
$$

where $\( \mathbf{r}(\mathbf{p}) \)$ is the residual vector between the measured and calculated TDOA values.

### Matrix Formulation:

1. **Distance Vector**: The distances $\( d_i(\mathbf{p}) \)$ for all satellites can be written as:

$$
\mathbf{d}(\mathbf{p}) = \begin{bmatrix}
\| \mathbf{p} - \mathbf{S}_1 \| \\
\| \mathbf{p} - \mathbf{S}_2 \| \\
\vdots \\
\| \mathbf{p} - \mathbf{S}_N \|
\end{bmatrix}
$$

2. **TDOA Computation**: The TDOA values for each pair $\( (i, j) \)$ can be written as:

$$
\mathbf{TDOA}_{ij} = \frac{1}{c} \left( \mathbf{d}_i(\mathbf{p}) - \mathbf{d}_j(\mathbf{p}) \right)
$$

3. **Objective Function**: The optimization function in matrix form is:

$$
f(\mathbf{p}) = \|\mathbf{r}(\mathbf{p})\|^2
$$

where the vector $\( \mathbf{r}(\mathbf{p}) \)$ is the residual between the measured and computed TDOA values.

### Solving the Optimization Problem:

Using least squares, the optimization problem becomes:

$$
\mathbf{p}^{*} = \arg\min_{\mathbf{p}} \|\mathbf{r}(\mathbf{p})\|^2
$$

This is a nonlinear least squares problem that can be solved using methods such as **Gauss-Newton** or **Levenberg-Marquardt**. In MATLAB, you can use the `lsqnonlin` function to minimize this objective.
