# AE370 Group Project 1

## The Objective

Satellites orbiting Earth are influenced by atmospheric drag, Earth's nonspherical gravitational field, and gravitational perturbations from the Moon. These effects can change a satellite's trajectory over time and may require corrective thruster maneuvers to maintain a desired orbit. We will numerically study a satellite's trajectory, investigate how lunar gravity changes the predicted motion, and estimate the corrective thrust needed to maintain its orbit.

## The Questions

1. How does atmospheric drag affect the satellite's orbital trajectory over time?
2. How does Earth's gravitational field, including the $J_2$ perturbation, affect the satellite's orbital trajectory?
3. How does the Moon's gravitational influence change the satellite's trajectory, and how does that influence compare with atmospheric drag and $J_2$ at different altitudes?
4. How does choosing the numerical method and time step affect the predicted trajectory, particularly when lunar gravity is included?
5. How much additional corrective thrust is required to maintain the desired orbit when the Moon's gravity is included?

## The Equation

We will model the satellite's acceleration using a second-order vector ordinary differential equation (ODE):

$$
\ddot{\vec r}=\vec a_g+\vec a_{\text{drag}}+\vec a_{J_2}+\vec a_{\text{Moon}}+\vec a_{\text{thrust}},
$$

where each term represents an acceleration due to a physical effect. The position $\vec r$ and velocity $\vec v=\dot{\vec r}$ are measured in an Earth-centered frame.

### Earth's Central Gravity

Newton's law of gravitation gives

$$
\vec a_g=-\frac{\mu_E}{r^3}\vec r,
$$

where $\mu_E=GM_E$ is Earth's gravitational parameter and $r=|\vec r|$ is the distance from Earth's center.

### Atmospheric Drag

Atmospheric drag is modeled by

$$
\vec a_{\text{drag}}=-\frac{C_D A\rho(h)}{2m}|\vec v_{\text{rel}}|\vec v_{\text{rel}},
$$

where $C_D$ is the drag coefficient, $A$ is the cross-sectional area, $m$ is the satellite's mass, and $\vec v_{\text{rel}}$ is its velocity relative to the atmosphere. Atmospheric density is approximated by

$$
\rho(h)=\rho_0 e^{-(h-h_0)/H},
$$

where $h=r-R_E$ is altitude, $\rho_0$ is the reference density at altitude $h_0$, and $H$ is the scale height. For a first approximation, the atmosphere may be treated as stationary in the Earth-centered frame; a more detailed model can account for atmospheric co-rotation.

### Earth's $J_2$ Perturbation

Earth's oblateness introduces the $J_2$ perturbation:

$$
\vec a_{J_2}=\frac{3J_2\mu_E R_E^2}{2r^5}
\begin{bmatrix}
x(5z^2/r^2-1)\\
y(5z^2/r^2-1)\\
z(5z^2/r^2-3)
\end{bmatrix},
$$

where $J_2\approx1.08263\times10^{-3}$ accounts for Earth's equatorial bulge. This expression assumes the $z$-axis is aligned with Earth's rotation axis.

### Lunar Gravitational Perturbation

The Moon introduces an additional, time-dependent gravitational acceleration. In an Earth-centered frame, the *third-body* perturbation is

$$
\boxed{\vec a_{\text{Moon}}=\mu_M\left(\frac{\vec r_M(t)-\vec r}{|\vec r_M(t)-\vec r|^3}-\frac{\vec r_M(t)}{|\vec r_M(t)|^3}\right)},
$$

where $\mu_M\approx4902.8\ \mathrm{km^3/s^2}$ is the Moon's gravitational parameter and $\vec r_M(t)$ is the position of the Moon relative to Earth's center. The first term represents the Moon's pull on the satellite; the second subtracts the Moon's acceleration of Earth's center.

For a simplified two-dimensional model, assume that the Moon moves in a circular orbit in the same plane as the satellite:

$$
\vec r_M(t)=D
\begin{bmatrix}
\cos(\omega_M t+\phi_0)\\
\sin(\omega_M t+\phi_0)
\end{bmatrix},
\qquad \omega_M=\frac{2\pi}{T_M},
$$

with the approximate values

| Parameter | Value |
| --- | --- |
| Mean Earth–Moon distance, $D$ | $384{,}400\ \mathrm{km}$ |
| Lunar orbital period, $T_M$ | $27.32\ \mathrm{days}$ |
| Lunar angular speed, $\omega_M$ | $2.662\times10^{-6}\ \mathrm{rad/s}$ |
| Initial lunar phase, $\phi_0$ | Chosen initial angle |

This is a simplified prescribed Moon orbit, not a numerical solution for the Moon's motion. It neglects the Moon's orbital eccentricity and inclination. A more realistic three-dimensional model could use a time-dependent lunar position from ephemeris data.

### Corrective Thrust

Corrective thrust is represented by

$$
\vec a_{\text{thrust}}=\frac{\vec F_{\text{thrust}}(t)}{m}=\vec u(t),
$$

where $\vec F_{\text{thrust}}(t)$ is the applied thrust force and $\vec u(t)$ is the control acceleration. To quantify orbit maintenance, we will specify a consistent control rule (for example, applying a correction when orbital altitude or another chosen orbital element deviates beyond a tolerance).

## The Two-Dimensional Numerical Model

For a first implementation, we will use the planar state vector

$$
\vec y(t)=\begin{bmatrix}x&y&v_x&v_y\end{bmatrix}^{T},
\qquad
\dot{\vec y}(t)=\begin{bmatrix}v_x&v_y&\ddot x&\ddot y\end{bmatrix}^{T}.
$$

Let $r=\sqrt{x^2+y^2}$, $\vec r_M=(x_M,y_M)^T$, and $d=\sqrt{(x_M-x)^2+(y_M-y)^2}$. In the planar equatorial approximation ($z=0$), the components become

$$
\begin{aligned}
\ddot x={}&-\frac{\mu_E x}{r^3}
-\frac{C_D A\rho(h)}{2m}|\vec v_{\text{rel}}|v_{\text{rel},x}
-\frac{3J_2\mu_E R_E^2}{2r^5}x\\
&+\mu_M\left(\frac{x_M-x}{d^3}-\frac{x_M}{D^3}\right)+u_x(t),
\end{aligned}
$$

$$
\begin{aligned}
\ddot y={}&-\frac{\mu_E y}{r^3}
-\frac{C_D A\rho(h)}{2m}|\vec v_{\text{rel}}|v_{\text{rel},y}
-\frac{3J_2\mu_E R_E^2}{2r^5}y\\
&+\mu_M\left(\frac{y_M-y}{d^3}-\frac{y_M}{D^3}\right)+u_y(t).
\end{aligned}
$$

These coupled, nonlinear, time-dependent ODEs can be integrated numerically. A three-dimensional extension would instead use the six-component state $(x,y,z,v_x,v_y,v_z)$ and the full vector $J_2$ and lunar-acceleration formulas above.

**Unit consistency:** The lunar parameters are listed in kilometers and seconds. Any simulation using SI units must convert them to meters and seconds, and all other quantities (including density, area, and mass) must be expressed in a compatible system.

## The Numerical Method

We will solve the first-order ODE system using a numerical integration method such as fourth-order Runge–Kutta (RK4). We will compare time steps to evaluate numerical convergence and separate numerical integration error from physical trajectory changes. For a selected initial orbit, we will simulate the following cases:

1. **Model A — Baseline:** Earth's central gravity, atmospheric drag, and $J_2$, without lunar gravity or corrective thrust.
2. **Model B — Lunar perturbation:** All Model A terms plus lunar third-body gravity, without corrective thrust.
3. **Model C — Orbit maintenance:** All Model B terms plus a defined corrective-thrust rule.

We will compare orbital altitude, trajectory displacement, and orbital elements such as eccentricity where appropriate. For Model C, we will also compare accumulated corrective $\Delta v$, calculated from the applied control acceleration as

$$
\Delta v_{\text{total}}=\int_0^{t_f}|\vec u(t)|\,dt
$$

for a continuous-thrust model (or the sum of the magnitudes of individual velocity corrections for impulsive maneuvers).

To isolate the effect of the Moon, we will compare Models A and B with identical initial conditions and numerical settings. We may repeat this comparison at multiple orbital altitudes: atmospheric drag is generally more important in low Earth orbit, while the relative importance of lunar gravitational perturbations increases for higher Earth orbits. Model C will show the orbit-maintenance cost under the chosen control rule.

## The Assumptions

- The initial simulation is two-dimensional and assumes an equatorial Earth orbit and a coplanar, circular lunar orbit.
- The Moon's position is prescribed as a function of time; the Moon's orbit is not integrated as an additional dynamical body.
- The effects of the Sun, solar radiation pressure, and other higher-order gravitational harmonics are neglected.
- The atmospheric density follows an exponential approximation; real density varies with solar activity and other conditions.
- A fixed reference orbit and control rule will be chosen before computing corrective thrust.
- An optional extension is a three-dimensional orbit model with lunar inclination and more realistic Moon positions.

This formulation extends the original project without changing its main emphasis on orbital dynamics, numerical integration, time-step sensitivity, and corrective thrust.
