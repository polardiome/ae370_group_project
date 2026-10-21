# AE370 Group Project 1

## The Objective
Satellites orbiting Earth are influenced by atmospheric drag and Earth’s nonspherical gravity, meaning corrective thruster maneuvers are required to maintain a desired orbit. We will study a satellite’s trajectory across time and the level of corrective thrust involved in maintaining its orbit. 

## The Questions
1. How does atmospheric drag affect the satellite’s orbital trajectory in time?

2. How does Earth’s gravitational field affect a satellite's orbital trajectory?

3. How does choosing the numerical method time step affect the predicted trajectory?

## The Equation
We will model the acceleration of the satellite using a second-order vector ODE

$$
\ddot{\vec{r}} = \vec{a}_{g} + \vec{a}_{\text{drag}} + \vec{a}_{J_2} + \vec{a}_{\text{thrust}},
$$

where each term represents a physical force. For gravity, Newton's law of gravitation gives

$$
\vec{a}_g=-\frac{\mu}{r^3}\vec{r},
$$

where =GM is Earth's gravitational parameter and $r=\|\vec r\|$ is the distance from Earth's center.

Atmospheric drag is modeled by

$$
\vec{a}_{\text{drag}}=-\frac{C_D A\rho(h)}{2m}\|\vec{v}_{\text{rel}}\|\vec{v}_{\text{rel}},
$$

where CD is the drag coefficient, A is the cross-sectional area, m is the satellite's mass, and vrel is its velocity relative to the atmosphere. Atmospheric density is approximated by

$$
\rho(h)=\rho_0 e^{-(h-h_0)/H},
$$

where h=r-RE is altitude, 0 is the reference density at altitude h0, and H is the scale height.

Earth's oblateness introduces the J2 perturbation:

$$
\vec{a}_{J_2}=\frac{3J_2\mu R_E^2}{2r^5}
\begin{bmatrix}
x(5z^2/r^2-1)\\
y(5z^2/r^2-1)\\
z(5z^2/r^2-3)
\end{bmatrix},
$$

where J21.0826310-3 accounts for Earth's equatorial bulge.

Finally, corrective thrust is represented by

$$
\vec{a}_{\text{thrust}}=\frac{\vec F_{\text{thrust}}(t)}{m},
$$

where Fthrust(t) is the applied thrust force.

Together, these terms form a coupled nonlinear ODE system with time as its only independent variable. The system generally requires numerical integration.

$$
\ddot{x} =
-\frac{\mu x}{(x^2+y^2)^{3/2}}
-\frac{C_D A \rho(h)}{2m}
v_{\mathrm{rel}}v_{\mathrm{rel},x}
+u_x(t)
$$

$$
\ddot{y} =
-\frac{\mu y}{(x^2+y^2)^{3/2}}
-\frac{C_D A \rho(h)}{2m}
v_{\mathrm{rel}}v_{\mathrm{rel},y}
+u_y(t)
$$

$$
\ddot{\mathbf{r}} =
-\frac{\mu}{r^3}\mathbf{r}
-\frac{C_D A \rho(h)}{2m}
\left|\mathbf{v}_{\mathrm{rel}}\right|
\mathbf{v}_{\mathrm{rel}}
+\mathbf{u}(t)
$$

https://docs.google.com/document/d/12XFq34vSOm7k19JTvr-Zwvom3pJGrt-xf3V_f78Y1HE/edit?tab=t.0
