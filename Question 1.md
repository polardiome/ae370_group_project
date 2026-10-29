How does atmospheric drag affect the satellite’s orbital trajectory in time?

What to say in the premise:

Atmospheric drag is an important factor that affects satellites in low Earth orbit. Although atmospheric density is relatively low at high altitudes, drag can remove mechanical energy from the satellite, causing its orbit to decay. 

What to say in the method:

The fourth-order Runge-Kutta (RK4) method is selected to numerically solve the system of ordinary differential equations describing how the satellite moves under the influence of Earth’s gravitational field and atmospheric drag. RK4 is an explicit one-step method that often uses four intermediate evaluations of the derivative function, allowing it to approximate the solution at the next time step. We chose this method because it provides having fourth order accuracy while avoiding the need to solve a nonlinear system of equations at each time step. The motion is modeled using a state vector containing the position and velocity in two dimensions: 


$$
\mathbf{u}(t)=
\begin{bmatrix}
x(t)\\
y(t)\\
v_x(t)\\
v_y(t)
\end{bmatrix}
$$

The 
