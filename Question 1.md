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

The governing system of the first-order differential equation is 

$$
\frac{d\mathbf{u}}{dt}=f(t,\mathbf{u})
$$


And where 

$$
f(t,\mathbf{u})=
\begin{bmatrix}
v_x\\
v_y\\
-\dfrac{\mu x}{r^3}-\dfrac{\rho C_D A}{2m}v v_x\\
-\dfrac{\mu y}{r^3}-\dfrac{\rho C_D A}{2m}v v_y
\end{bmatrix}
$$

With 

$$
r=\sqrt{x^2+y^2},\qquad
v=\sqrt{v_x^2+v_y^2}
$$

Here, $\mu=GM_E$ is representative of the Earth’s gravitational parameter, $\rho$ represents the atmospheric density, and $C_D$ is about the drag coefficient, $A$ is the satellite's cross-sectional area, and $m$ is the satellite's mass. NOTE: The drag model assumes that the atmosphere is stationary. 

In order to approximate the satellite position and velocity at each time step, using the RK4 method will calculate the four intermediate slopes:

$$
\begin{aligned}
k_1 &= f(t_n,\mathbf{u}_n)\\[6pt]
k_2 &= f\left(t_n+\frac{\Delta t}{2},
\mathbf{u}_n+\frac{\Delta t}{2}k_1\right)\\[6pt]
k_3 &= f\left(t_n+\frac{\Delta t}{2},
\mathbf{u}_n+\frac{\Delta t}{2}k_2\right)\\[6pt]
k_4 &= f\left(t_n+\Delta t,
\mathbf{u}_n+\Delta t k_3\right)
\end{aligned}
$$

The solution is advanced as the weighted average of these slopes are:

 $$
\boxed{
\mathbf{u}_{n+1}
=
\mathbf{u}_n+
\frac{\Delta t}{6}
(k_1+2k_2+2k_3+k_4)
}
$$

The four slopes can approximate the rate of change of the satellite’s state in different points within the time step. Their weighted average will provide an estimate of the satellite’s updated position and velocity. This will allow the trajectory to be calculated over the desired simulation period. The RK4 method has a truncation error of $O(\Delta t^5)$ and a global error of $O(\Delta t^4)$, this is assuming sufficiently smooth solutions and appropriate numerical stability. So, decreasing the time step will also be required to improve the accuracy of the computed trajectory.  RK4 will require 4 evaluations of the governing equations per time step, making it more expensive per step than Forward Euler or Heun’s method, but will be potentially more efficient for a specified accuracy. 

For evaluating numerical accuracy, a convergence study will be performed using progressively smaller time steps. Because RK4 is conditionally stable, the selected time step will also be checked to ensure that numerical errors do not significantly affect the computed orbital trajectory.

Background information: 



## Satellite Parameters and Initial Conditions

The satellite is modeled in a two-dimensional orbit around Earth. The simulation considers Earth's gravitational acceleration and atmospheric drag. The following parameters are used to establish the initial conditions and calculate the satellite's orbital trajectory.

| Parameter | Symbol | Value | Units |
|---|---|---|---|
| Earth's radius | RE | 6371 | km |
| Earth's gravitational parameter |  | 398600 | km³/s² |
| Initial satellite altitude | h0 | 400 | km |
| Satellite mass | m | 500 | kg |
| Drag coefficient | CD | 2.2 | Dimensionless |
| Satellite cross-sectional area | A | 4 | m² |
| Initial orbital velocity | v0 | 7.67 | km/s |
| Reference atmospheric density at 400 km | ref | 410-12 | kg/m³ |
| Atmospheric scale height | H | 60 | km |

The initial orbital radius is calculated by adding Earth's radius to the satellite's initial altitude:

$$
r_0 = R_E + h_0
$$

$$
r_0 = 6371 + 400 = 6771\ \text{km}
$$

Assuming an initially circular orbit, the required orbital velocity is

$$
v_0 = \sqrt{\frac{\mu}{r_0}}
$$

$$
v_0 = \sqrt{\frac{398600}{6771}}
\approx 7.67\ \text{km/s}
$$

The initial state vector is

$$
\mathbf{u}(0)=
\begin{bmatrix}
6771\\
0\\
0\\
\sqrt{398600/6771}
\end{bmatrix}
$$

where the position components are in kilometers and the velocity components are in kilometers per second.

The atmospheric density is approximated using an exponential model:

$$
\rho(h)=\rho_{ref}\exp\left(-\frac{h-h_0}{H}\right)
$$

The numerical implementation converts all dimensional quantities to SI units to maintain consistency.



CODE: 
import numpy as np
import matplotlib.pyplot as plt

mu = 3.986e14
earth_radius = 6371000
height = 400000
mass = 500
cd = 2.2
area = 4
rho0 = 4e-12
H = 60000

r0 = earth_radius + height
v0 = np.sqrt(mu / r0)
initial = np.array([r0, 0, 0, v0])

period = 2 * np.pi * np.sqrt(r0**3 / mu)


def satellite(t, u, drag_on):
    x, y, vx, vy = u

    r = np.sqrt(x**2 + y**2)
    v = np.sqrt(vx**2 + vy**2)

    ax = -mu * x / r**3
    ay = -mu * y / r**3

    if drag_on:
        altitude = r - earth_radius
        rho = rho0 * np.exp(-(altitude - height) / H)

        ax -= 0.5 * rho * cd * area / mass * v * vx
        ay -= 0.5 * rho * cd * area / mass * v * vy

    return np.array([vx, vy, ax, ay])


def rk4(t, u, dt, drag_on):
    k1 = satellite(t, u, drag_on)
    k2 = satellite(t + dt/2, u + dt*k1/2, drag_on)
    k3 = satellite(t + dt/2, u + dt*k2/2, drag_on)
    k4 = satellite(t + dt, u + dt*k3, drag_on)

    return u + dt/6 * (k1 + 2*k2 + 2*k3 + k4)


def simulate(days, dt, drag_on):
    total_time = days * 24 * 3600
    steps = int(np.ceil(total_time / dt))
    dt = total_time / steps

    times = np.linspace(0, total_time, steps + 1)
    u = np.zeros((steps + 1, 4))
    u[0] = initial

    for i in range(steps):
        u[i+1] = rk4(times[i], u[i], dt, drag_on)

        r = np.sqrt(u[i+1, 0]**2 + u[i+1, 1]**2)
        if r <= earth_radius:
            return times[:i+2], u[:i+2]

    return times, u


t1, u1 = simulate(7, 20, False)
t2, u2 = simulate(7, 20, True)

h1 = np.sqrt(u1[:, 0]**2 + u1[:, 1]**2) - earth_radius
h2 = np.sqrt(u2[:, 0]**2 + u2[:, 1]**2) - earth_radius

print("Initial speed:", round(v0/1000, 4), "km/s")
print("Orbital period:", round(period/60, 2), "minutes")
print("Final altitude without drag:", round(h1[-1]/1000, 4), "km")
print("Final altitude with drag:", round(h2[-1]/1000, 4), "km")


plt.figure()

plt.plot(t1/86400, h1/1000, label="Without drag")
plt.plot(t2/86400, h2/1000, label="With drag")

plt.xlabel("Time (days)")
plt.ylabel("Altitude (km)")
plt.title("Effect of Atmospheric Drag on Altitude")
plt.legend()
plt.grid()
plt.show()


plt.figure()

plt.plot(u1[:, 0]/1000, u1[:, 1]/1000, label="Without drag")
plt.plot(u2[:, 0]/1000, u2[:, 1]/1000, label="With drag")

theta = np.linspace(0, 2*np.pi, 200)
plt.plot(earth_radius*np.cos(theta)/1000,
         earth_radius*np.sin(theta)/1000, color="gray")

plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title("Satellite Orbital Trajectory")
plt.axis("equal")
plt.legend()
plt.show()


def energy(u):
    r = np.sqrt(u[:, 0]**2 + u[:, 1]**2)
    v2 = u[:, 2]**2 + u[:, 3]**2
    return 0.5*v2 - mu/r


e1 = energy(u1)
e2 = energy(u2)

print("Energy change without drag:", e1[-1] - e1[0])
print("Energy change with drag:", e2[-1] - e2[0])


print("\nRK4 convergence test")
print("Time step     Position error")

for dt in [80, 40, 20, 10]:
    steps = int(np.ceil(period / dt))
    step = period / steps
    u = initial.copy()

    for i in range(steps):
        u = rk4(i*step, u, step, False)

    error = np.sqrt((u[0]-initial[0])**2 +
                    (u[1]-initial[1])**2)

    print(dt, "seconds:", round(error, 6), "m")


