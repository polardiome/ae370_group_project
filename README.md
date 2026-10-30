# AE 370 Group Project 1: Numerical Study of a Perturbed Earth Orbit

Jiho Seo, Pranshu Bhardwaj, Adi Jha, Dmitry Moryakov. Fall 2026.

We simulate a satellite in a planar, equatorial Earth orbit with a from-scratch RK4 integrator. We then measure how atmospheric drag, Earth's oblateness ($J_2$), and the Moon's gravity change the trajectory, and how much corrective $\Delta v$ is needed to hold the orbit. The full write-up is [`report.pdf`](report.pdf); its source is [`report.tex`](report.tex).

## Questions

1. How does atmospheric drag affect the satellite's orbital trajectory over time?
2. How does Earth's gravitational field, including the $J_2$ perturbation, affect the satellite's orbital trajectory?
3. How does the Moon's gravitational influence change the satellite's trajectory, and how does it compare with drag and $J_2$ at different altitudes?
4. How do the numerical method and time step affect the predicted trajectory, particularly when lunar gravity is included?
5. How much additional corrective thrust is required to maintain the desired orbit when the Moon's gravity is included?

## Repository contents

| File | Contents |
| --- | --- |
| `simulation.py` | Equations of motion, RK4 and Heun integrators, validation checks, all experiments, and figure generation |
| `results.json` | Every number quoted in the report, written by `simulation.py` |
| `figures/` | The five report figures (stability, validation, Earth effects, lunar effect, maintenance) |
| `report.tex`, `report.pdf` | Report source and compiled PDF |
| `requirements.txt` | Python dependencies (NumPy and Matplotlib) |

## How to run

```bash
pip install -r requirements.txt
python simulation.py
```

This regenerates `results.json` and everything in `figures/` in about one minute. Then compile `report.tex` with LaTeX (for example `pdflatex report.tex`, run twice so cross-references resolve) to rebuild the PDF.

## Model

The state is $\vec u=(x,y,v_x,v_y)^T$ in an Earth-centered inertial frame, with

$$
\ddot{\vec r}=\vec a_g+\vec a_{\text{drag}}+\vec a_{J_2}+\vec a_{\text{Moon}}.
$$

- **Central gravity:** $\vec a_g=-\mu_E\vec r/r^3$.
- **Drag:** $\vec a_{\text{drag}}=-\frac{C_DA\rho(h)}{2m}|\vec v|\vec v$. It uses a stationary atmosphere and the exponential density $\rho=\rho_0e^{-(h-h_0)/H}$.
- **$J_2$:** in the equatorial plane, $\vec a_{J_2}=-\frac{3J_2\mu_ER_E^2}{2r^5}\vec r$.
- **Moon:** a third-body perturbation, $\vec a_{\text{Moon}}=\mu_M\left(\frac{\vec r_M-\vec r}{|\vec r_M-\vec r|^3}-\frac{\vec r_M}{D^3}\right)$. The Moon follows a prescribed circular orbit in the same plane as the satellite, $\vec r_M=D(\cos(\omega_Mt+\phi_0),\sin(\omega_Mt+\phi_0))$.

| Parameter | Value | Parameter | Value |
| --- | --- | --- | --- |
| $\mu_E$ | $3.986004418\times10^{14}$ m³/s² | $R_E$ | $6378.137$ km |
| $J_2$ | $1.08263\times10^{-3}$ | $\mu_M$ | $4.9028\times10^{12}$ m³/s² |
| $D$ | $384{,}400$ km | $T_M$ | $27.32$ days |
| $m$, $A$, $C_D$ | $100$ kg, $1$ m², $2.2$ | $\rho_0$, $h_0$, $H$ | $4\times10^{-12}$ kg/m³, $400$ km, $60$ km |

All calculations use SI units. Each run starts with the central-gravity circular state at one of two altitudes: **300 km** (low orbit) and **35,786 km** (geostationary altitude). Runs last 7 days, with $\Delta t=5$ s at low altitude and $120$ s at high altitude.

### Runs

| Run | Terms included |
| --- | --- |
| central | central gravity |
| drag | central gravity and drag |
| j2 | central gravity and $J_2$ |
| baseline | central gravity, drag, and $J_2$ |
| moon | baseline plus the Moon |
| baseline / moon controlled | the same runs with the maintenance rule below |

**Maintenance rule.** Every 600 s the velocity is reset to the $J_2$-corrected circular speed, plus a radial term that returns the satellite to its initial radius within 3600 s. The cost is the sum of the velocity jumps, $\Delta v_{\text{total}}=\sum_j|\vec v^+-\vec v^-|$.

## Validation

- **Exact two-body solution:** the observed convergence orders are 2.03 for Heun and 4.10 for RK4. RK4 conserves energy and angular momentum to about $10^{-12}$.
- **$J_2$ checks:**
  - The $J_2$-only runs conserve the $J_2$ energy and the angular momentum to about $10^{-10}$.
  - The radial drop matches the analytic value $3J_2R_E^2/r_0$: 19.785 km at low altitude and 3.134 km at high altitude.
- **Drag check:** an analytic decay estimate predicts a 7-day altitude loss of 16.649 km; the simulation gives 16.643 km.
- **Full model:** refining the step to $\Delta t/2$ and $\Delta t/4$ gives observed RK4 orders of 4.28 (low) and 4.21 (high).

## Key results (7 days)

| Quantity | Low orbit (300 km) | High orbit (35,786 km) |
| --- | --- | --- |
| Altitude loss, drag only | 16.643 km | negligible |
| Maximum radial drop, $J_2$ only | 19.785 km | 3.134 km |
| Final altitude, baseline | 270.427 km | 35,785.988 km |
| Maximum lunar position separation | 0.497 km | 55.457 km |
| Maximum lunar altitude difference | 5.5 m | 3.660 km |
| Maintenance $\Delta v$ without / with Moon | 20.837 / 20.739 m/s | 0.058 / 3.482 m/s |
| Lunar change in $\Delta v$ over 8 lunar phases | −0.0987 to −0.0971 m/s | +3.396 to +3.465 m/s |

- **Low orbit:** drag and $J_2$ dominate, and the Moon's effect is small.
- **High orbit:** drag vanishes, and the Moon's radial effect is comparable to $J_2$.
- **Heun at the same step sizes:** its errors (50–530 km) are larger than the lunar signal, so the choice of method matters.
- **Maintenance costs:** these depend on the control rule. At high altitude, much of the extra cost comes from the rule cancelling the Moon's periodic radial oscillation.

## Assumptions and limitations

- The model is planar and equatorial. The Moon's orbit is circular and coplanar, while the real one is inclined 18.3–28.6° to Earth's equator, so out-of-plane effects are not modeled.
- The Moon's position is prescribed, not integrated.
- The Sun, solar radiation pressure, higher gravity harmonics, atmospheric rotation, and mass loss are neglected.
- One exponential density law is used at every altitude.
- Only one maintenance rule is used, and its $\Delta v$ is not an optimal station-keeping cost.
