import json
from pathlib import Path
from time import perf_counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

MU = 3.986004418e14
RE = 6.378137e6
J2 = 1.08263e-3
MM = 4.9028e12
DIST = 384400e3
WM = 2 * np.pi / (27.32 * 86400)
MASS, AREA, CD = 100.0, 1.0, 2.2
RHO, HREF, SCALE = 4e-12, 400e3, 60e3
DURATION = 7 * 86400.0
CHECK = 600.0
RETURN_TIME = 3600.0
PHASES = np.arange(8)*np.pi/4
FIG = Path(__file__).resolve().parent / "figures"


def rhs(t, state, drag=False, j2=False, moon=False, phase=0.0):
    x, y, vx, vy = state
    r2 = x*x + y*y
    r = np.sqrt(r2)
    factor = -MU / (r2*r)
    if j2:
        factor -= 1.5 * J2 * MU * RE**2 / r**5
    ax, ay = factor*x, factor*y
    if drag:
        density = RHO * np.exp(-((r-RE)-HREF)/SCALE)
        b = -CD * AREA * density * np.hypot(vx, vy) / (2*MASS)
        ax += b*vx
        ay += b*vy
    if moon:
        mx, my = DIST*np.cos(WM*t+phase), DIST*np.sin(WM*t+phase)
        dx, dy = mx-x, my-y
        d3 = (dx*dx + dy*dy)**1.5
        ax += MM*(dx/d3 - mx/DIST**3)
        ay += MM*(dy/d3 - my/DIST**3)
    return np.array([vx, vy, ax, ay])


def rk4(f, t, u, dt):
    k1 = f(t, u)
    k2 = f(t+dt/2, u+dt*k1/2)
    k3 = f(t+dt/2, u+dt*k2/2)
    k4 = f(t+dt, u+dt*k3)
    return u + dt*(k1+2*k2+2*k3+k4)/6


def heun(f, t, u, dt):
    k1 = f(t, u)
    return u + dt*(k1 + f(t+dt, u+dt*k1))/2


def initial(altitude):
    r = RE + altitude
    return np.array([r, 0.0, 0.0, np.sqrt(MU/r)])


def integrate(f, u0, end, dt, method=rk4, controlled=False):
    n = int(round(end/dt))
    if not np.isclose(n*dt, end):
        raise ValueError("Duration must be an integer multiple of dt")
    stride = int(round(CHECK/dt))
    if controlled and not np.isclose(stride*dt, CHECK):
        raise ValueError("Maneuvers must fall on grid points")
    t = np.arange(n+1)*dt
    u = np.empty((n+1, 4))
    u[0] = u0
    dv = np.zeros(n+1)
    target = np.linalg.norm(u0[:2])
    for k in range(n):
        u[k+1] = method(f, t[k], u[k], dt)
        dv[k+1] = dv[k]
        if controlled and (k+1) % stride == 0:
            r = np.linalg.norm(u[k+1, :2])
            er = u[k+1, :2]/r
            et = np.array([-er[1], er[0]])
            desired = np.sqrt(MU/r*(1+1.5*J2*(RE/r)**2))*et - (r-target)/RETURN_TIME*er
            dv[k+1] += np.linalg.norm(desired-u[k+1, 2:])
            u[k+1, 2:] = desired
    if not np.all(np.isfinite(u)):
        raise RuntimeError("Nonfinite solution")
    if np.min(np.linalg.norm(u[:, :2], axis=1)) <= RE:
        raise RuntimeError("Trajectory intersects Earth: model no longer valid")
    return t, u, dv


def radius(u):
    return np.linalg.norm(u[:, :2], axis=1)


def j2_energy(u):
    r = radius(u)
    return np.sum(u[:, 2:]**2, axis=1)/2 - MU/r - MU*J2*RE**2/(2*r**3)


def angular_momentum(u):
    return u[:, 0]*u[:, 3]-u[:, 1]*u[:, 2]


def j2_turning_radius(r0):
    # Inner turning point of the equatorial J2 effective potential for v0 = sqrt(MU/r0).
    def V(r):
        return MU*r0/(2*r*r) - MU/r - MU*J2*RE**2/(2*r**3)
    lo, hi = r0*(1-0.02), r0*(1-1e-7)
    for _ in range(200):
        mid = (lo+hi)/2
        if V(mid) > V(r0):
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2


def independent_checks(alt, runs):
    r0 = RE + alt
    k = 1.5*J2*(RE/r0)**2
    u = runs['j2'][1]
    energy, angular = j2_energy(u), angular_momentum(u)
    checks = dict(
        j2_energy_relative=float(np.max(abs(energy/energy[0]-1))),
        j2_angular_relative=float(np.max(abs(angular/angular[0]-1))),
        j2_drop_linear_km=2*r0*k/1000,
        j2_drop_exact_km=(r0-j2_turning_radius(r0))/1000,
        j2_drop_simulated_km=float(r0-np.min(radius(u)))/1000)
    # Circular-orbit decay da/dt = -rho*sqrt(MU*a)*CD*A/m, with density following the altitude.
    c = RHO*np.sqrt(MU*r0)*CD*AREA/MASS
    checks['drag_rate_km_day'] = c*np.exp(-(alt-HREF)/SCALE)*86400/1000
    checks['drag_loss_constant_density_km'] = checks['drag_rate_km_day']*DURATION/86400
    inner = np.exp((alt-HREF)/SCALE) - c*DURATION/SCALE
    if inner > 0:
        checks['drag_loss_varying_density_km'] = (alt-HREF-SCALE*np.log(inner))/1000
    checks['drag_loss_simulated_km'] = float(r0-radius(runs['drag'][1])[-1])/1000
    return checks


def savefig(name):
    plt.savefig(FIG / (name + ".pdf"), bbox_inches="tight")
    plt.close()


def validation():
    u0 = initial(300e3)
    r0 = u0[0]
    period = 2*np.pi*np.sqrt(r0**3/MU)
    out = {}
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.2))
    for method, stages in [(heun, 2), (rk4, 4)]:
        rows = []
        for per_orbit in [80, 160, 320, 640, 1280]:
            dt = period/per_orbit
            start = perf_counter()
            t, u, _ = integrate(rhs, u0, 3*period, dt, method)
            runtime = perf_counter()-start
            exact = r0*np.column_stack((np.cos(2*np.pi*t/period),
                                       np.sin(2*np.pi*t/period)))
            error = np.max(np.linalg.norm(u[:, :2]-exact, axis=1))
            energy = np.sum(u[:, 2:]**2, axis=1)/2-MU/radius(u)
            angular = u[:, 0]*u[:, 3]-u[:, 1]*u[:, 2]
            rows.append(dict(dt=dt, error_m=error, evaluations=stages*(len(t)-1),
                             runtime_s=runtime,
                             energy_relative=float(np.max(abs(energy/energy[0]-1))),
                             angular_relative=float(np.max(abs(angular/angular[0]-1)))))
        order = np.log2(rows[-2]['error_m']/rows[-1]['error_m'])
        assert abs(order-stages) < 0.3, 'Unexpected circular-orbit convergence order'
        out[method.__name__] = rows
        axes[0].loglog([r['dt'] for r in rows], [r['error_m'] for r in rows], 'o-', label=method.__name__.upper())
        axes[1].loglog([r['evaluations'] for r in rows], [r['error_m'] for r in rows], 'o-', label=method.__name__.upper())
    axes[0].set_xlabel('Time step (s)')
    axes[1].set_xlabel('Right-hand-side evaluations')
    for ax in axes:
        ax.set_ylabel('Maximum position error (m)')
        ax.grid(True, which='both', alpha=.25)
        ax.legend()
    fig.tight_layout()
    savefig('validation')
    out['period_s'] = period
    return out


def experiments():
    results, tracks = {}, {}
    cases = {'central': {}, 'drag': {'drag': True}, 'j2': {'j2': True},
             'baseline': {'drag': True, 'j2': True},
             'moon': {'drag': True, 'j2': True, 'moon': True}}
    for label, alt, dt in [('Low', 300e3, 5.0), ('High', 35786e3, 120.0)]:
        print(f'Running {label} orbit', flush=True)
        u0 = initial(alt)
        runs, data = {}, {}
        for name, flags in cases.items():
            f = lambda t, u, flags=flags: rhs(t, u, **flags)
            start = perf_counter()
            t, u, dv = integrate(f, u0, DURATION, dt)
            runtime = perf_counter()-start
            runs[name] = (t, u, dv)
            data[name] = dict(final_altitude_km=(radius(u)[-1]-RE)/1000,
                              max_altitude_change_km=float(np.max(abs(radius(u)-u0[0]))/1000),
                              steps=len(t)-1, runtime_s=runtime)
        for name in ['baseline', 'moon']:
            flags = cases[name]
            f = lambda t, u, flags=flags: rhs(t, u, **flags)
            t, u, dv = integrate(f, u0, DURATION, dt, controlled=True)
            runs[name+'_controlled'] = (t, u, dv)
            data[name+'_controlled'] = dict(delta_v_m_s=dv[-1],
                max_altitude_error_km=float(np.max(abs(radius(u)-u0[0]))/1000))
        # Refine the two matched runs and both maintenance runs independently.
        fine = {}
        for name in ['baseline', 'moon', 'baseline_controlled', 'moon_controlled']:
            flags = cases[name.replace('_controlled', '')]
            f = lambda t, u, flags=flags: rhs(t, u, **flags)
            fine[name] = integrate(f, u0, DURATION, dt/2,
                                    controlled=name.endswith('_controlled'))
        # Compare methods with lunar gravity included, against fine-step RK4.
        f = lambda t, u: rhs(t, u, drag=True, j2=True, moon=True)
        heun_errors = []
        for comparison_dt, ref in [(dt, fine['moon'][1][::2]), (dt/2, fine['moon'][1])]:
            _, uh, _ = integrate(f, u0, DURATION, comparison_dt, method=heun)
            heun_errors.append(dict(dt_s=comparison_dt,
                max_position_difference_m=float(np.max(np.linalg.norm(uh[:, :2]-ref[:, :2], axis=1)))))
        data['heun_lunar_comparison'] = heun_errors
        ub, um = runs['baseline'][1], runs['moon'][1]
        fb, fm = fine['baseline'][1][::2], fine['moon'][1][::2]
        displacement = np.linalg.norm(um[:, :2]-ub[:, :2], axis=1)
        lunar_alt = radius(um)-radius(ub)
        convergence = {}
        for name in fine:
            coarse = runs[name]
            refined = fine[name]
            convergence[name] = dict(
                max_position_difference_m=float(np.max(np.linalg.norm(coarse[1][:, :2]-refined[1][::2, :2], axis=1))),
                max_altitude_difference_m=float(np.max(abs(radius(coarse[1])-radius(refined[1][::2])))),
                delta_v_difference_m_s=float(abs(coarse[2][-1]-refined[2][-1])))
        data['lunar'] = dict(max_displacement_km=float(max(displacement)/1000),
            max_altitude_difference_km=float(max(abs(lunar_alt))/1000),
            paired_displacement_refinement_m=float(np.max(np.linalg.norm((um[:, :2]-ub[:, :2])-(fm[:, :2]-fb[:, :2]), axis=1))),
            additional_delta_v_m_s=float(runs['moon_controlled'][2][-1]-runs['baseline_controlled'][2][-1]),
            fine_additional_delta_v_m_s=float(fine['moon_controlled'][2][-1]-fine['baseline_controlled'][2][-1]))
        # Observed order on the full lunar model from dt, dt/2 and dt/4.
        _, uq, _ = integrate(lambda t, u: rhs(t, u, drag=True, j2=True, moon=True), u0, DURATION, dt/4)
        e1 = np.max(np.linalg.norm(um[:, :2]-fine['moon'][1][::2, :2], axis=1))
        e2 = np.max(np.linalg.norm(fine['moon'][1][:, :2]-uq[::2, :2], axis=1))
        convergence['moon_order'] = dict(dt_s=[dt, dt/2, dt/4],
            difference_dt_dt2_m=float(e1), difference_dt2_dt4_m=float(e2),
            observed_order=float(np.log2(e1/e2)))
        data['convergence'] = convergence
        # Lunar maintenance cost for several initial lunar phases.
        sweep = []
        for phase in PHASES:
            f = lambda t, u, phase=phase: rhs(t, u, drag=True, j2=True, moon=True, phase=phase)
            _, u, dv = integrate(f, u0, DURATION, dt, controlled=True)
            sweep.append(dict(phase_deg=float(np.degrees(phase)), delta_v_m_s=dv[-1],
                              additional_delta_v_m_s=float(dv[-1]-runs['baseline_controlled'][2][-1]),
                              max_altitude_error_km=float(np.max(abs(radius(u)-u0[0]))/1000)))
        data['phase_sweep'] = sweep
        data['checks'] = independent_checks(alt, runs)
        data['dt_s'] = dt
        data['altitude_km'] = alt/1000
        results[label] = data
        tracks[label] = runs
    return results, tracks


def figures(tracks):
    fig, axes = plt.subplots(2, 2, figsize=(9, 6))
    for row, (label, runs) in enumerate(tracks.items()):
        t, uc, _ = runs['central']
        for name in ['drag', 'j2', 'baseline']:
            u = runs[name][1]
            axes[row, 0].plot(t/86400, (radius(u)-radius(uc))/1000, label=name)
            axes[row, 1].plot(t/86400, np.linalg.norm(u[:, :2]-uc[:, :2], axis=1)/1000, label=name)
        axes[row, 0].set_ylabel(f'{label}: altitude difference (km)')
        axes[row, 1].set_ylabel(f'{label}: position difference (km)')
        for ax in axes[row]:
            ax.set_xlabel('Time (days)')
            ax.grid(alpha=.25)
            ax.legend()
    fig.tight_layout()
    savefig('earth_effects')
    fig, axes = plt.subplots(2, 2, figsize=(9, 6))
    for row, (label, runs) in enumerate(tracks.items()):
        t, ub, _ = runs['baseline']
        um = runs['moon'][1]
        axes[row, 0].plot(t/86400, (radius(um)-radius(ub))/1000)
        axes[row, 1].plot(t/86400, np.linalg.norm(um[:, :2]-ub[:, :2], axis=1)/1000)
        axes[row, 0].set_ylabel(f'{label}: lunar altitude difference (km)')
        axes[row, 1].set_ylabel(f'{label}: lunar position difference (km)')
        for ax in axes[row]:
            ax.set_xlabel('Time (days)')
            ax.grid(alpha=.25)
    fig.tight_layout()
    savefig('lunar_effect')
    fig, axes = plt.subplots(2, 2, figsize=(9, 6))
    for row, (label, runs) in enumerate(tracks.items()):
        for name, text in [('baseline_controlled', 'without Moon'), ('moon_controlled', 'with Moon')]:
            t, u, dv = runs[name]
            axes[row, 0].plot(t/86400, (radius(u)-radius(u)[0])/1000, label=text)
            axes[row, 1].plot(t/86400, dv, label=text)
        axes[row, 0].set_ylabel(f'{label}: altitude error (km)')
        axes[row, 1].set_ylabel(f'{label}: accumulated delta-v (m/s)')
        for ax in axes[row]:
            ax.set_xlabel('Time (days)')
            ax.grid(alpha=.25)
            ax.legend()
    fig.tight_layout()
    savefig('maintenance')
    w = np.linspace(-3.3, 1, 401)[None, :] + 1j*np.linspace(-3.2, 3.2, 401)[:, None]
    fig, ax = plt.subplots(figsize=(5, 3.5))
    for amp, color, label in [(abs(1+w+w*w/2), 'C1', 'Heun'),
                              (abs(1+w+w*w/2+w**3/6+w**4/24), 'C0', 'RK4')]:
        ax.contour(w.real, w.imag, amp, levels=[1], colors=[color])
        ax.plot([], [], color=color, label=label)
    ax.axhline(0, color='gray', lw=.5)
    ax.axvline(0, color='gray', lw=.5)
    ax.set(xlabel='Re(w)', ylabel='Im(w)', title='Absolute stability boundaries, w = dt × lambda')
    ax.legend()
    ax.grid(alpha=.25)
    fig.tight_layout()
    savefig('stability')


def main():
    FIG.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size': 9, 'axes.titlesize': 10})
    report = {'validation': validation()}
    report['experiments'], tracks = experiments()
    figures(tracks)
    (FIG.parent/'results.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
