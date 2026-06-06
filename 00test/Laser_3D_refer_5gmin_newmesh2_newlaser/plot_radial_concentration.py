#!/usr/bin/env python3
"""
Radial particle mass-flux density plot (trajectory-based).

Method (replaces the old instantaneous-slice approach):
  1. Build full particle trajectories from processorN/<time>/lagrangian/
     using (origProcId, origId) as the persistent particle key.
  2. For each x-plane, detect exact crossing events via cubic Lagrange
     interpolation (falls back to quadratic/linear when needed).
  3. Bin crossing masses by radial position r = sqrt(y^2 + z^2) into
     annular rings.
  4. Divide by ring area and statistical time duration to get mass flux
     density [kg/m^2/s].

Usage: python plot_radial_concentration.py
"""

import os, sys, glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Reuse trajectory/crossing infrastructure from analyze_crossing_of
from analyze_crossing_of import build_trajectories, find_plane_crossings

# ---------------------------------------------------------------------------
# Configuration  (edit here, or leave as None to be prompted)
# ---------------------------------------------------------------------------

CASE_DIR   = "."
CLOUD_NAME = "kinematicCloud"

# Coordinate system (all in metres after internal conversion):
#   x-axis    : axial (laser propagation direction, +x downstream)
#   Origin    : nozzle exit / laser origin
#   powderInlet : x = -0.0627 m
#   gasInlet    : x ~  0 m
#   focal pt    : x = +0.0085 m  (typical analysis plane)
#   outlet      : x = +0.070 m

T_START   = None   # start time [s]        -- prompted if None
T_END     = None   # end time   [s]        -- prompted if None
X_TARGETS = None   # list of x planes [m]  -- prompted if None (input in mm)
R_MAX     = None   # max radial extent [m] -- prompted if None (input in mm)
N_RBINS   = 60     # number of radial bins

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def prompt_float(msg, default=None):
    s = input(f"{msg}" + (f" [{default}]: " if default is not None else ": ")).strip()
    return float(s) if s else float(default)

def prompt_floatlist(msg, default=None):
    s = input(f"{msg}" + (f" [{default}]: " if default is not None else ": ")).strip()
    if not s and default is not None:
        s = default
    return [float(v.strip()) for v in s.replace(',', ' ').split()]

def prompt_int(msg, default=None):
    s = input(f"{msg}" + (f" [{default}]: " if default is not None else ": ")).strip()
    return int(s) if s else int(default)

def _try_float(s):
    try: float(s); return True
    except: return False

# ---------------------------------------------------------------------------
# Discover processors and available time steps
# ---------------------------------------------------------------------------

proc_dirs = sorted(glob.glob(os.path.join(CASE_DIR, 'processor[0-9]*')))
if not proc_dirs:
    sys.exit(f"No processorN directories found in {CASE_DIR}")

n_procs = len(proc_dirs)
print(f"Found {n_procs} processor directories: processor0 ... processor{n_procs-1}")

time_set = set()
for pdir in proc_dirs:
    for name in os.listdir(pdir):
        lag = os.path.join(pdir, name, 'lagrangian', CLOUD_NAME, 'positions')
        if os.path.exists(lag):
            time_set.add(name)

time_entries = sorted([(float(s), s) for s in time_set if _try_float(s)])
if not time_entries:
    sys.exit("No time steps with particle data found.")

print(f"Available time steps : {len(time_entries)}")
print(f"Time range           : {time_entries[0][0]:.5f}s - {time_entries[-1][0]:.5f}s\n")

# ---------------------------------------------------------------------------
# Interactive parameter input
# ---------------------------------------------------------------------------

if T_START is None:
    T_START = prompt_float("Time range start [s]", 0.0)
if T_END is None:
    T_END = prompt_float("Time range end   [s]", time_entries[-1][0])
if X_TARGETS is None:
    X_TARGETS = [v / 1e3 for v in prompt_floatlist(
        "x-plane positions [mm], space/comma-separated", "5 7 8.5"
    )]
if R_MAX is None:
    R_MAX = prompt_float("r_max (max radial position) [mm]", 5.0) / 1e3
N_RBINS = prompt_int("Number of radial bins", N_RBINS)

t_duration = T_END - T_START
if t_duration <= 0:
    sys.exit("T_END must be greater than T_START.")

print(f"\nSettings:")
print(f"  Time range   : [{T_START:.5f}, {T_END:.5f}] s  (duration {t_duration:.5f} s)")
print(f"  x-planes     : {[round(x * 1e3, 3) for x in X_TARGETS]} mm")
print(f"  r_max        : {R_MAX * 1e3:.3f} mm")
print(f"  radial bins  : {N_RBINS}")

# ---------------------------------------------------------------------------
# Build particle trajectories (one pass, shared across all x-planes)
# ---------------------------------------------------------------------------

trajectories = build_trajectories(
    CASE_DIR, time_entries, n_procs, T_START, T_END, cloud=CLOUD_NAME
)

# ---------------------------------------------------------------------------
# For each x-plane: detect crossings -> bin by r -> compute flux density
#
# flux(r) = sum(m_crossing_in_bin) / ring_area / t_duration   [kg/m^2/s]
# ---------------------------------------------------------------------------

r_edges    = np.linspace(0, R_MAX, N_RBINS + 1)
r_centers  = 0.5 * (r_edges[:-1] + r_edges[1:])
ring_areas = np.pi * (r_edges[1:]**2 - r_edges[:-1]**2)

flux_profiles = np.zeros((len(X_TARGETS), N_RBINS))
n_crossings   = np.zeros(len(X_TARGETS), dtype=int)
total_mdot    = np.zeros(len(X_TARGETS))   # [kg/s]

for xi, x_tgt in enumerate(X_TARGETS):
    crossings = find_plane_crossings(trajectories, x_tgt, direction='positive')
    n_crossings[xi] = len(crossings)
    if not crossings:
        print(f"  x={x_tgt*1e3:.2f} mm : no crossing events found")
        continue

    # find_plane_crossings stores y and z in mm -> convert back to m
    y_c = np.array([c['y'] for c in crossings]) * 1e-3   # mm -> m
    z_c = np.array([c['z'] for c in crossings]) * 1e-3
    m_c = np.array([c['mass'] for c in crossings])        # kg
    r_c = np.sqrt(y_c**2 + z_c**2)

    mass_bin, _ = np.histogram(r_c, bins=r_edges, weights=m_c)

    with np.errstate(invalid='ignore', divide='ignore'):
        flux_profiles[xi] = np.where(
            ring_areas > 0,
            mass_bin / ring_areas / t_duration,
            0.0
        )

    total_mdot[xi] = np.sum(mass_bin) / t_duration   # kg/s

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(7, 5))

markers = ['o', 's', '^', 'D', 'v', 'p']
for xi, x_tgt in enumerate(X_TARGETS):
    label = f"x={x_tgt*1e3:.1f} mm"
    ax.plot(r_centers * 1e3, flux_profiles[xi],
            marker=markers[xi % len(markers)],
            markersize=4, linewidth=1.2, label=label)

ax.set_xlabel("Radial Position (mm)", fontsize=12)
ax.set_ylabel("Mass Flux Density (kg/m\u00b2/s)", fontsize=12)
ax.set_title(
    f"Particle mass flux density at x-planes\n"
    f"t=[{T_START:.4f}, {T_END:.4f}] s  (trajectory interpolation)",
    fontsize=11
)
ax.legend(fontsize=10)
ax.set_xlim(0, R_MAX * 1e3)
ax.set_ylim(bottom=0)
ax.grid(True, linestyle='--', alpha=0.4)

out_file = "radial_concentration.png"
fig.savefig(out_file, dpi=150, bbox_inches='tight')
print(f"\nPlot saved to {out_file}")

# ---------------------------------------------------------------------------
# Summary tables
# ---------------------------------------------------------------------------

print(f"\n{'x (mm)':>10} | {'crossings':>10} | {'peak (kg/m\u00b2/s)':>16} | "
      f"{'r_peak (mm)':>12} | {'mdot (g/min)':>12}")
print("-" * 70)
for xi, x_tgt in enumerate(X_TARGETS):
    pk = np.argmax(flux_profiles[xi])
    print(f"{x_tgt*1e3:>10.2f} | {n_crossings[xi]:>10} | "
          f"{flux_profiles[xi][pk]:>16.4f} | "
          f"{r_centers[pk]*1e3:>12.3f} | "
          f"{total_mdot[xi]*1e3*60:>12.4f}")

print(f"\nNote: mdot = integral of flux x ring_area over all bins  [kg/s -> g/min]")
