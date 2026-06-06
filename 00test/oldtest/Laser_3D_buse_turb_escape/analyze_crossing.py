#!/usr/bin/env python3
"""
通过粒子轨迹插值计算粒子穿越观察面的精确位置和时间
Track particle trajectories and interpolate exact plane-crossing events.

算法流程:
1. 读取所有时间步的 VTP 文件，以 (origProcId, origId) 为键构建每个粒子的完整轨迹
2. 对每个粒子轨迹，遍历相邻时间步对，检测 x 坐标是否跨越观察面 x=x_target
3. 对穿越事件作线性插值，精确求出穿越时刻、位置 (y, z) 及粒子属性
4. 统计并可视化：YOZ 平面质量分布、径向分布、粒子轨迹图、温度分布
"""

import json
import base64
import xml.etree.ElementTree as ET
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import subprocess


# ---------------------------------------------------------------------------
# VTP parsing
# ---------------------------------------------------------------------------

def parse_vtp(filename):
    """
    解析 vtkPolyData XML (.vtp) 文件，无需 vtk 库。
    返回 (points [N,3], fields dict)，其中 fields 包含 origId/origProcId 等字段。
    """
    tree = ET.parse(filename)
    root = tree.getroot()

    vtk_root = root if root.tag in ('VTKFile', '{VTK}VTKFile') else root
    hdr_type = vtk_root.get('header_type', 'UInt32')
    hdr_size = 8 if hdr_type == 'UInt64' else 4
    byte_order = vtk_root.get('byte_order', 'LittleEndian')
    endian = '<' if 'Little' in byte_order else '>'

    piece = root.find('.//{*}Piece') or root.find('.//Piece')
    if piece is None:
        return None, None

    n_points = int(piece.get('NumberOfPoints', 0))
    if n_points == 0:
        return None, None

    dtype_map = {
        'Float32': endian + 'f4', 'Float64': endian + 'f8',
        'Int32':   endian + 'i4', 'Int64':   endian + 'i8',
        'UInt8':   'u1',          'UInt32':   endian + 'u4', 'UInt64': endian + 'u8'
    }

    def decode_array(da_elem):
        dtype_str = dtype_map.get(da_elem.get('type', 'Float32'), endian + 'f4')
        ncomp = int(da_elem.get('NumberOfComponents', 1))
        fmt = da_elem.get('format', 'ascii')
        text = (da_elem.text or '').strip()
        if fmt == 'ascii':
            arr = np.fromstring(text, dtype=np.dtype(dtype_str), sep=' ')
        else:
            raw = base64.b64decode(text)
            arr = np.frombuffer(raw[hdr_size:], dtype=np.dtype(dtype_str))
        if ncomp > 1 and arr.size % ncomp == 0:
            arr = arr.reshape(-1, ncomp)
        return arr

    # Points
    pts_node = piece.find('.//{*}Points') or piece.find('.//Points')
    pts = None
    if pts_node is not None:
        da = pts_node.find('.//{*}DataArray') or pts_node.find('.//DataArray')
        if da is not None:
            pts = decode_array(da).flatten().astype(np.float64).reshape(-1, 3)

    if pts is None:
        return None, None

    # PointData fields
    fields = {}
    pd_node = piece.find('.//{*}PointData') or piece.find('.//PointData')
    if pd_node is not None:
        das = pd_node.findall('.//{*}DataArray') or pd_node.findall('.//DataArray')
        for da in das:
            name = da.get('Name')
            if name and name not in fields:
                fields[name] = decode_array(da)

    return pts, fields


# ---------------------------------------------------------------------------
# Trajectory building
# ---------------------------------------------------------------------------

def build_trajectories(series_entries, t_min, t_max):
    """
    从所有时间步读取 VTP，按 (origProcId, origId) 构建粒子轨迹字典。

    Returns
    -------
    dict pid -> {
        'times': ndarray, 'x': ndarray, 'y': ndarray, 'z': ndarray,
        'd': ndarray, 'rho': ndarray, 'nParticle': ndarray, 'T': ndarray
    }
    """
    filtered = sorted((t, p) for t, p in series_entries if t_min <= t <= t_max)
    if not filtered:
        raise ValueError(f"No time steps found in [{t_min}, {t_max}]")

    print(f"\nBuilding trajectories: {len(filtered)} time steps [{t_min:.4f}s – {t_max:.4f}s]")

    # Accumulate as lists for efficiency
    raw = defaultdict(lambda: {k: [] for k in ('times', 'x', 'y', 'z', 'd', 'rho', 'nParticle', 'T')})

    for i, (t, filepath) in enumerate(filtered):
        pts, fields = parse_vtp(filepath)
        if pts is None:
            continue

        N = len(pts)
        orig_id   = fields.get('origId',     np.zeros(N, dtype=int)).astype(int).ravel()
        orig_proc = fields.get('origProcId', np.zeros(N, dtype=int)).astype(int).ravel()
        d_f   = fields.get('d',         np.ones(N) * np.nan).ravel()
        rho_f = fields.get('rho',       np.ones(N) * np.nan).ravel()
        nP_f  = fields.get('nParticle', np.ones(N)).ravel()
        T_f   = fields.get('T',         np.full(N, np.nan)).ravel()

        for j in range(N):
            pid = (int(orig_proc[j]), int(orig_id[j]))
            tr = raw[pid]
            tr['times'].append(t)
            tr['x'].append(float(pts[j, 0]))
            tr['y'].append(float(pts[j, 1]))
            tr['z'].append(float(pts[j, 2]))
            tr['d'].append(float(d_f[j]))
            tr['rho'].append(float(rho_f[j]))
            tr['nParticle'].append(float(nP_f[j]))
            tr['T'].append(float(T_f[j]))

        if (i + 1) % 20 == 0 or i == len(filtered) - 1:
            print(f"  [{i+1:4d}/{len(filtered)}] t={t:.5f}s  {N:4d} particles  "
                  f"(unique tracked: {len(raw)})")

    # Convert to numpy arrays and sort by time
    trajectories = {}
    for pid, tr in raw.items():
        order = np.argsort(tr['times'])
        trajectories[pid] = {k: np.array(tr[k])[order] for k in tr}

    print(f"Total unique particles tracked: {len(trajectories)}")
    return trajectories


# ---------------------------------------------------------------------------
# Plane crossing detection (linear interpolation)
# ---------------------------------------------------------------------------

def find_plane_crossings(trajectories, x_target, direction='positive'):
    """
    对每个粒子轨迹，找出 x = x_target 的精确穿越事件（线性插值）。

    Parameters
    ----------
    direction : 'positive'  – 仅统计从左往右穿越（x 增大方向）
                'negative'  – 仅统计从右往左穿越
                'both'      – 两个方向都统计

    Returns
    -------
    list of dicts, each with keys:
        pid, t, y (mm), z (mm), d, rho, nParticle, T, mass, direction
    """
    crossings = []

    for pid, tr in trajectories.items():
        xs = tr['x']
        ts = tr['times']
        if len(ts) < 2:
            continue

        for i in range(len(ts) - 1):
            x0, x1 = xs[i], xs[i + 1]
            fwd = (x0 < x_target) and (x1 >= x_target)
            bwd = (x0 >= x_target) and (x1 < x_target)

            if direction == 'positive' and not fwd:
                continue
            if direction == 'negative' and not bwd:
                continue
            if direction == 'both' and not (fwd or bwd):
                continue

            dx = x1 - x0
            if abs(dx) < 1e-15:
                continue

            # Linear interpolation factor α ∈ [0, 1]
            alpha = (x_target - x0) / dx

            def lerp(arr):
                return float(arr[i]) + alpha * (float(arr[i + 1]) - float(arr[i]))

            t_c   = lerp(ts)
            y_c   = lerp(tr['y'])
            z_c   = lerp(tr['z'])
            d_c   = lerp(tr['d'])
            rho_c = lerp(tr['rho'])
            nP_c  = lerp(tr['nParticle'])
            T_c   = lerp(tr['T'])

            vol  = (4.0 / 3.0) * np.pi * (d_c / 2.0) ** 3
            mass = vol * rho_c * nP_c

            crossings.append({
                'pid':       pid,
                't':         t_c,
                'y':         y_c * 1000,   # → mm
                'z':         z_c * 1000,   # → mm
                'd':         d_c,
                'rho':       rho_c,
                'nParticle': nP_c,
                'T':         T_c,
                'mass':      mass,
                'direction': 'forward' if fwd else 'backward',
            })

    return crossings


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _lim(arr, pad=1.1):
    return max(np.abs(arr).max(), 0.1) * pad


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    case_dir = os.getcwd()
    vtk_dir  = os.path.join(case_dir, 'VTK', 'lagrangian', 'kinematicCloud')

    # --- optional foamToVTK ---
    run_vtk = input("Run foamToVTK? [1=yes / 0=skip]: ").strip()
    if run_vtk == '1':
        print("Running foamToVTK ...")
        res = subprocess.run(
            'foamToVTK 2>&1 | tail -30',
            shell=True, cwd=case_dir, capture_output=True, text=True
        )
        print(res.stdout)
        if res.returncode != 0:
            print("foamToVTK stderr:", res.stderr)

    # --- build complete file list from disk + series time mapping ---
    series_file = os.path.join(vtk_dir, 'kinematicCloud.vtp.series')
    with open(series_file) as f:
        series = json.load(f)
    # series lookup: basename -> time
    series_time = {e['name']: e['time'] for e in series['files']}

    # scan all .vtp files actually present on disk
    import glob as _glob
    disk_files = sorted(_glob.glob(os.path.join(vtk_dir, 'kinematicCloud_*.vtp')))
    all_entries = []
    missing_time = []
    for fpath in disk_files:
        bname = os.path.basename(fpath)
        if bname in series_time:
            all_entries.append((series_time[bname], fpath))
        else:
            missing_time.append(bname)

    all_entries.sort(key=lambda x: x[0])

    print("\n" + "=" * 70)
    print("Particle Trajectory Plane-Crossing Analysis")
    print("=" * 70)
    print(f"VTK files on disk : {len(disk_files)}")
    print(f"Matched with times: {len(all_entries)}")
    if missing_time:
        print(f"WARNING: {len(missing_time)} file(s) on disk have no time in series file "
              f"and will be skipped: {missing_time[:5]}{'...' if len(missing_time)>5 else ''}")
    print(f"Available time range: {all_entries[0][0]:.5f}s – {all_entries[-1][0]:.5f}s")

    t_min        = float(input("Start time (s) [e.g. 0.0]:    ").strip() or "0.0")
    t_max        = float(input("End   time (s) [e.g. 0.06]:   ").strip() or "0.06")
    x_target_mm  = float(input("Observation plane x (mm) [e.g. 37.5 correspond to 11.5mm in EXP ]: ").strip() or "31.5")
    x_target     = x_target_mm / 1000.0

    # --- build trajectories ---
    trajectories = build_trajectories(all_entries, t_min, t_max)

    # --- find crossings ---
    crossings = find_plane_crossings(trajectories, x_target, direction='positive')

    # unique parcels that crossed
    crossing_pids_set = {c['pid'] for c in crossings}
    n_total_tracked   = len(trajectories)
    n_crossed         = len(crossing_pids_set)
    n_not_crossed     = n_total_tracked - n_crossed

    print(f"\n=== Trajectory Reconstruction Summary ===")
    print(f"Unique parcels tracked (reconstructed) : {n_total_tracked}")
    print(f"  ├─ crossed x = {x_target_mm:.1f} mm                : {n_crossed}"
          f"  ({100*n_crossed/n_total_tracked:.1f}%)")
    print(f"  └─ did NOT cross (still flying / escaped elsewhere): {n_not_crossed}"
          f"  ({100*n_not_crossed/n_total_tracked:.1f}%)")

    print(f"\n=== Crossing Events at x = {x_target_mm:.1f} mm ===")
    print(f"Total crossing events : {len(crossings)}"
          f"  (multi-cross by same parcel: {len(crossings) - n_crossed})")

    if not crossings:
        print("No crossing events found. Check x_target or time range.")
        return

    y_arr = np.array([c['y']    for c in crossings])
    z_arr = np.array([c['z']    for c in crossings])
    m_arr = np.array([c['mass'] for c in crossings])
    t_arr = np.array([c['t']    for c in crossings])
    T_arr = np.array([c['T']    for c in crossings])
    d_arr = np.array([c['d']    for c in crossings])

    # physical particle count (nParticle weighted, each parcel once)
    total_physical = sum(
        next(c['nParticle'] for c in crossings if c['pid'] == pid)
        for pid in crossing_pids_set
    )
    print(f"Physical particles (nParticle sum)     : {total_physical:.0f}")
    print(f"Total crossing mass   : {np.sum(m_arr)*1e9:.4f} ng")
    print(f"Y range (mm)          : {y_arr.min():.3f} – {y_arr.max():.3f}")
    print(f"Z range (mm)          : {z_arr.min():.3f} – {z_arr.max():.3f}")
    print(f"Diameter range (μm)   : {d_arr.min()*1e6:.2f} – {d_arr.max()*1e6:.2f}")
    valid_T = ~np.isnan(T_arr)
    if np.any(valid_T):
        T_mean = np.average(T_arr[valid_T], weights=m_arr[valid_T])
        print(f"T range (K)           : {np.nanmin(T_arr):.1f} – {np.nanmax(T_arr):.1f}")
        print(f"Mass-weighted mean T  : {T_mean:.1f} K")

    xy_lim = _lim(np.concatenate([y_arr, z_arr]))

    # =========================================================
    # Figure 1 – YOZ plane crossing distribution
    # =========================================================
    fig1, axes1 = plt.subplots(2, 2, figsize=(14, 12))
    fig1.suptitle(
        f'Particle Plane-Crossing Distribution at x = {x_target_mm:.1f} mm\n'
        f'(Trajectory Interpolation, t = [{t_min:.3f}s – {t_max:.3f}s], '
        f'{len(crossings)} events)',
        fontsize=12
    )

    # 1.1 scatter, colour = crossing time, size ∝ mass
    sc1 = axes1[0, 0].scatter(
        y_arr, z_arr,
        c=t_arr * 1000,
        s=m_arr / m_arr.max() * 80 + 2,
        cmap='viridis', alpha=0.75, linewidths=0
    )
    plt.colorbar(sc1, ax=axes1[0, 0], label='Crossing time (ms)')
    axes1[0, 0].set(xlabel='Y (mm)', ylabel='Z (mm)',
                    title='Crossing positions (colour = time, size ∝ mass)',
                    xlim=(-xy_lim, xy_lim), ylim=(-xy_lim, xy_lim), aspect='equal')
    axes1[0, 0].grid(alpha=0.3)

    # 1.2 2-D mass density heatmap
    bins2d = np.linspace(-xy_lim, xy_lim, 41)
    H, ye, ze = np.histogram2d(y_arr, z_arr, bins=[bins2d, bins2d], weights=m_arr)
    H_pct = H / np.sum(m_arr) * 100
    im12 = axes1[0, 1].imshow(
        H_pct.T, origin='lower', aspect='auto',
        extent=[ye[0], ye[-1], ze[0], ze[-1]], cmap='YlOrRd'
    )
    plt.colorbar(im12, ax=axes1[0, 1], label='Mass fraction per bin (%)')
    axes1[0, 1].set(xlabel='Y (mm)', ylabel='Z (mm)',
                    title='Mass-weighted crossing density [% of total]',
                    xlim=(-xy_lim, xy_lim), ylim=(-xy_lim, xy_lim))

    # 1.3 Y marginal
    yh, ye2 = np.histogram(y_arr, bins=50, weights=m_arr)
    yc = (ye2[:-1] + ye2[1:]) / 2
    axes1[1, 0].bar(yc, yh * 1e9, width=ye2[1] - ye2[0], alpha=0.7, color='steelblue')
    axes1[1, 0].set(xlabel='Y (mm)', ylabel='Crossing mass (ng)',
                    title='Mass distribution along Y')
    axes1[1, 0].grid(alpha=0.3)

    # 1.4 Z marginal
    zh, ze2 = np.histogram(z_arr, bins=50, weights=m_arr)
    zc = (ze2[:-1] + ze2[1:]) / 2
    axes1[1, 1].bar(zc, zh * 1e9, width=ze2[1] - ze2[0], alpha=0.7, color='seagreen')
    axes1[1, 1].set(xlabel='Z (mm)', ylabel='Crossing mass (ng)',
                    title='Mass distribution along Z')
    axes1[1, 1].grid(alpha=0.3)

    plt.tight_layout()
    tag = f'x{x_target_mm:.1f}mm_t{t_min*1000:.1f}-{t_max*1000:.1f}ms'
    out1 = os.path.join(case_dir, f'crossing_yoz_{tag}.png')
    fig1.savefig(out1, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {out1}")
    plt.close(fig1)

    # =========================================================
    # Figure 2 – radial distribution + mass rate vs time
    # =========================================================
    fig2, axes2 = plt.subplots(1, 3, figsize=(18, 5))

    y_ctr = np.average(y_arr, weights=m_arr)
    z_ctr = np.average(z_arr, weights=m_arr)
    r = np.sqrt((y_arr - y_ctr) ** 2 + (z_arr - z_ctr) ** 2)

    r_bins = np.linspace(0, r.max() * 1.02, 40)
    rh, re = np.histogram(r, bins=r_bins, weights=m_arr)
    rc = (re[:-1] + re[1:]) / 2
    ring_areas = np.pi * (re[1:] ** 2 - re[:-1] ** 2)
    r_density = np.where(ring_areas > 0, rh / ring_areas * 1e9, 0)

    axes2[0].bar(rc, r_density, width=re[1] - re[0], alpha=0.7, color='purple')
    axes2[0].set(xlabel='Radial distance (mm)', ylabel='Mass density (ng/mm²)',
                 title=f'Radial mass density\n'
                       f'(centroid Y={y_ctr:.3f}, Z={z_ctr:.3f} mm)')
    axes2[0].grid(alpha=0.3)

    cum_r = np.cumsum(rh)
    axes2[1].plot(rc, cum_r / cum_r[-1] * 100, 'b-', lw=2)
    r50 = rc[np.argmax(cum_r >= 0.50 * cum_r[-1])]
    r90 = rc[np.argmax(cum_r >= 0.90 * cum_r[-1])]
    axes2[1].axhline(50, color='r', ls='--', alpha=0.7, label=f'50%  r={r50:.3f}mm')
    axes2[1].axhline(90, color='g', ls='--', alpha=0.7, label=f'90%  r={r90:.3f}mm')
    axes2[1].axvline(r50, color='r', ls=':', alpha=0.5)
    axes2[1].axvline(r90, color='g', ls=':', alpha=0.5)
    axes2[1].set(xlabel='Radial distance (mm)', ylabel='Cumulative mass (%)',
                 title=f'Cumulative radial mass\n(r50={r50:.3f}mm, r90={r90:.3f}mm)')
    axes2[1].legend(fontsize=8)
    axes2[1].grid(alpha=0.3)
    print(f"\nMass-weighted centroid : Y={y_ctr:.3f} mm, Z={z_ctr:.3f} mm")
    print(f"50% mass radius        : {r50:.3f} mm")
    print(f"90% mass radius        : {r90:.3f} mm")

    # Crossing mass rate vs time
    t_bins = np.linspace(t_arr.min(), t_arr.max(), min(40, len(crossings) + 1))
    if len(t_bins) < 3:
        t_bins = np.linspace(t_min, t_max, 40)
    mh_t, te_t = np.histogram(t_arr, bins=t_bins, weights=m_arr)
    nh_t, _    = np.histogram(t_arr, bins=t_bins)
    tc_t = (te_t[:-1] + te_t[1:]) / 2
    dt_t = te_t[1] - te_t[0]

    ax2r = axes2[2].twinx()
    axes2[2].bar(tc_t * 1000, mh_t / dt_t * 1e9,
                 width=dt_t * 1000, alpha=0.5, color='tomato', label='Mass rate')
    ax2r.plot(tc_t * 1000, nh_t, 'b-o', ms=3, label='Count')
    axes2[2].set(xlabel='Time (ms)', ylabel='Mass flow rate  (ng/s)')
    axes2[2].yaxis.label.set_color('tomato')
    ax2r.set_ylabel('Crossing count per bin', color='b')
    axes2[2].set_title('Crossing mass rate & count vs time')
    axes2[2].grid(alpha=0.3)
    from matplotlib.lines import Line2D
    axes2[2].legend(
        [Line2D([0], [0], color='tomato', lw=4, alpha=0.5),
         Line2D([0], [0], color='b', marker='o')],
        ['Mass rate (ng/s)', 'Crossing count'],
        loc='upper left', fontsize=8
    )

    plt.tight_layout()
    out2 = os.path.join(case_dir, f'crossing_radial_{tag}.png')
    fig2.savefig(out2, dpi=150, bbox_inches='tight')
    print(f"Saved: {out2}")
    plt.close(fig2)

    # =========================================================
    # Figure 3 – Particle trajectories
    # =========================================================
    crossing_pids = list({c['pid'] for c in crossings})
    max_plot = 300
    plot_pids = crossing_pids[:max_plot]

    fig3, axes3 = plt.subplots(1, 3, figsize=(18, 6))
    fig3.suptitle(
        f'Particle trajectories crossing x = {x_target_mm:.1f} mm\n'
        f'({len(plot_pids)} of {len(crossing_pids)} crossing particles shown)',
        fontsize=11
    )

    # Colour each trajectory by its crossing time for visual clarity
    c_times = []
    for pid in plot_pids:
        # take the first crossing time for this pid
        ct = next((c['t'] for c in crossings if c['pid'] == pid), 0.0)
        c_times.append(ct)
    c_times = np.array(c_times)
    norm_c  = plt.Normalize(c_times.min(), c_times.max())
    cmap_c  = plt.cm.viridis

    # 3.1 XY projection
    for k, pid in enumerate(plot_pids):
        tr = trajectories[pid]
        axes3[0].plot(
            np.array(tr['x']) * 1000, np.array(tr['y']) * 1000,
            lw=0.6, alpha=0.4, color=cmap_c(norm_c(c_times[k]))
        )
    axes3[0].axvline(x_target_mm, color='red', lw=1.5, ls='--', label=f'x={x_target_mm}mm')
    axes3[0].set(xlabel='X (mm)', ylabel='Y (mm)', title='XY trajectory projection')
    axes3[0].legend(fontsize=8)
    axes3[0].grid(alpha=0.3)

    # 3.2 XZ projection
    for k, pid in enumerate(plot_pids):
        tr = trajectories[pid]
        axes3[1].plot(
            np.array(tr['x']) * 1000, np.array(tr['z']) * 1000,
            lw=0.6, alpha=0.4, color=cmap_c(norm_c(c_times[k]))
        )
    axes3[1].axvline(x_target_mm, color='red', lw=1.5, ls='--', label=f'x={x_target_mm}mm')
    axes3[1].set(xlabel='X (mm)', ylabel='Z (mm)', title='XZ trajectory projection')
    axes3[1].legend(fontsize=8)
    axes3[1].grid(alpha=0.3)

    # 3.3 YZ projection + interpolated crossing points
    for k, pid in enumerate(plot_pids):
        tr = trajectories[pid]
        axes3[2].plot(
            np.array(tr['y']) * 1000, np.array(tr['z']) * 1000,
            lw=0.5, alpha=0.25, color='gray'
        )
    sc3 = axes3[2].scatter(
        y_arr, z_arr,
        c=m_arr / m_arr.max(),
        s=20, cmap='hot_r', alpha=0.85,
        zorder=5, linewidths=0
    )
    plt.colorbar(sc3, ax=axes3[2], label='Relative mass')
    axes3[2].set(xlabel='Y (mm)', ylabel='Z (mm)',
                 title=f'YZ trajectories + crossing points at x={x_target_mm:.1f}mm',
                 xlim=(-xy_lim, xy_lim), ylim=(-xy_lim, xy_lim), aspect='equal')
    axes3[2].grid(alpha=0.3)

    sm3 = plt.cm.ScalarMappable(norm=norm_c, cmap=cmap_c)
    sm3.set_array([])
    cb3 = plt.colorbar(sm3, ax=axes3[0])
    cb3.set_label('Crossing time (s)')
    cb3b = plt.colorbar(sm3, ax=axes3[1])
    cb3b.set_label('Crossing time (s)')

    plt.tight_layout()
    out3 = os.path.join(case_dir, f'particle_trajectories_{tag}.png')
    fig3.savefig(out3, dpi=150, bbox_inches='tight')
    print(f"Saved: {out3}")
    plt.close(fig3)

    # =========================================================
    # Figure 4 – Temperature distribution at crossing (optional)
    # =========================================================
    if np.any(valid_T):
        fig4, axes4 = plt.subplots(1, 2, figsize=(14, 5))
        fig4.suptitle(f'Particle temperature at plane crossing (x={x_target_mm:.1f}mm)')

        sc4 = axes4[0].scatter(
            y_arr[valid_T], z_arr[valid_T],
            c=T_arr[valid_T],
            s=m_arr[valid_T] / m_arr.max() * 80 + 2,
            cmap='hot', alpha=0.75, linewidths=0
        )
        plt.colorbar(sc4, ax=axes4[0], label='Temperature (K)')
        axes4[0].set(xlabel='Y (mm)', ylabel='Z (mm)',
                     title='Temperature at crossing (size ∝ mass)',
                     xlim=(-xy_lim, xy_lim), ylim=(-xy_lim, xy_lim), aspect='equal')
        axes4[0].grid(alpha=0.3)

        Th, Te = np.histogram(T_arr[valid_T], bins=50, weights=m_arr[valid_T])
        Tc = (Te[:-1] + Te[1:]) / 2
        axes4[1].bar(Tc, Th * 1e9, width=Te[1] - Te[0], alpha=0.7, color='orangered')
        axes4[1].set(xlabel='Temperature (K)', ylabel='Mass (ng)',
                     title='Mass-weighted temperature distribution at crossing')
        axes4[1].grid(alpha=0.3)
        T_mean_marked = np.average(T_arr[valid_T], weights=m_arr[valid_T])
        axes4[1].axvline(T_mean_marked, color='k', ls='--', lw=1.5,
                         label=f'Mean = {T_mean_marked:.0f} K')
        axes4[1].legend()

        plt.tight_layout()
        out4 = os.path.join(case_dir, f'crossing_temperature_{tag}.png')
        fig4.savefig(out4, dpi=150, bbox_inches='tight')
        print(f"Saved: {out4}")
        plt.close(fig4)

    print("\n✓ Analysis complete.")


if __name__ == '__main__':
    main()
