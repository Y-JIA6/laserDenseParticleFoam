#!/usr/bin/env python3
"""
直接从 OpenFOAM processorN 目录读取拉格朗日粒子数据，
无需 foamToVTK / VTK 库。

算法与 analyze_crossing.py 相同：
1. 扫描所有 processorN/<time>/lagrangian/kinematicCloud/ 目录，
   收集时间步列表
2. 逐时间步从所有 proc 合并粒子数据，以 (origProcId, origId) 为键
   构建完整轨迹
3. 线性插值检测粒子穿越观察面 x = x_target 的精确事件
4. 统计并输出 4 张图表（与原脚本相同）
"""

import re
import os
import glob

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict


# ---------------------------------------------------------------------------
# OpenFOAM lagrangian file parsers
# ---------------------------------------------------------------------------

# Matches both '// * * * * //'
# and '// ****...* //'
_HEADER_SEP = re.compile(r'// [\* ]+//')


def _body(filepath):
    """Return the data section of an OpenFOAM field file.

    OpenFOAM files have the structure:
        /* header */
        FoamFile { ... }
        // * * * ... * //      <- first separator
        <data>
        // ****...*** //      <- second separator
    We return the text between the first and last separator.
    """
    with open(filepath, 'r') as fh:
        content = fh.read()
    parts = _HEADER_SEP.split(content)
    # parts: [header, data, trailing_newline_or_empty]
    if len(parts) >= 2:
        return parts[1].strip()
    return content.strip()


def parse_of_positions(filepath):
    """
    解析 OpenFOAM lagrangian positions 文件。
    格式:
        N
        (
        (x y z) cellId
        ...
        )
    返回 ndarray shape (N, 3)，若 N==0 返回空数组。
    """
    body = _body(filepath)
    lines = [l.strip() for l in body.splitlines() if l.strip()]
    if not lines:
        return np.empty((0, 3), dtype=np.float64)

    try:
        n = int(lines[0])
    except ValueError:
        return np.empty((0, 3), dtype=np.float64)

    if n == 0:
        return np.empty((0, 3), dtype=np.float64)

    pts = []
    vec_re = re.compile(r'\(([^)]+)\)')
    for line in lines[1:]:
        m = vec_re.match(line)
        if m:
            vals = [float(v) for v in m.group(1).split()]
            if len(vals) == 3:
                pts.append(vals)
    return np.array(pts, dtype=np.float64)


def parse_of_scalar(filepath):
    """
    解析 OpenFOAM scalarField / labelField 文件。
    支持三种格式：
        N{value}             —— 均匀（所有粒子相同值）
        N(v1 v2 ... vN)      —— 非均匀，内联（单行）
        N                    —— 非均匀，多行：值在 ( ... ) 块里每行一个
        (
        v1
        v2
        ...
        )
    返回 ndarray shape (N,)，失败返回 None。
    """
    body = _body(filepath)
    # 去掉内嵌注释行
    body = re.sub(r'//[^\n]*', '', body).strip()

    # 均匀格式: N{value}
    m = re.match(r'^(\d+)\s*\{([^}]+)\}', body)
    if m:
        n = int(m.group(1))
        val = float(m.group(2).strip())
        return np.full(n, val, dtype=np.float64)

    # 非均匀内联格式: N(v1 v2 ...)
    m = re.match(r'^(\d+)\s*\((.+)\)\s*$', body, re.DOTALL)
    if m:
        inner = m.group(2).strip()
        return np.array(inner.split(), dtype=np.float64)

    return None


def parse_of_vector(filepath):
    """
    解析 OpenFOAM vectorField 文件。
    支持均匀格式: N{(vx vy vz)}
    和多行格式:   N\n(\n(v1x v1y v1z)\n...\n)
    返回 ndarray shape (N, 3)，失败返回 None。
    """
    body = _body(filepath)
    body = re.sub(r'//[^\n]*', '', body).strip()

    # 均匀格式: N{(vx vy vz)}
    m = re.match(r'^(\d+)\s*\{\s*\(([^)]+)\)\s*\}', body)
    if m:
        n = int(m.group(1))
        vals = [float(v) for v in m.group(2).split()]
        return np.tile(vals, (n, 1)).astype(np.float64)

    # 多行格式
    n_match = re.match(r'^(\d+)', body)
    if not n_match:
        return None
    n = int(n_match.group(1))
    if n == 0:
        return np.empty((0, 3), dtype=np.float64)
    vec_re = re.compile(r'\(([^()]+)\)')
    result = []
    for v in vec_re.findall(body):
        parts = v.split()
        if len(parts) == 3:
            result.append([float(x) for x in parts])
    if len(result) == n:
        return np.array(result, dtype=np.float64)
    return None


# ---------------------------------------------------------------------------
# Read one time step from all processors
# ---------------------------------------------------------------------------

def read_timestep_of(case_dir, time_str, n_procs, cloud='kinematicCloud'):
    """
    从所有 processorN 目录读取指定时间步的粒子数据并合并。

    Returns
    -------
    pts    : ndarray (N_total, 3) or None
    fields : dict of ndarray (N_total,)  or None
    """
    pts_list   = []
    fld_lists  = defaultdict(list)
    FIELDS     = ('origId', 'origProcId', 'd', 'rho', 'nParticle', 'T')

    for proc in range(n_procs):
        lag_dir  = os.path.join(case_dir, f'processor{proc}',
                                time_str, 'lagrangian', cloud)
        pos_file = os.path.join(lag_dir, 'positions')
        if not os.path.exists(pos_file):
            continue

        pts = parse_of_positions(pos_file)
        if pts is None or len(pts) == 0:
            continue

        N = len(pts)
        pts_list.append(pts)

        for field in FIELDS:
            fpath = os.path.join(lag_dir, field)
            if os.path.exists(fpath):
                arr = parse_of_scalar(fpath)
                fld_lists[field].append(
                    arr if (arr is not None and len(arr) == N)
                    else np.full(N, np.nan)
                )
            else:
                fld_lists[field].append(np.full(N, np.nan))

        # 矢量场 U → Ux, Uy, Uz
        u_fpath = os.path.join(lag_dir, 'U')
        if os.path.exists(u_fpath):
            U_arr = parse_of_vector(u_fpath)
            if U_arr is not None and len(U_arr) == N:
                fld_lists['Ux'].append(U_arr[:, 0])
                fld_lists['Uy'].append(U_arr[:, 1])
                fld_lists['Uz'].append(U_arr[:, 2])
            else:
                for k in ('Ux', 'Uy', 'Uz'):
                    fld_lists[k].append(np.full(N, np.nan))
        else:
            for k in ('Ux', 'Uy', 'Uz'):
                fld_lists[k].append(np.full(N, np.nan))

    if not pts_list:
        return None, None

    all_pts = np.vstack(pts_list)
    all_fields = {k: np.concatenate(v) for k, v in fld_lists.items()}
    return all_pts, all_fields


# ---------------------------------------------------------------------------
# Trajectory building
# ---------------------------------------------------------------------------

def build_trajectories(case_dir, time_entries, n_procs, t_min, t_max,
                       cloud='kinematicCloud'):
    """
    从所有时间步、所有 processor 读取粒子数据，
    以 (origProcId, origId) 为键构建完整轨迹字典。

    Parameters
    ----------
    time_entries : list of (float, str)   — (时间值, 时间目录名)

    Returns
    -------
    dict pid -> {
        'times': ndarray, 'x': ndarray, 'y': ndarray, 'z': ndarray,
        'd': ndarray, 'rho': ndarray, 'nParticle': ndarray, 'T': ndarray
    }
    """
    filtered = sorted((t, s) for t, s in time_entries if t_min <= t <= t_max)
    if not filtered:
        raise ValueError(f"在 [{t_min}, {t_max}] 范围内没有找到时间步")

    print(f"\n构建轨迹：{len(filtered)} 个时间步  [{t_min:.4f}s – {t_max:.4f}s]")

    raw = defaultdict(lambda: {k: [] for k in
                               ('times', 'x', 'y', 'z', 'd', 'rho', 'nParticle', 'T', 'Ux', 'Uy', 'Uz')})

    for i, (t, tstr) in enumerate(filtered):
        pts, fields = read_timestep_of(case_dir, tstr, n_procs, cloud)
        if pts is None:
            continue

        N = len(pts)
        orig_id   = fields.get('origId',   np.zeros(N)).astype(int).ravel()
        orig_proc = fields.get('origProcId', np.zeros(N)).astype(int).ravel()
        d_f   = fields.get('d',         np.full(N, np.nan)).ravel()
        rho_f = fields.get('rho',       np.full(N, np.nan)).ravel()
        nP_f  = fields.get('nParticle', np.ones(N)).ravel()
        T_f   = fields.get('T',         np.full(N, np.nan)).ravel()
        Ux_f  = fields.get('Ux',        np.full(N, np.nan)).ravel()
        Uy_f  = fields.get('Uy',        np.full(N, np.nan)).ravel()
        Uz_f  = fields.get('Uz',        np.full(N, np.nan)).ravel()

        for j in range(N):
            pid = (int(orig_proc[j]), int(orig_id[j]))
            tr  = raw[pid]
            tr['times'].append(t)
            tr['x'].append(float(pts[j, 0]))
            tr['y'].append(float(pts[j, 1]))
            tr['z'].append(float(pts[j, 2]))
            tr['d'].append(float(d_f[j]))
            tr['rho'].append(float(rho_f[j]))
            tr['nParticle'].append(float(nP_f[j]))
            tr['T'].append(float(T_f[j]))
            tr['Ux'].append(float(Ux_f[j]))
            tr['Uy'].append(float(Uy_f[j]))
            tr['Uz'].append(float(Uz_f[j]))

        if (i + 1) % 20 == 0 or i == len(filtered) - 1:
            print(f"  [{i+1:4d}/{len(filtered)}] t={t:.5f}s  {N:4d} 粒子  "
                  f"(已追踪唯一粒子: {len(raw)})")

    trajectories = {}
    for pid, tr in raw.items():
        order = np.argsort(tr['times'])
        trajectories[pid] = {k: np.array(tr[k])[order] for k in tr}

    print(f"追踪到的唯一粒子总数: {len(trajectories)}")
    return trajectories


# ---------------------------------------------------------------------------
# Plane crossing detection (linear interpolation)
# ---------------------------------------------------------------------------

def find_plane_crossings(trajectories, x_target, direction='positive'):
    """
    对每个粒子轨迹找出 x = x_target 的精确穿越事件（三次 Lagrange 插值）。

    穿越发生在相邻时间步 [t_i, t_{i+1}] 之间，选取 4 点模板
    (i-1, i, i+1, i+2)，即穿越区间两侧各 2 个已保存时间步，
    拟合三次多项式 x(tau)，用 np.roots 解出精确穿越时刻 tau_c，
    再用同一 Lagrange 权重插值所有场量。
    时间步不足 4 个时退化为二次（3 点）；不足 3 个退化为线性。

    direction : 'positive' 左→右 | 'negative' 右→左 | 'both' 两方向
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

            n_pts = len(ts)

            # --- Stencil selection: cubic (4-pt) preferred; quadratic (3-pt)
            #     when n_pts < 4; linear when n_pts < 3. ---
            # Crossing is between i and i+1; stencil spans i-1..i+2 so that
            # there are 2 saved steps on each side of the crossing interval.
            interp_deg = 0
            jj = None

            if n_pts >= 4:
                j0 = max(0, min(i - 1, n_pts - 4))   # clamp to valid range
                jj = tuple(range(j0, j0 + 4))
                interp_deg = 3
            elif n_pts == 3:
                jj = (0, 1, 2)
                interp_deg = 2

            poly_ok = False
            if jj is not None:
                t_sten = np.array([ts[k] for k in jj])
                dt_ref = t_sten[-1] - t_sten[0]
                poly_ok = dt_ref > 1e-15

            if poly_ok:
                t0r     = t_sten[0]
                tau_arr = (t_sten - t0r) / dt_ref   # normalised to [0, 1]
                tau_i   = (ts[i]     - t0r) / dt_ref
                tau_i1  = (ts[i + 1] - t0r) / dt_ref

                # Fit polynomial x(tau) and solve x(tau_c) = x_target
                x_sten       = np.array([xs[k] for k in jj])
                cx           = np.polyfit(tau_arr, x_sten, interp_deg)
                cx_shifted   = cx.copy()
                cx_shifted[-1] -= x_target          # p(tau) - x_target = 0
                roots = np.roots(cx_shifted)

                lo, hi = min(tau_i, tau_i1), max(tau_i, tau_i1)
                real_roots = [r.real for r in roots
                              if abs(r.imag) < 1e-9
                              and lo - 1e-9 <= r.real <= hi + 1e-9]

                if real_roots:
                    if len(real_roots) == 1:
                        tau_c = real_roots[0]
                    else:
                        tau_lin = tau_i + (x_target - x0) / dx * (tau_i1 - tau_i)
                        tau_c   = min(real_roots, key=lambda r: abs(r - tau_lin))
                else:
                    # No root in segment: fall back to linear-in-tau
                    tau_c = tau_i + (x_target - x0) / dx * (tau_i1 - tau_i)

                t_c = tau_c * dt_ref + t0r

                # Lagrange basis weights at tau_c (generic for 3- or 4-point stencil)
                n_j = len(jj)
                L   = np.ones(n_j)
                for k in range(n_j):
                    for m in range(n_j):
                        if m != k:
                            L[k] *= (tau_c - tau_arr[m]) / (tau_arr[k] - tau_arr[m])

                def qeval(arr, _j=jj, _L=L):
                    return sum(float(arr[_j[k]]) * _L[k] for k in range(len(_j)))

                y_c   = qeval(tr['y'])
                z_c   = qeval(tr['z'])
                d_c   = qeval(tr['d'])
                rho_c = qeval(tr['rho'])
                nP_c  = qeval(tr['nParticle'])
                T_c   = max(qeval(tr['T']), 300.0)   # clamp: cubic interp undershoot
                Ux_c  = qeval(tr['Ux'])
                Uy_c  = qeval(tr['Uy'])
                Uz_c  = qeval(tr['Uz'])

            else:
                # Linear fallback (trajectory has fewer than 3 time steps)
                alpha = (x_target - x0) / dx

                def lerp(arr, _i=i, _a=alpha):
                    return float(arr[_i]) + _a * (float(arr[_i + 1]) - float(arr[_i]))

                t_c   = lerp(ts)
                y_c   = lerp(tr['y'])
                z_c   = lerp(tr['z'])
                d_c   = lerp(tr['d'])
                rho_c = lerp(tr['rho'])
                nP_c  = lerp(tr['nParticle'])
                T_c   = max(lerp(tr['T']), 300.0)   # clamp: interp undershoot
                Ux_c  = lerp(tr['Ux'])
                Uy_c  = lerp(tr['Uy'])
                Uz_c  = lerp(tr['Uz'])

            speed_c = float(np.sqrt(Ux_c**2 + Uy_c**2 + Uz_c**2))

            vol  = (4.0 / 3.0) * np.pi * (d_c / 2.0) ** 3
            mass = vol * rho_c * nP_c

            crossings.append({
                'pid':       pid,
                't':         t_c,
                'y':         y_c * 1000,   # m → mm
                'z':         z_c * 1000,
                'd':         d_c,
                'rho':       rho_c,
                'nParticle': nP_c,
                'T':         T_c,
                'speed':     speed_c,
                'mass':      mass,
                'direction': 'forward' if fwd else 'backward',
            })

    return crossings


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _lim(arr, pad=1.1):
    return max(np.abs(arr).max(), 0.1) * pad


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    case_dir = os.getcwd()

    # --- 扫描 processorN 目录，收集处理器数 & 时间步 ---
    proc_dirs = sorted(glob.glob(os.path.join(case_dir, 'processor[0-9]*')))
    if not proc_dirs:
        raise RuntimeError(f"在 {case_dir} 下未找到 processorN 目录")

    n_procs = len(proc_dirs)
    print(f"找到 {n_procs} 个处理器目录: processor0 … processor{n_procs-1}")

    # 从所有 proc 中收集存在粒子数据的时间步
    time_set = set()
    cloud = 'kinematicCloud'
    for pdir in proc_dirs:
        for tdir in os.listdir(pdir):
            lag = os.path.join(pdir, tdir, 'lagrangian', cloud, 'positions')
            if os.path.exists(lag):
                time_set.add(tdir)

    if not time_set:
        raise RuntimeError("没有找到任何含粒子数据的时间步，请确认 cloud 名称正确。")

    # 将时间目录名转为浮点数排序
    def to_float(s):
        try:
            return float(s)
        except ValueError:
            return None

    time_entries = sorted(
        [(to_float(s), s) for s in time_set if to_float(s) is not None]
    )

    print(f"\n可用时间步数: {len(time_entries)}")
    print(f"时间范围: {time_entries[0][0]:.5f}s – {time_entries[-1][0]:.5f}s")

    t_min       = float(input("起始时间 (s) [如 0.0]:    ").strip() or "0.0")
    t_max       = float(input("结束时间 (s) [如 0.06]:   ").strip() or "0.06")
    x_target_mm = float(input("观察面 x 坐标 (mm) [如 7]: ").strip() or "7.0")
    x_target    = x_target_mm / 1000.0
    r_max_mm    = float(input("径向统计最大半径 (mm, 5=自动): ").strip() or "5") or None
    r_n_bins    = int(input("径向分段数 [默认 35]:          ").strip() or "35")
    smooth_sigma = float(input("平滑 sigma (mm, 0=不平滑) [默认 0.05]: ").strip() or "0.05")

    # --- 构建轨迹 ---
    trajectories = build_trajectories(case_dir, time_entries, n_procs,
                                      t_min, t_max, cloud)

    # --- 检测穿越事件 ---
    crossings = find_plane_crossings(trajectories, x_target, direction='positive')

    crossing_pids_set = {c['pid'] for c in crossings}
    n_total_tracked   = len(trajectories)
    n_crossed         = len(crossing_pids_set)
    n_not_crossed     = n_total_tracked - n_crossed

    print(f"\n=== 轨迹重建统计 ===")
    print(f"追踪到的唯一粒子 (parcel) 总数 : {n_total_tracked}")
    print(f"  ├─ 穿越 x = {x_target_mm:.1f} mm          : {n_crossed}"
          f"  ({100*n_crossed/max(n_total_tracked,1):.1f}%)")
    print(f"  └─ 未穿越（仍在飞行/逃出其他方向）: {n_not_crossed}"
          f"  ({100*n_not_crossed/max(n_total_tracked,1):.1f}%)")

    print(f"\n=== x = {x_target_mm:.1f} mm 处穿越事件 ===")
    print(f"总穿越事件数 : {len(crossings)}"
          f"  (同一粒子多次穿越: {len(crossings) - n_crossed})")

    if not crossings:
        print("没有找到穿越事件，请检查 x_target 或时间范围。")
        return

    y_arr = np.array([c['y']    for c in crossings])
    z_arr = np.array([c['z']    for c in crossings])
    m_arr = np.array([c['mass'] for c in crossings])
    t_arr = np.array([c['t']    for c in crossings])
    T_arr = np.array([c['T']    for c in crossings])
    d_arr = np.array([c['d']    for c in crossings])

    total_physical = sum(
        next(c['nParticle'] for c in crossings if c['pid'] == pid)
        for pid in crossing_pids_set
    )
    print(f"物理粒子数 (nParticle 之和)     : {total_physical:.0f}")
    print(f"穿越总质量                       : {np.sum(m_arr)*1e3:.6f} g")
    print(f"Y 范围 (mm)                      : {y_arr.min():.3f} – {y_arr.max():.3f}")
    print(f"Z 范围 (mm)                      : {z_arr.min():.3f} – {z_arr.max():.3f}")
    print(f"直径范围 (μm)                    : {d_arr.min()*1e6:.2f} – {d_arr.max()*1e6:.2f}")

    valid_T = ~np.isnan(T_arr)
    if np.any(valid_T):
        T_mean = np.average(T_arr[valid_T], weights=m_arr[valid_T])
        print(f"温度范围 (K)                     : {np.nanmin(T_arr):.1f} – {np.nanmax(T_arr):.1f}")
        print(f"质量加权平均温度                  : {T_mean:.1f} K")

    xy_lim = _lim(np.concatenate([y_arr, z_arr]))

    tag = f'x{x_target_mm:.1f}mm_t{t_min*1000:.1f}-{t_max*1000:.1f}ms'

    # Mass-weighted centroid and radial distance (reused in fig1 and fig2)
    y_ctr = np.average(y_arr, weights=m_arr)
    z_ctr = np.average(z_arr, weights=m_arr)
    r_arr = np.sqrt(y_arr ** 2 + z_arr ** 2)   # radial distance from origin
    r_centroid = np.sqrt(y_ctr ** 2 + z_ctr ** 2)   # centroid's radial distance from origin
    r_upper = (r_max_mm if r_max_mm is not None else r_arr.max() * 1.02)

    # =========================================================
    # Figure 1 – YOZ 平面穿越分布
    # =========================================================
    fig1, axes1 = plt.subplots(2, 2, figsize=(14, 12))
    # Replace bottom-left 2-D axes with a 3-D axes
    axes1[1, 0].remove()
    ax_3d = fig1.add_subplot(2, 2, 3, projection='3d')
    fig1.suptitle(
        f'Particle Plane-Crossing Distribution  x = {x_target_mm:.1f} mm\n'
        f'(Trajectory Interpolation, t = [{t_min:.3f}s - {t_max:.3f}s], {len(crossings)} events)',
        fontsize=12
    )

    sc1 = axes1[0, 0].scatter(
        y_arr, z_arr, c=d_arr * 1e6,
        s=m_arr / m_arr.max() * 80 + 2,
        cmap='plasma', alpha=0.75, linewidths=0
    )
    plt.colorbar(sc1, ax=axes1[0, 0], label='Particle diameter (μm)')
    axes1[0, 0].set(xlabel='Y (mm)', ylabel='Z (mm)',
                    title='Crossing positions (colour = diameter)',

                    xlim=(-xy_lim, xy_lim), ylim=(-xy_lim, xy_lim), aspect='equal')
    axes1[0, 0].grid(alpha=0.3)

    bins2d = np.linspace(-xy_lim, xy_lim, 41)
    H, ye, ze = np.histogram2d(y_arr, z_arr, bins=[bins2d, bins2d], weights=m_arr)
    H_pct = H / np.sum(m_arr) * 100
    im12 = axes1[0, 1].imshow(
        H_pct.T, origin='lower', aspect='auto',
        extent=[ye[0], ye[-1], ze[0], ze[-1]], cmap='YlOrRd'
    )
    plt.colorbar(im12, ax=axes1[0, 1], label='Mass fraction (%/bin)')
    axes1[0, 1].set(xlabel='Y (mm)', ylabel='Z (mm)',
                    title='Mass-weighted crossing density [% of total]',

                    xlim=(-xy_lim, xy_lim), ylim=(-xy_lim, xy_lim))

    # --- Bottom-left: 3D surface of mass distribution on the YZ observation plane ---
    n3d = 25
    bins3d = np.linspace(-xy_lim, xy_lim, n3d + 1)
    H3d, ye3d, ze3d = np.histogram2d(y_arr, z_arr, bins=[bins3d, bins3d], weights=m_arr)
    yc3d = (ye3d[:-1] + ye3d[1:]) / 2
    zc3d = (ze3d[:-1] + ze3d[1:]) / 2
    YC3d, ZC3d = np.meshgrid(yc3d, zc3d, indexing='ij')
    surf = ax_3d.plot_surface(YC3d, ZC3d, H3d * 1e3,
                              cmap='hot_r', edgecolor='none', alpha=0.9)
    fig1.colorbar(surf, ax=ax_3d, shrink=0.45, pad=0.12, label='Mass (g)')
    ax_3d.set_xlabel('Y (mm)', labelpad=6)
    ax_3d.set_ylabel('Z (mm)', labelpad=6)
    ax_3d.set_zlabel('Mass (g)', labelpad=6)
    ax_3d.set_title('3D mass distribution on YZ plane')
    ax_3d.view_init(elev=30, azim=-60)

    # --- Bottom-right: radial mass distribution along 8 cuts every 45° ---
    phi_arr = np.degrees(np.arctan2(z_arr, y_arr))  # -180 … +180
    r_bins_cut = np.linspace(0, r_upper, r_n_bins + 1)
    rc_cut = (r_bins_cut[:-1] + r_bins_cut[1:]) / 2
    colors_8 = plt.cm.tab10(np.linspace(0, 0.8, 8))
    for idx, ang in enumerate(range(0, 360, 45)):
        dphi = ((phi_arr - ang + 180) % 360) - 180
        mask = np.abs(dphi) <= 22.5
        if mask.sum() < 2:
            continue
        rh_cut, _ = np.histogram(r_arr[mask], bins=r_bins_cut, weights=m_arr[mask])
        axes1[1, 1].plot(rc_cut, rh_cut * 1e3, color=colors_8[idx],
                         lw=1.5, label=f'{ang}\u00b0')
    axes1[1, 1].set(
        xlabel='Radial distance (mm)',
        ylabel='Crossing mass (g)',
        title=f'Radial mass: 8 cuts \u00d7 45\u00b0\n'
              f'(centroid Y={y_ctr:.2f} mm, Z={z_ctr:.2f} mm)'
    )
    axes1[1, 1].legend(fontsize=7, ncol=2, loc='upper right')
    axes1[1, 1].grid(alpha=0.3)

    plt.tight_layout()
    out1 = os.path.join(case_dir, f'crossing_yoz_{tag}.png')
    fig1.savefig(out1, dpi=150, bbox_inches='tight')
    print(f"\n已保存: {out1}")
    plt.close(fig1)

    # =========================================================
    # Figure 2 – 径向分布 + 平滑密度曲线 + 质量流量随时间
    # =========================================================
    from scipy.ndimage import gaussian_filter1d

    fig2, axes2 = plt.subplots(1, 4, figsize=(24, 5))

    # y_ctr, z_ctr, r_upper already computed before fig1
    r = r_arr
    t_duration = t_max - t_min
    r_bins = np.linspace(0, r_upper, r_n_bins + 1)
    rh, re = np.histogram(r, bins=r_bins, weights=m_arr)
    rc = (re[:-1] + re[1:]) / 2
    ring_areas = np.pi * (re[1:] ** 2 - re[:-1] ** 2)
    with np.errstate(invalid='ignore', divide='ignore'):
        r_density = np.where(ring_areas > 0, rh / ring_areas * 1e3 / t_duration, 0)

    axes2[0].bar(rc, r_density, width=re[1] - re[0], alpha=0.7, color='purple')
    axes2[0].axvline(r_centroid, color='k', ls='--', lw=1.5, label=f'centroid r={r_centroid:.3f} mm')
    axes2[0].legend(fontsize=8)
    axes2[0].set(xlabel='Radial distance (mm)', ylabel='Mass flux density (g/mm\u00b2/s)',
                 title=f'Radial mass flux density\n(centroid r={r_centroid:.3f} mm)')
    axes2[0].grid(alpha=0.3)

    cum_r = np.cumsum(rh)
    cum_total = cum_r[-1] if cum_r[-1] > 0 else 1.0
    axes2[1].plot(rc, cum_r / cum_total * 100, 'b-', lw=2)
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
    print(f"\n质量加权质心 (beam centre) : Y={y_ctr:.3f} mm, Z={z_ctr:.3f} mm")
    print(f"质心径向距离 r_centroid     : {r_centroid:.3f} mm")
    print(f"50% 质量半径  : {r50:.3f} mm")
    print(f"90% 质量半径  : {r90:.3f} mm")

    # --- 平滑密度曲线 ---
    # 用更密的采样点做 KDE 式平滑：先高密度插值后 Gaussian 滤波
    rc_fine = np.linspace(0, r_upper, max(300, r_n_bins * 10))
    dr_fine = rc_fine[1] - rc_fine[0]
    rh_fine, re_fine = np.histogram(r, bins=np.append(rc_fine - dr_fine / 2,
                                                       rc_fine[-1] + dr_fine / 2),
                                    weights=m_arr)
    ring_fine = np.pi * ((rc_fine + dr_fine / 2) ** 2 - (rc_fine - dr_fine / 2) ** 2)
    with np.errstate(invalid='ignore', divide='ignore'):
        r_density_fine = np.where(ring_fine > 0, rh_fine / ring_fine * 1e3 / t_duration, 0)

    if smooth_sigma > 0:
        # convert sigma from mm to index units
        sigma_idx = smooth_sigma / dr_fine
        r_density_smooth = gaussian_filter1d(r_density_fine, sigma=sigma_idx)
        smooth_label = f'Smoothed (sigma={smooth_sigma:.2f} mm)'
    else:
        r_density_smooth = r_density_fine
        smooth_label = 'Raw (no smoothing)'

    axes2[2].bar(rc, r_density, width=re[1] - re[0],
                 alpha=0.3, color='purple', label='Histogram')
    axes2[2].plot(rc_fine, r_density_smooth, 'r-', lw=2, label=smooth_label)
    axes2[2].axvline(r_centroid, color='k', ls='--', lw=1.5, label=f'centroid r={r_centroid:.3f} mm')
    axes2[2].set(xlabel='Radial distance (mm)', ylabel='Mass flux density (g/mm\u00b2/s)',
                 title=f'Smoothed radial mass flux density\n(centroid r={r_centroid:.3f} mm, sigma={smooth_sigma:.2f} mm)')
    axes2[2].legend(fontsize=8)
    axes2[2].grid(alpha=0.3)

    t_bins = np.linspace(t_arr.min(), t_arr.max(), min(40, len(crossings) + 1))
    if len(t_bins) < 3:
        t_bins = np.linspace(t_min, t_max, 40)
    mh_t, te_t = np.histogram(t_arr, bins=t_bins, weights=m_arr)
    nh_t, _    = np.histogram(t_arr, bins=t_bins)
    tc_t = (te_t[:-1] + te_t[1:]) / 2
    dt_t = te_t[1] - te_t[0]

    ax2r = axes2[3].twinx()
    axes2[3].bar(tc_t * 1000, mh_t / dt_t * 1e3,
                 width=dt_t * 1000, alpha=0.5, color='tomato', label='Mass rate')
    ax2r.plot(tc_t * 1000, nh_t, 'b-o', ms=3, label='Count')
    axes2[3].set(xlabel='Time (ms)', ylabel='Mass flow rate (g/s)')
    axes2[3].yaxis.label.set_color('tomato')
    ax2r.set_ylabel('Crossing count per bin', color='b')
    axes2[3].set_title('Crossing mass rate & count vs time')
    axes2[3].grid(alpha=0.3)
    from matplotlib.lines import Line2D
    axes2[3].legend(
        [Line2D([0], [0], color='tomato', lw=4, alpha=0.5),
         Line2D([0], [0], color='b', marker='o')],
        ['Mass rate (g/s)', 'Count'],
        loc='upper left', fontsize=8
    )

    plt.tight_layout()
    out2 = os.path.join(case_dir, f'crossing_radial_{tag}.png')
    fig2.savefig(out2, dpi=150, bbox_inches='tight')
    print(f"已保存: {out2}")
    plt.close(fig2)

    # =========================================================
    # Figure 3 – 粒子轨迹图
    # =========================================================
    crossing_pids = list({c['pid'] for c in crossings})
    max_plot = 300
    plot_pids = crossing_pids[:max_plot]

    fig3, axes3 = plt.subplots(1, 3, figsize=(30, 6))
    fig3.suptitle(
        f'Particle trajectories crossing x = {x_target_mm:.1f} mm\n'
        f'({len(plot_pids)} of {len(crossing_pids)} crossing particles shown)',
        fontsize=11
    )

    c_times = np.array([
        next((c['t'] for c in crossings if c['pid'] == pid), 0.0)
        for pid in plot_pids
    ])
    norm_c = plt.Normalize(c_times.min(), c_times.max())
    cmap_c = plt.cm.viridis

    for k, pid in enumerate(plot_pids):
        tr = trajectories[pid]
        col = cmap_c(norm_c(c_times[k]))
        axes3[0].plot(np.array(tr['x']) * 1000, np.array(tr['y']) * 1000,
                      lw=0.6, alpha=0.4, color=col)
        axes3[1].plot(np.array(tr['x']) * 1000, np.array(tr['z']) * 1000,
                      lw=0.6, alpha=0.4, color=col)
        axes3[2].plot(np.array(tr['y']) * 1000, np.array(tr['z']) * 1000,
                      lw=0.5, alpha=0.25, color='gray')

    for ax, ylabel, title in [
        (axes3[0], 'Y (mm)', 'XY trajectory projection'),
        (axes3[1], 'Z (mm)', 'XZ trajectory projection'),
    ]:
        ax.axvline(x_target_mm, color='red', lw=1.5, ls='--',
                   label=f'x={x_target_mm:.1f}mm')
        ax.set(xlabel='X (mm)', ylabel=ylabel, title=title)
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(6))
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    sc3 = axes3[2].scatter(
        y_arr, z_arr, c=m_arr / m_arr.max(),
        s=20, cmap='hot_r', alpha=0.85, zorder=5, linewidths=0
    )
    plt.colorbar(sc3, ax=axes3[2], label='Relative mass')
    axes3[2].set(xlabel='Y (mm)', ylabel='Z (mm)',
                 title=f'YZ trajectories + crossing points (x={x_target_mm:.1f}mm)',
                 xlim=(-xy_lim, xy_lim), ylim=(-xy_lim, xy_lim), aspect='equal')
    axes3[2].grid(alpha=0.3)

    sm3 = plt.cm.ScalarMappable(norm=norm_c, cmap=cmap_c)
    sm3.set_array([])
    plt.colorbar(sm3, ax=axes3[0]).set_label('Crossing time (s)')
    plt.colorbar(sm3, ax=axes3[1]).set_label('Crossing time (s)')

    plt.tight_layout()
    out3 = os.path.join(case_dir, f'particle_trajectories_{tag}.png')
    fig3.savefig(out3, dpi=150, bbox_inches='tight')
    print(f"已保存: {out3}")
    plt.close(fig3)

    # =========================================================
    # Figure 4 – 穿越时温度分布（可选）
    # =========================================================
    if np.any(valid_T):
        fig4, axes4 = plt.subplots(1, 3, figsize=(21, 5))
        fig4.suptitle(f'Particle temperature at plane crossing (x={x_target_mm:.1f}mm)')

        sc4 = axes4[0].scatter(
            y_arr[valid_T], z_arr[valid_T], c=T_arr[valid_T],
            s=m_arr[valid_T] / m_arr.max() * 80 + 2,
            cmap='hot', alpha=0.75, linewidths=0
        )
        plt.colorbar(sc4, ax=axes4[0], label='Temperature (K)')
        axes4[0].set(xlabel='Y (mm)', ylabel='Z (mm)',
                     title='Temperature at crossing (size ~ mass)',
                     xlim=(-xy_lim, xy_lim), ylim=(-xy_lim, xy_lim), aspect='equal')
        axes4[0].grid(alpha=0.3)

        Th, Te = np.histogram(T_arr[valid_T], bins=50, weights=m_arr[valid_T])
        Tc = (Te[:-1] + Te[1:]) / 2
        axes4[1].bar(Tc, Th * 1e3, width=Te[1] - Te[0], alpha=0.7, color='orangered')
        T_mean_v = np.average(T_arr[valid_T], weights=m_arr[valid_T])
        axes4[1].axvline(T_mean_v, color='k', ls='--', lw=1.5,
                         label=f'Mean = {T_mean_v:.0f} K')
        axes4[1].set(xlabel='Temperature (K)', ylabel='Mass (g)',
                     title='Mass-weighted temperature distribution at crossing')
        axes4[1].legend()
        axes4[1].grid(alpha=0.3)

        # --- Temperature vs radial position scatter ---
        r_T = r_arr[valid_T]
        axes4[2].scatter(r_T, T_arr[valid_T],
                         s=6, alpha=0.5, color='steelblue', linewidths=0)
        axes4[2].set(xlabel='Radial position (mm)', ylabel='Temperature (K)',
                     title=f'Temperature vs radial position\n(x={x_target_mm:.1f} mm)')
        axes4[2].grid(alpha=0.3)

        plt.tight_layout()
        out4 = os.path.join(case_dir, f'crossing_temperature_{tag}.png')
        fig4.savefig(out4, dpi=150, bbox_inches='tight')
        print(f"已保存: {out4}")
        plt.close(fig4)

    # =========================================================
    # Figure 5 – 多截面平均温度和速度对比 (类似 Fig.11)
    # =========================================================
    x_planes_mm = [11, 13, 15, 17, 19, 20.5, 22, 25, 27, 30.5, 32, 34, 36, 38, 40.5]
    print(f"\n计算多截面统计 ({x_planes_mm} mm) ...")

    T_simple   = []
    T_mw       = []
    spd_simple = []
    spd_mw     = []

    for xp_mm in x_planes_mm:
        xp  = xp_mm / 1000.0
        crs = find_plane_crossings(trajectories, xp, direction='positive')
        if not crs:
            T_simple.append(np.nan);   T_mw.append(np.nan)
            spd_simple.append(np.nan); spd_mw.append(np.nan)
            print(f"  x={xp_mm} mm: 无穿越事件")
            continue

        T_p   = np.array([c['T']     for c in crs])
        spd_p = np.array([c['speed'] for c in crs])
        m_p   = np.array([c['mass']  for c in crs])

        vT   = ~np.isnan(T_p)
        vspd = ~np.isnan(spd_p) & (spd_p > 0)

        T_s  = float(np.nanmean(T_p))                              if np.any(vT)   else np.nan
        T_w  = float(np.average(T_p[vT],   weights=m_p[vT]))      if np.any(vT)   else np.nan
        sp_s = float(np.nanmean(spd_p))                            if np.any(vspd) else np.nan
        sp_w = float(np.average(spd_p[vspd], weights=m_p[vspd]))  if np.any(vspd) else np.nan

        T_simple.append(T_s);   T_mw.append(T_w)
        spd_simple.append(sp_s); spd_mw.append(sp_w)
        print(f"  x={xp_mm} mm: "
              f"T_simple={T_s:.1f} K  T_mw={T_w:.1f} K  "
              f"spd_simple={sp_s:.3f} m/s  spd_mw={sp_w:.3f} m/s  "
              f"({len(crs)} crossings)")

    x_labels = [f'{xp:.1f}' for xp in x_planes_mm]
    x_pos    = np.arange(len(x_planes_mm))
    bar_w    = 0.35

    fig5, axes5 = plt.subplots(1, 2, figsize=(12, 5))
    fig5.suptitle('Averaged temperatures and velocities at transverse planes',
                  fontsize=13)

    # (a) 温度
    ax = axes5[0]
    ax.bar(x_pos - bar_w/2, T_simple, bar_w,
           label='Simple average',       color='steelblue', alpha=0.85)
    ax.bar(x_pos + bar_w/2, T_mw,     bar_w,
           label='Mass-weighted average', color='tomato',    alpha=0.85)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('Averaged temperature (K)')
    ax.set_title('(a) Averaged temperatures')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    ax.set_ylim(bottom=0)

    # (b) 速度
    ax = axes5[1]
    ax.bar(x_pos - bar_w/2, spd_simple, bar_w,
           label='Simple average',       color='steelblue', alpha=0.85)
    ax.bar(x_pos + bar_w/2, spd_mw,     bar_w,
           label='Mass-weighted average', color='tomato',    alpha=0.85)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('Averaged velocity (m/s)')
    ax.set_title('(b) Averaged velocities')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    out5 = os.path.join(case_dir, f'avg_T_speed_planes_{tag}.png')
    fig5.savefig(out5, dpi=150, bbox_inches='tight')
    print(f"已保存: {out5}")
    plt.close(fig5)

    print("\n✓ 分析完成。")


if __name__ == '__main__':
    main()
