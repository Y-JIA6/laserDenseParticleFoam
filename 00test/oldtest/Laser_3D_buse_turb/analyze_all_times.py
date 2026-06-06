#!/usr/bin/env python3
"""
统计 0-0.099 秒内所有时间步在 x=37.5mm 处 YOZ 平面的粒子累积质量分布
Analyze cumulative particle mass distribution at x=37.5mm for all time steps
"""

import json
import struct
import base64
import xml.etree.ElementTree as ET
import zlib

import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import subprocess

def parse_vtp(filename):
    """Parse VTP (vtkPolyData XML) file without vtk library"""
    tree = ET.parse(filename)
    root = tree.getroot()

    # Read global attributes from VTKFile element
    vtk_root = root if root.tag in ('VTKFile', '{VTK}VTKFile') else root.find('.//{*}VTKFile') or root
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

    dtype_map = {'Float32': endian + 'f4', 'Float64': endian + 'f8',
                 'Int32': endian + 'i4', 'Int64': endian + 'i8',
                 'UInt8': 'u1', 'UInt32': endian + 'u4', 'UInt64': endian + 'u8'}

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
            flat = decode_array(da).flatten().astype(np.float64)
            pts = flat.reshape(-1, 3)

    if pts is None:
        return None, None

    # PointData fields
    fields = {}
    pd_node = piece.find('.//{*}PointData') or piece.find('.//PointData')
    if pd_node is not None:
        for da in pd_node.findall('.//{*}DataArray') or pd_node.findall('.//DataArray'):
            name = da.get('Name')
            if name:
                fields[name] = decode_array(da)

    return pts, fields

def calculate_particle_mass(d, rho, nParticle=None):
    """计算粒子质量 m = (4/3) * pi * (d/2)^3 * rho * nParticle"""
    volume = (4.0/3.0) * np.pi * (d/2.0)**3
    mass = volume * rho
    if nParticle is not None:
        mass = mass * nParticle
    return mass

def main():
    case_dir = os.path.dirname(os.path.abspath(__file__))
    vtk_dir = os.path.join(case_dir, 'VTK/lagrangian/kinematicCloud')
    
    # Run foamToVTK to convert OpenFOAM data to VTK format
    run_vtk = input("Run foamToVTK? [1=yes / 0=skip]: ").strip()
    if run_vtk == '1':
        print("Running foamToVTK...")
        result = subprocess.run(
            'foamToVTK 2>&1 | tail -30',
            shell=True,
            cwd=case_dir,
            capture_output=True,
            text=True
        )
        print(result.stdout)
        if result.returncode != 0:
            print(f"foamToVTK exited with code {result.returncode}")
            print(result.stderr)
    else:
        print("foamToVTK skipped.")
    
    print("=" * 70)
    print("Cumulative Particle Mass Distribution Analysis (0-0.099s)")
    print("=" * 70)
    
    # 从 .vtp.series 文件读取精确的时间-文件名映射
    series_file = os.path.join(vtk_dir, 'kinematicCloud.vtp.series')
    with open(series_file, 'r') as f:
        series = json.load(f)
    all_entries = [(e['time'], os.path.join(vtk_dir, e['name'])) for e in series['files']]
    t_min = float(input("Start time (s) [e.g. 0.0]: ").strip() or "0.0")
    t_max = float(input("End   time (s) [e.g. 0.1]: ").strip() or "0.1")
    filtered_files = sorted((t, p) for t, p in all_entries if t_min < t <= t_max)
    print(f"Found {len(filtered_files)} time steps from series file")
    
    # 累积数据
    x_target = 0.0315  # 37.5 mm
    x_tolerance = 0.001  # ±3 mm
    
    all_y = []
    all_z = []
    all_mass = []
    all_T = []
    all_times = []
    
    stats_per_time = []
    
    for time, vtk_file in filtered_files:
        points, fields = parse_vtp(vtk_file)
        
        if points is None or len(points) == 0:
            continue
        
        x = points[:, 0]
        y = points[:, 1]
        z = points[:, 2]
        
        # 筛选x位置
        mask = np.abs(x - x_target) < x_tolerance
        n_selected = np.sum(mask)
        
        if n_selected == 0:
            continue
        
        y_sel = y[mask]
        z_sel = z[mask]
        
        # 计算质量
        if 'd' in fields and 'rho' in fields:
            d_sel = fields['d'][mask]
            rho_sel = fields['rho'][mask]
            nP = fields.get('nParticle', np.ones_like(d_sel))
            if 'nParticle' in fields:
                nP_sel = nP[mask]
            else:
                nP_sel = np.ones_like(d_sel)
            mass_sel = calculate_particle_mass(d_sel, rho_sel, nP_sel)
        else:
            mass_sel = np.ones(n_selected)
        
        # 温度
        T_sel = None
        if 'T_particle' in fields:
            T_sel = fields['T_particle'][mask]
        
        all_y.extend(y_sel * 1000)  # mm
        all_z.extend(z_sel * 1000)  # mm
        all_mass.extend(mass_sel)
        if T_sel is not None:
            all_T.extend(T_sel)
        all_times.extend([time] * n_selected)
        
        stats_per_time.append({
            'time': time,
            'n_particles': n_selected,
            'total_mass': np.sum(mass_sel)
        })
    
    # 转换为numpy数组
    all_y = np.array(all_y)
    all_z = np.array(all_z)
    all_mass = np.array(all_mass)
    all_T = np.array(all_T) if len(all_T) > 0 else None
    all_times = np.array(all_times)
    
    print(f"\n=== Cumulative Statistics at x = {x_target*1000:.1f} mm ===")
    print(f"Total particle-time samples: {len(all_y)}")
    print(f"Total cumulative mass: {np.sum(all_mass)*1e9:.3f} ng")
    print(f"Y range: {all_y.min():.2f} - {all_y.max():.2f} mm")
    print(f"Z range: {all_z.min():.2f} - {all_z.max():.2f} mm")
    if all_T is not None and len(all_T) > 0:
        print(f"Temperature range: {all_T.min():.1f} - {all_T.max():.1f} K")
    
    # ========== 绘图 ==========
    
    # 图1: 累积质量分布（2x2）
    fig1 = plt.figure(figsize=(16, 14))
    
    # 1.1 散点图（按时间着色）
    ax1 = fig1.add_subplot(2, 2, 1)
    scatter = ax1.scatter(all_y, all_z, c=all_times*1000, s=all_mass/all_mass.max()*50+2, 
                         cmap='viridis', alpha=0.5)
    cbar = plt.colorbar(scatter, ax=ax1)
    cbar.set_label('Time (ms)')
    ax1.set_xlabel('Y (mm)')
    ax1.set_ylabel('Z (mm)')
    ax1.set_title(f'All Particle Positions (0-99ms) at x={x_target*1000:.1f}mm\n({len(all_y)} samples)')
    ax1.set_aspect('equal')
    lim = max(np.abs(all_y).max(), np.abs(all_z).max()) * 1.05
    ax1.set_xlim(-lim, lim)
    ax1.set_ylim(-lim, lim)
    ax1.grid(True, alpha=0.3)
    
    # 1.2 2D质量密度直方图
    ax2 = fig1.add_subplot(2, 2, 2)
    n_bins = 40
    y_bins = np.linspace(-lim, lim, n_bins)
    z_bins = np.linspace(-lim, lim, n_bins)
    H, ye, ze = np.histogram2d(all_y, all_z, bins=[y_bins, z_bins], weights=all_mass)
    H_pct = H / np.sum(all_mass) * 100  # percentage of total mass per bin
    
    im = ax2.imshow(H_pct.T, origin='lower', aspect='auto',
                   extent=[ye[0], ye[-1], ze[0], ze[-1]], cmap='YlOrRd')
    cbar2 = plt.colorbar(im, ax=ax2)
    cbar2.set_label('Mass fraction per bin (%)')
    ax2.set_xlabel('Y (mm)')
    ax2.set_ylabel('Z (mm)')
    ax2.set_xlim(-lim, lim)
    ax2.set_ylim(-lim, lim)
    ax2.set_title('Cumulative Mass Distribution (0-99ms) [% of total]')
    
    # 1.3 Y方向质量分布
    ax3 = fig1.add_subplot(2, 2, 3)
    y_hist, y_edges = np.histogram(all_y, bins=50, weights=all_mass)
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2
    ax3.bar(y_centers, y_hist*1e9, width=y_edges[1]-y_edges[0], alpha=0.7, color='steelblue')
    ax3.set_xlabel('Y (mm)')
    ax3.set_ylabel('Cumulative Mass (ng)')
    ax3.set_title('Cumulative Mass Distribution along Y axis')
    ax3.set_xlim(-lim, lim)
    ax3.grid(True, alpha=0.3)
    
    # 1.4 Z方向质量分布
    ax4 = fig1.add_subplot(2, 2, 4)
    z_hist, z_edges = np.histogram(all_z, bins=50, weights=all_mass)
    z_centers = (z_edges[:-1] + z_edges[1:]) / 2
    ax4.bar(z_centers, z_hist*1e9, width=z_edges[1]-z_edges[0], alpha=0.7, color='seagreen')
    ax4.set_xlabel('Z (mm)')
    ax4.set_ylabel('Cumulative Mass (ng)')
    ax4.set_title('Cumulative Mass Distribution along Z axis')
    ax4.set_xlim(-lim, lim)
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    out1 = 'cumulative_yoz_distribution_0-99ms.png'
    plt.savefig(out1, dpi=150)
    print(f"\nSaved: {out1}")
    plt.close()
    
    # 图2: 径向分布
    fig2, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    y_center = np.average(all_y, weights=all_mass)
    z_center = np.average(all_z, weights=all_mass)
    r = np.sqrt((all_y - y_center)**2 + (all_z - z_center)**2)
    
    # 2.1 径向质量密度
    r_bins = np.linspace(0, r.max(), 40)
    r_hist, r_edges = np.histogram(r, bins=r_bins, weights=all_mass)
    r_centers = (r_edges[:-1] + r_edges[1:]) / 2
    ring_areas = np.pi * (r_edges[1:]**2 - r_edges[:-1]**2)
    r_density = r_hist / ring_areas * 1e9
    
    axes[0].bar(r_centers, r_density, width=r_edges[1]-r_edges[0], alpha=0.7, color='purple')
    axes[0].set_xlabel('Radial distance (mm)')
    axes[0].set_ylabel('Mass density (ng/mm²)')
    axes[0].set_title(f'Radial Mass Density\n(Center: Y={y_center:.2f}, Z={z_center:.2f} mm)')
    axes[0].grid(True, alpha=0.3)
    
    # 2.2 累积径向分布
    cumulative_mass = np.cumsum(r_hist)
    axes[1].plot(r_centers, cumulative_mass/cumulative_mass[-1]*100, 'b-', linewidth=2)
    axes[1].axhline(y=50, color='r', linestyle='--', alpha=0.7, label='50%')
    axes[1].axhline(y=90, color='g', linestyle='--', alpha=0.7, label='90%')
    axes[1].set_xlabel('Radial distance (mm)')
    axes[1].set_ylabel('Cumulative mass (%)')
    axes[1].set_title('Cumulative Radial Mass Distribution')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    # 找50%和90%半径
    r50_idx = np.argmax(cumulative_mass >= 0.5 * cumulative_mass[-1])
    r90_idx = np.argmax(cumulative_mass >= 0.9 * cumulative_mass[-1])
    print(f"Mass-weighted centroid: Y={y_center:.3f}mm, Z={z_center:.3f}mm")
    print(f"50% mass radius: {r_centers[r50_idx]:.3f} mm")
    print(f"90% mass radius: {r_centers[r90_idx]:.3f} mm")
    
    # 2.3 每个时间步的粒子数和质量
    times = [s['time']*1000 for s in stats_per_time]
    n_particles = [s['n_particles'] for s in stats_per_time]
    masses = [s['total_mass']*1e9 for s in stats_per_time]
    
    ax3_twin = axes[2].twinx()
    l1 = axes[2].plot(times, n_particles, 'b-o', markersize=2, label='Particles')
    l2 = ax3_twin.plot(times, masses, 'r-s', markersize=2, label='Mass (ng)')
    axes[2].set_xlabel('Time (ms)')
    axes[2].set_ylabel('Number of particles', color='b')
    ax3_twin.set_ylabel('Mass (ng)', color='r')
    axes[2].set_title('Particle Count & Mass vs Time')
    axes[2].grid(True, alpha=0.3)
    
    lines = l1 + l2
    labels = [l.get_label() for l in lines]
    axes[2].legend(lines, labels, loc='upper left')
    
    plt.tight_layout()
    out2 = 'cumulative_radial_distribution_0-99ms.png'
    plt.savefig(out2, dpi=150)
    print(f"Saved: {out2}")
    plt.close()
    
    # 图3: 温度分布（如果有）
    if all_T is not None and len(all_T) > 0:
        fig3, axes3 = plt.subplots(1, 2, figsize=(14, 5))
        
        # 按温度着色的散点图
        scatter3 = axes3[0].scatter(all_y, all_z, c=all_T, s=all_mass/all_mass.max()*50+2,
                                    cmap='hot', alpha=0.5)
        cbar3 = plt.colorbar(scatter3, ax=axes3[0])
        cbar3.set_label('Temperature (K)')
        axes3[0].set_xlabel('Y (mm)')
        axes3[0].set_ylabel('Z (mm)')
        axes3[0].set_title('Particle Temperature Distribution')
        axes3[0].set_aspect('equal')
        lim_t = max(np.abs(all_y).max(), np.abs(all_z).max()) * 1.05
        axes3[0].set_xlim(-lim_t, lim_t)
        axes3[0].set_ylim(-lim_t, lim_t)
        axes3[0].grid(True, alpha=0.3)
        
        # 温度直方图
        T_hist, T_edges = np.histogram(all_T, bins=50, weights=all_mass)
        T_centers = (T_edges[:-1] + T_edges[1:]) / 2
        axes3[1].bar(T_centers, T_hist*1e9, width=T_edges[1]-T_edges[0], alpha=0.7, color='orangered')
        axes3[1].set_xlabel('Temperature (K)')
        axes3[1].set_ylabel('Mass (ng)')
        axes3[1].set_title('Mass-weighted Temperature Distribution')
        axes3[1].grid(True, alpha=0.3)
        
        T_mean = np.average(all_T, weights=all_mass)
        print(f"Mass-weighted mean temperature: {T_mean:.1f} K")
        
        plt.tight_layout()
        out3 = 'cumulative_temperature_distribution_0-99ms.png'
        plt.savefig(out3, dpi=150)
        print(f"Saved: {out3}")
        plt.close()
    
    print("\nAnalysis complete!")

if __name__ == "__main__":
    main()
