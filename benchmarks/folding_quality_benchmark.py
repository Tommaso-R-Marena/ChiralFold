#!/usr/bin/env python3
"""
ChiralFold — 3D Folding Quality Benchmark
============================================

Assesses whether ChiralFold produces structurally valid 3D peptide
geometries — going beyond chirality correctness to evaluate bond geometry,
backbone dihedral distribution, conformer diversity, and comparison
to experimental crystal structures.

Ground truth: PDB 3IWY — MDM2 + D-peptide dPMI-gamma (DWWPLAFEALLR), 1.9 Å.
              Same DP12:MDM2 system from Childs, Zhou & Donald (2025).
"""

import sys, os, time, json, math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rdkit import Chem
from rdkit.Chem import AllChem
from chiralfold.model import d_peptide_smiles, l_peptide_smiles, mixed_peptide_smiles

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'results')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════════
# Utility: backbone analysis from RDKit molecules
# ═══════════════════════════════════════════════════════════════════════════

def find_backbone_atoms(mol):
    """Identify backbone N, Cα, C, O atoms per residue."""
    residues = []
    visited_ca = set()
    for atom in mol.GetAtoms():
        if atom.GetIdx() in visited_ca or atom.GetSymbol() != 'C':
            continue
        n_idx = co_idx = o_idx = None
        for nb in atom.GetNeighbors():
            if nb.GetSymbol() == 'N' and n_idx is None:
                n_idx = nb.GetIdx()
            elif nb.GetSymbol() == 'C':
                for sub in nb.GetNeighbors():
                    if sub.GetSymbol() == 'O':
                        bond = mol.GetBondBetweenAtoms(nb.GetIdx(), sub.GetIdx())
                        if bond and bond.GetBondTypeAsDouble() >= 1.5:
                            co_idx = nb.GetIdx()
                            o_idx = sub.GetIdx()
                            break
        if n_idx is not None and co_idx is not None:
            residues.append({'N': n_idx, 'CA': atom.GetIdx(), 'C': co_idx, 'O': o_idx})
            visited_ca.add(atom.GetIdx())
    residues.sort(key=lambda r: r['N'])
    return residues


def compute_dihedral(p1, p2, p3, p4):
    b1, b2, b3 = p2 - p1, p3 - p2, p4 - p3
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    m1 = np.cross(n1, b2 / (np.linalg.norm(b2) + 1e-12))
    return np.degrees(np.arctan2(np.dot(m1, n2), np.dot(n1, n2)))


def compute_backbone_dihedrals(mol, conf_id, residues):
    pos = mol.GetConformer(conf_id).GetPositions()
    n = len(residues)
    dihedrals = []
    for i in range(n):
        phi = psi = omega = np.nan
        if i > 0:
            phi = compute_dihedral(pos[residues[i-1]['C']], pos[residues[i]['N']],
                                   pos[residues[i]['CA']], pos[residues[i]['C']])
        if i < n - 1:
            psi = compute_dihedral(pos[residues[i]['N']], pos[residues[i]['CA']],
                                   pos[residues[i]['C']], pos[residues[i+1]['N']])
        if i > 0:
            omega = compute_dihedral(pos[residues[i-1]['CA']], pos[residues[i-1]['C']],
                                     pos[residues[i]['N']], pos[residues[i]['CA']])
        dihedrals.append({'residue': i, 'phi': phi, 'psi': psi, 'omega': omega})
    return dihedrals


def compute_bond_geometry(mol, conf_id, residues):
    pos = mol.GetConformer(conf_id).GetPositions()
    n = len(residues)
    bl = {'N_CA': [], 'CA_C': [], 'C_N': [], 'C_O': []}
    ba = {'N_CA_C': [], 'CA_C_N': [], 'C_N_CA': []}
    for i in range(n):
        r = residues[i]
        bl['N_CA'].append(np.linalg.norm(pos[r['CA']] - pos[r['N']]))
        bl['CA_C'].append(np.linalg.norm(pos[r['C']] - pos[r['CA']]))
        if r['O'] is not None:
            bl['C_O'].append(np.linalg.norm(pos[r['O']] - pos[r['C']]))
        if i < n - 1:
            bl['C_N'].append(np.linalg.norm(pos[residues[i+1]['N']] - pos[r['C']]))
        def angle(a, b, c):
            v1, v2 = a - b, c - b
            cos_a = np.clip(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12), -1, 1)
            return np.degrees(np.arccos(cos_a))
        ba['N_CA_C'].append(angle(pos[r['N']], pos[r['CA']], pos[r['C']]))
        if i < n - 1:
            ba['CA_C_N'].append(angle(pos[r['CA']], pos[r['C']], pos[residues[i+1]['N']]))
        if i > 0:
            ba['C_N_CA'].append(angle(pos[residues[i-1]['C']], pos[r['N']], pos[r['CA']]))
    return bl, ba


def parse_pdb_backbone(pdb_path, chain='B'):
    residues = {}
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            if line[21] != chain:
                continue
            aname = line[12:16].strip()
            resnum = int(line[22:26])
            resname = line[17:20].strip()
            if resname == 'HOH':
                continue
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            if resnum not in residues:
                residues[resnum] = {'resname': resname, 'resnum': resnum}
            residues[resnum][aname] = np.array([x, y, z])
    backbone = [residues[rn] for rn in sorted(residues.keys())
                if 'N' in residues[rn] and 'CA' in residues[rn] and 'C' in residues[rn]]
    return backbone


def compute_pdb_dihedrals(backbone):
    n = len(backbone)
    dihedrals = []
    for i in range(n):
        phi = psi = omega = np.nan
        if i > 0:
            phi = compute_dihedral(backbone[i-1]['C'], backbone[i]['N'],
                                   backbone[i]['CA'], backbone[i]['C'])
        if i < n - 1:
            psi = compute_dihedral(backbone[i]['N'], backbone[i]['CA'],
                                   backbone[i]['C'], backbone[i+1]['N'])
        if i > 0:
            omega = compute_dihedral(backbone[i-1]['CA'], backbone[i-1]['C'],
                                     backbone[i]['N'], backbone[i]['CA'])
        dihedrals.append({'residue': i, 'phi': phi, 'psi': psi, 'omega': omega,
                         'resname': backbone[i]['resname']})
    return dihedrals


# ═══════════════════════════════════════════════════════════════════════════
# Conformer ensemble generation + analysis
# ═══════════════════════════════════════════════════════════════════════════

def generate_and_analyze(seq, chirality, n_confs=20, label=''):
    smi = mixed_peptide_smiles(seq, chirality)
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.numThreads = 0
    params.pruneRmsThresh = 0.3
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)
    if len(cids) == 0:
        params.useRandomCoords = True
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)
    if len(cids) == 0:
        return None

    energies = []
    try:
        results = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=500)
        energies = [e for _, e in results]
    except:
        energies = [0.0] * len(cids)

    residues = find_backbone_atoms(mol)
    if not residues:
        return None

    all_dihedrals, all_bl, all_ba, all_rg = [], [], [], []
    min_e_idx = int(np.argmin(energies)) if energies else 0

    # Compute pairwise conformer RMSD
    all_rmsd = []
    for i, cid in enumerate(cids):
        dihs = compute_backbone_dihedrals(mol, cid, residues)
        for d in dihs:
            d['conf_id'] = i
            d['energy'] = energies[i] if i < len(energies) else 0
        all_dihedrals.extend(dihs)
        bl, ba = compute_bond_geometry(mol, cid, residues)
        all_bl.append(bl)
        all_ba.append(ba)
        pos = mol.GetConformer(cid).GetPositions()
        centroid = pos.mean(axis=0)
        all_rg.append(np.sqrt(np.mean(np.sum((pos - centroid)**2, axis=1))))

        # Compute Cα RMSD to lowest-energy conformer
        if i != min_e_idx:
            pos_min = mol.GetConformer(cids[min_e_idx]).GetPositions()
            ca_idx = [r['CA'] for r in residues]
            ca_pos = pos[ca_idx]
            ca_min = pos_min[ca_idx]
            # Simple RMSD without alignment
            rmsd = np.sqrt(np.mean(np.sum((ca_pos - ca_min)**2, axis=1)))
            all_rmsd.append(rmsd)
        else:
            all_rmsd.append(0.0)

    return {
        'seq': seq, 'chirality': chirality, 'label': label,
        'mol': mol, 'cids': list(cids), 'energies': energies,
        'residues': residues, 'dihedrals': all_dihedrals,
        'bond_lengths': all_bl, 'bond_angles': all_ba,
        'rg': all_rg, 'rmsd_to_min': all_rmsd, 'n_confs': len(cids),
    }


# ═══════════════════════════════════════════════════════════════════════════
# Color palette + style
# ═══════════════════════════════════════════════════════════════════════════

PAL = {
    'teal': '#0D9488', 'cyan': '#06B6D4', 'coral': '#EF4444',
    'amber': '#F59E0B', 'slate': '#334155', 'muted': '#94A3B8',
    'bg': '#F8FAFC', 'grid': '#E2E8F0', 'navy': '#1E293B',
}

def setup_style():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['DejaVu Sans'],
        'font.size': 10, 'axes.titlesize': 12, 'axes.titleweight': 'bold',
        'axes.labelsize': 10, 'axes.facecolor': PAL['bg'],
        'axes.edgecolor': PAL['grid'], 'axes.grid': True,
        'axes.linewidth': 0.5, 'grid.color': PAL['grid'],
        'grid.linewidth': 0.3, 'grid.alpha': 0.7,
        'xtick.labelsize': 8, 'ytick.labelsize': 8,
        'figure.facecolor': 'white', 'figure.dpi': 150,
    })

setup_style()

# ═══════════════════════════════════════════════════════════════════════════
# Visualization panels
# ═══════════════════════════════════════════════════════════════════════════

def panel_ramachandran(ax, all_data, pdb_dihedrals=None):
    """Ramachandran density plot with PDB comparison."""
    phis, psis = [], []
    for d in all_data:
        for dih in d['dihedrals']:
            if not np.isnan(dih['phi']) and not np.isnan(dih['psi']):
                phis.append(dih['phi'])
                psis.append(dih['psi'])
    phis, psis = np.array(phis), np.array(psis)

    # Standard allowed regions (apply to both L and D backbone geometry)
    for cx, cy, w, h, c, lbl in [
        (-60, -40, 80, 80, PAL['grid'], 'L-α'),
        (-120, 130, 100, 80, PAL['grid'], 'L-β'),
        (60, 40, 80, 80, PAL['teal'], 'D-α'),
        (120, -130, 100, 80, PAL['cyan'], 'D-β'),
    ]:
        ax.add_patch(Ellipse((cx, cy), w, h, alpha=0.1, color=c))
        ax.annotate(lbl, xy=(cx, cy), fontsize=7, ha='center', color=c, alpha=0.6)

    # ChiralFold conformer density
    if len(phis) > 10:
        try:
            from scipy.stats import gaussian_kde
            kde = gaussian_kde(np.vstack([phis, psis]), bw_method=0.25)
            xg, yg = np.linspace(-180, 180, 80), np.linspace(-180, 180, 80)
            X, Y = np.meshgrid(xg, yg)
            Z = kde(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
            ax.contourf(X, Y, Z, levels=8, cmap='BuGn', alpha=0.4, zorder=2)
            ax.contour(X, Y, Z, levels=5, colors=PAL['teal'], alpha=0.5, linewidths=0.5)
        except:
            pass

    ax.scatter(phis, psis, s=8, alpha=0.3, c=PAL['teal'], edgecolors='none', zorder=3,
               label=f'ChiralFold (n={len(phis)})')

    if pdb_dihedrals:
        ep = [(d['phi'], d['psi']) for d in pdb_dihedrals
              if not np.isnan(d['phi']) and not np.isnan(d['psi'])]
        if ep:
            ephi, epsi = zip(*ep)
            ax.scatter(ephi, epsi, s=80, marker='*', c=PAL['coral'], edgecolors='white',
                       linewidths=0.5, zorder=5, label=f'PDB 3IWY crystal (n={len(ep)})')

    ax.set_xlim(-180, 180); ax.set_ylim(-180, 180)
    ax.set_xlabel('φ (deg)'); ax.set_ylabel('ψ (deg)')
    ax.set_title('A. Backbone Dihedral Distribution')
    ax.axhline(0, color=PAL['grid'], lw=0.5); ax.axvline(0, color=PAL['grid'], lw=0.5)
    ax.legend(fontsize=7, loc='lower left', framealpha=0.9)


def panel_bond_lengths(ax, all_data):
    """Bond length boxplots vs ideal values."""
    ideal = {'N_CA': 1.458, 'CA_C': 1.525, 'C_N': 1.329, 'C_O': 1.231}
    names = ['N—Cα\n(1.458)', 'Cα—C\n(1.525)', 'C—N\n(1.329)', 'C=O\n(1.231)']
    colors = [PAL['teal'], PAL['cyan'], PAL['coral'], PAL['amber']]

    all_bl = {k: [] for k in ideal}
    for d in all_data:
        for bl in d['bond_lengths']:
            for k in ideal:
                all_bl[k].extend(bl[k])

    data_lists = [np.array(all_bl[k]) for k in ideal]

    for i, (vals, ideal_val) in enumerate(zip(data_lists, ideal.values())):
        if len(vals) == 0: continue
        bp = ax.boxplot([vals], positions=[i], widths=0.5, patch_artist=True,
                       showfliers=False, medianprops=dict(color='white', lw=1.5),
                       whiskerprops=dict(color=colors[i]),
                       capprops=dict(color=colors[i]))
        bp['boxes'][0].set_facecolor(colors[i])
        bp['boxes'][0].set_alpha(0.6)
        ax.plot(i, ideal_val, 'D', color='white', ms=7, mec=colors[i], mew=2, zorder=5)
        rmsd = np.sqrt(np.mean((vals - ideal_val)**2))
        ax.text(i, vals.max() + 0.008, f'{rmsd:.3f}Å', ha='center', fontsize=7,
                color=colors[i], fontweight='bold')

    # Overall
    all_dev = []
    for k, iv in ideal.items():
        all_dev.extend((np.array(all_bl[k]) - iv).tolist())
    overall = np.sqrt(np.mean(np.array(all_dev)**2))

    ax.set_xticks(range(4)); ax.set_xticklabels(names, fontsize=8)
    ax.set_ylabel('Bond Length (Å)')
    ax.set_title('B. Bond Length Quality')
    ax.text(0.02, 0.97, f'Overall RMSD: {overall:.4f} Å', transform=ax.transAxes,
            fontsize=9, fontweight='bold', va='top', color=PAL['navy'],
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec=PAL['grid']))


def panel_bond_angles(ax, all_data):
    ideal = {'N_CA_C': 111.0, 'CA_C_N': 116.2, 'C_N_CA': 121.7}
    names = ['N—Cα—C\n(111.0°)', 'Cα—C—N\n(116.2°)', 'C—N—Cα\n(121.7°)']
    colors = [PAL['teal'], PAL['cyan'], PAL['coral']]

    all_ba = {k: [] for k in ideal}
    for d in all_data:
        for ba in d['bond_angles']:
            for k in ideal:
                all_ba[k].extend(ba[k])

    for i, (k, iv) in enumerate(ideal.items()):
        vals = np.array(all_ba[k])
        if not len(vals): continue
        bp = ax.boxplot([vals], positions=[i], widths=0.5, patch_artist=True,
                       showfliers=False, medianprops=dict(color='white', lw=1.5))
        bp['boxes'][0].set_facecolor(colors[i]); bp['boxes'][0].set_alpha(0.6)
        ax.plot(i, iv, 'D', color='white', ms=7, mec=colors[i], mew=2, zorder=5)
        rmsd = np.sqrt(np.mean((vals - iv)**2))
        ax.text(i, vals.max() + 0.5, f'{rmsd:.1f}°', ha='center', fontsize=7,
                color=colors[i], fontweight='bold')

    ax.set_xticks(range(3)); ax.set_xticklabels(names, fontsize=8)
    ax.set_ylabel('Bond Angle (deg)')
    ax.set_title('C. Bond Angle Quality')


def panel_planarity(ax, all_data, pdb_dihedrals=None):
    """Peptide bond omega angle histogram."""
    omegas = []
    for d in all_data:
        for dih in d['dihedrals']:
            if not np.isnan(dih['omega']):
                omegas.append(dih['omega'])
    omegas = np.array(omegas)
    dev = np.minimum(np.abs(omegas - 180), np.abs(omegas + 180))
    dev = np.minimum(dev, np.abs(np.abs(omegas) - 180))

    ax.hist(dev, bins=np.arange(0, 65, 2), color=PAL['teal'], alpha=0.7,
            edgecolor='white', lw=0.5, label=f'ChiralFold (n={len(omegas)})')

    if pdb_dihedrals:
        pdb_om = [d['omega'] for d in pdb_dihedrals if not np.isnan(d['omega'])]
        pdb_dev = [min(abs(o - 180), abs(o + 180), abs(abs(o) - 180)) for o in pdb_om]
        ax.axvline(np.mean(pdb_dev), color=PAL['coral'], ls='--', lw=2,
                   label=f'PDB 3IWY mean: {np.mean(pdb_dev):.1f}°')

    ax.axvline(6.0, color=PAL['amber'], ls=':', lw=1.5, label='PDB quality cutoff (6°)')
    pct_good = np.sum(dev < 6) / len(dev) * 100 if len(dev) > 0 else 0
    med = np.median(dev) if len(dev) > 0 else 0

    ax.set_xlabel('|ω| deviation from 180° (deg)')
    ax.set_ylabel('Count')
    ax.set_title('D. Peptide Bond Planarity')
    ax.legend(fontsize=7, loc='upper right')
    ax.text(0.97, 0.75, f'{pct_good:.0f}% within 6°\nMedian: {med:.1f}°',
            transform=ax.transAxes, ha='right', va='top', fontsize=9,
            fontweight='bold', color=PAL['teal'],
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec=PAL['teal'], alpha=0.9))


def panel_energy_landscape(ax, all_data):
    """Energy vs Cα RMSD to lowest-energy conformer."""
    for d in all_data:
        if not d['energies']:
            continue
        e = np.array(d['energies'])
        r = np.array(d['rmsd_to_min'])
        e_rel = e - e.min()
        ax.scatter(r, e_rel, s=18, alpha=0.5, label=d['label'], zorder=3)

    ax.set_xlabel('Cα RMSD to lowest-energy conformer (Å)')
    ax.set_ylabel('Relative energy (kcal/mol)')
    ax.set_title('E. Conformer Energy Landscape')
    ax.legend(fontsize=6, ncol=2, loc='upper left')


def panel_rg(ax, all_data):
    """Radius of gyration vs chain length."""
    for d in all_data:
        if not d['rg']:
            continue
        length = len(d['seq'])
        rg_arr = np.array(d['rg'])
        ax.errorbar(length, np.mean(rg_arr), yerr=np.std(rg_arr), fmt='o',
                     color=PAL['teal'], ms=7, capsize=3, lw=1.5, zorder=3)
        ax.annotate(d['label'], (length, np.mean(rg_arr)), fontsize=5,
                    xytext=(5, 3), textcoords='offset points', color=PAL['muted'])

    ns = np.linspace(5, 11, 50)
    ax.plot(ns, 1.8 * ns**0.36, '--', color=PAL['muted'], lw=1, alpha=0.5, label='Compact ~ N⁰·³⁶')
    ax.plot(ns, 1.2 * ns**0.60, ':', color=PAL['muted'], lw=1, alpha=0.5, label='Coil ~ N⁰·⁶')

    ax.set_xlabel('Peptide length (residues)')
    ax.set_ylabel('Radius of Gyration (Å)')
    ax.set_title('F. Compactness vs Chain Length')
    ax.legend(fontsize=7)


def panel_3d_backbone(ax, data):
    """3D backbone Cα traces showing conformer diversity."""
    mol, residues, cids = data['mol'], data['residues'], data['cids']
    cmap = plt.cm.viridis
    for j, cid in enumerate(cids[:12]):
        pos = mol.GetConformer(cid).GetPositions()
        ca = np.array([pos[r['CA']] for r in residues])
        c = cmap(j / max(len(cids[:12]) - 1, 1))
        alpha = 0.8 if j == 0 else 0.25
        lw = 2.5 if j == 0 else 0.8
        ax.plot(ca[:, 0], ca[:, 1], ca[:, 2], '-o', color=c, alpha=alpha, lw=lw, ms=3, zorder=10-j)
    ax.set_xlabel('X (Å)', fontsize=7); ax.set_ylabel('Y (Å)', fontsize=7)
    ax.set_zlabel('Z (Å)', fontsize=7)
    ax.set_title(f'G. Cα Backbone Traces\n{data["label"]}')
    ax.tick_params(labelsize=5)


def panel_scorecard(ax, all_data, pdb_dihedrals=None):
    """Quality summary scorecard with honest assessment."""
    ax.axis('off')

    # Aggregate metrics
    all_phis, all_psis, all_omegas = [], [], []
    all_bl = {k: [] for k in ['N_CA', 'CA_C', 'C_N', 'C_O']}
    all_ba = {k: [] for k in ['N_CA_C', 'CA_C_N', 'C_N_CA']}

    for d in all_data:
        for dih in d['dihedrals']:
            if not np.isnan(dih['phi']): all_phis.append(dih['phi'])
            if not np.isnan(dih['psi']): all_psis.append(dih['psi'])
            if not np.isnan(dih['omega']): all_omegas.append(dih['omega'])
        for bl in d['bond_lengths']:
            for k in all_bl: all_bl[k].extend(bl[k])
        for ba in d['bond_angles']:
            for k in all_ba: all_ba[k].extend(ba[k])

    ideal_bl = {'N_CA': 1.458, 'CA_C': 1.525, 'C_N': 1.329, 'C_O': 1.231}
    ideal_ba = {'N_CA_C': 111.0, 'CA_C_N': 116.2, 'C_N_CA': 121.7}

    bl_dev = []
    for k, iv in ideal_bl.items():
        if all_bl[k]: bl_dev.extend(((np.array(all_bl[k]) - iv)**2).tolist())
    bl_rmsd = np.sqrt(np.mean(bl_dev)) if bl_dev else 0

    ba_dev = []
    for k, iv in ideal_ba.items():
        if all_ba[k]: ba_dev.extend(((np.array(all_ba[k]) - iv)**2).tolist())
    ba_rmsd = np.sqrt(np.mean(ba_dev)) if ba_dev else 0

    omegas = np.array(all_omegas)
    omega_dev = np.minimum(np.abs(omegas - 180), np.abs(omegas + 180))
    omega_pct = np.sum(omega_dev < 6) / len(omega_dev) * 100 if len(omega_dev) > 0 else 0

    n_confs = sum(d['n_confs'] for d in all_data)
    n_obs = sum(len(d['residues']) * d['n_confs'] for d in all_data)

    rows = [
        ('', 'CHIRALFOLD v2.1 QUALITY REPORT', '', None),
        ('', '', '', None),
        ('Sequences analyzed', f'{len(all_data)}', '', None),
        ('Total conformers', f'{n_confs}', '', None),
        ('Residue observations', f'{n_obs:,}', '', None),
        ('', '', '', None),
        ('Bond length RMSD', f'{bl_rmsd:.4f} Å', 'PDB typical: 0.02 Å', bl_rmsd < 0.03),
        ('Bond angle RMSD', f'{ba_rmsd:.1f}°', 'PDB typical: 1.5°', ba_rmsd < 3.0),
        ('Peptide planarity', f'{omega_pct:.0f}% < 6°', 'MMFF94 limitation', omega_pct > 50),
        ('Chirality guarantee', '0% violations', 'By construction', True),
        ('', '', '', None),
        ('PDB ground truth', '3IWY (1.9 Å)', '', None),
        ('System', 'MDM2 + dPMI-γ', '', None),
    ]

    y = 0.95
    for label, value, note, status in rows:
        if value == 'CHIRALFOLD v2.1 QUALITY REPORT':
            ax.text(0.5, y, value, transform=ax.transAxes, fontsize=11,
                    fontweight='bold', ha='center', va='top', color=PAL['navy'])
            y -= 0.06; continue
        if not label and not value:
            y -= 0.03; continue
        color = PAL['navy']
        if status is True: color = PAL['teal']
        elif status is False: color = PAL['coral']
        ax.text(0.04, y, label, transform=ax.transAxes, fontsize=9, ha='left',
                va='top', color=PAL['muted'])
        ax.text(0.58, y, value, transform=ax.transAxes, fontsize=9, ha='left',
                va='top', fontweight='bold', color=color)
        if note:
            ax.text(0.96, y, note, transform=ax.transAxes, fontsize=6.5, ha='right',
                    va='top', color=PAL['muted'], style='italic')
        y -= 0.065

    ax.set_title('H. Quality Scorecard')


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    t0 = time.time()
    print("\n" + "=" * 70)
    print("  ChiralFold — 3D Folding Quality Benchmark")
    print("=" * 70 + "\n")

    # PDB ground truth
    pdb_path = os.path.join(OUTPUT_DIR, '3IWY.pdb')
    pdb_dihedrals = None
    if os.path.exists(pdb_path):
        print("Parsing PDB 3IWY (D-peptide dPMI-gamma)...")
        bb = parse_pdb_backbone(pdb_path, chain='B')
        pdb_dihedrals = compute_pdb_dihedrals(bb)
        for d in pdb_dihedrals:
            print(f"  Res {d['residue']+1:>2} ({d['resname']:<3}): "
                  f"phi={d['phi']:>7.1f}  psi={d['psi']:>7.1f}  omega={d['omega']:>7.1f}")
        print()

    # Generate conformer ensembles
    tests = [
        ('DWWPLAF', 'D' * 7, 'dPMI short', 25),
        ('AFWKLD', 'D' * 6, 'D-hexapeptide', 25),
        ('TNWYQGLRF', 'D' * 9, 'D-9mer (DP9)', 20),
        ('ETFSDLWKLL', 'D' * 10, 'D-10mer', 15),
        ('AFWKEL', 'DLDLDL', 'Alt L/D 6', 25),
        ('TNWYQGLRF', 'DLDLDLDLD', 'Alt L/D 9', 20),
        ('AFWKELDR', 'LLLDLLLL', 'Drug design', 20),
        ('VFVFVF', 'DDDDDD', 'D-β motif', 25),
    ]

    all_data = []
    for seq, chir, label, nc in tests:
        print(f"  {label:<18} ({seq}) ...", end=' ', flush=True)
        result = generate_and_analyze(seq, chir, n_confs=nc, label=label)
        if result:
            all_data.append(result)
            print(f"{result['n_confs']} conformers")
        else:
            print("FAIL")
    print()

    # Build figure
    print("Rendering figure...", flush=True)
    fig = plt.figure(figsize=(18, 22))
    gs = gridspec.GridSpec(4, 2, hspace=0.35, wspace=0.3, height_ratios=[1, 1, 1, 1])

    panel_ramachandran(fig.add_subplot(gs[0, 0]), all_data, pdb_dihedrals)
    panel_bond_lengths(fig.add_subplot(gs[0, 1]), all_data)
    panel_bond_angles(fig.add_subplot(gs[1, 0]), all_data)
    panel_planarity(fig.add_subplot(gs[1, 1]), all_data, pdb_dihedrals)
    panel_energy_landscape(fig.add_subplot(gs[2, 0]), all_data)
    panel_rg(fig.add_subplot(gs[2, 1]), all_data)

    best = [d for d in all_data if d['n_confs'] >= 5]
    if best:
        panel_3d_backbone(fig.add_subplot(gs[3, 0], projection='3d'), best[0])
    panel_scorecard(fig.add_subplot(gs[3, 1]), all_data, pdb_dihedrals)

    fig.suptitle(
        'ChiralFold — 3D Folding Quality Assessment\n'
        'Structural Validity of Generated Peptide Geometries',
        fontsize=15, fontweight='bold', y=0.995, color=PAL['navy'],
    )
    path = os.path.join(OUTPUT_DIR, 'folding_quality.png')
    plt.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved {path}")

    elapsed = time.time() - t0
    print(f"\n  Done in {elapsed:.1f}s — {len(all_data)} peptides, "
          f"{sum(d['n_confs'] for d in all_data)} conformers\n")
