"""Regenerate the precomputed trajectories in data/trajectories/ (and the erlotinib parameters in data/ligands/).

The analysis notebooks (09, 10, 11, 12) read these files so that they run in seconds on Colab and in CI.
Everything here uses the same force fields and settings that the notebooks teach. Run it inside the
conda environment from environment.yml (it needs AmberTools for the ligand charges):

    python scripts/generate_trajectories.py all            # everything (about 1 hour on an 8-core CPU)
    python scripts/generate_trajectories.py alanine        # alanine dipeptide, implicit solvent (notebook 10)
    python scripts/generate_trajectories.py chignolin      # chignolin CLN025 in water, 3 replicas (notebooks 08, 09, 12)

Lengths are short on purpose (teaching data that can be rebuilt on a laptop). Make them longer with
CHIGNOLIN_NS (default 0.25), CHIGNOLIN_EQ_PS (25), EGFR_PS (60), EGFR_EQ_PS (20). Set CHIGNOLIN_REPLICA=i to
run one replica (e.g. in parallel), then `python scripts/generate_trajectories.py merge_chignolin`.
    python scripts/generate_trajectories.py ligand         # erlotinib GAFF2/AM1-BCC parameters (notebook 11)
    python scripts/generate_trajectories.py egfr           # EGFR kinase-erlotinib complex in water (notebook 11)

A GPU makes this much faster: set OPENMM_PLATFORM=CUDA (or OpenCL).
"""
import datetime as dt
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import openmm
from openmm import app, unit

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
STRUCT = DATA / "structures"
TRAJ = DATA / "trajectories"
LIG = DATA / "ligands"
PLATFORM = os.environ.get("OPENMM_PLATFORM", "CPU")
ERLOTINIB_SMILES = "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC"


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def provenance(key, value):
    f = DATA / "provenance.json"
    meta = json.loads(f.read_text()) if f.exists() else {}
    value = dict(value, generated=dt.date.today().isoformat(), openmm=openmm.__version__)
    meta[key] = value
    f.write_text(json.dumps(meta, indent=2) + "\n")


def make_sim(topology, system, integrator, seed=None):
    platform = openmm.Platform.getPlatformByName(PLATFORM)
    if seed is not None and hasattr(integrator, "setRandomNumberSeed"):
        integrator.setRandomNumberSeed(seed)
    return app.Simulation(topology, system, integrator, platform)


def add_restraints(system, positions, atoms, k=400.0):
    """Harmonic position restraints (kJ/mol/nm^2) on the given atom indices, periodic-safe."""
    force = openmm.CustomExternalForce("0.5*k*periodicdistance(x, y, z, x0, y0, z0)^2")
    force.addGlobalParameter("k", k)
    for p in ("x0", "y0", "z0"):
        force.addPerParticleParameter(p)
    for i in atoms:
        force.addParticle(i, positions[i].value_in_unit(unit.nanometer))
    return system.addForce(force)


# --------------------------------------------------------------------------- alanine dipeptide
def alanine(n_rep=4, ns=10.0, save_ps=1.0):
    import mdtraj as md
    pdb = app.PDBFile(str(STRUCT / "alanine_dipeptide.pdb"))
    ff = app.ForceField("amber14-all.xml", "implicit/obc2.xml")
    dt_ps, rows = 0.002, []
    for rep in range(n_rep):
        system = ff.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
        integ = openmm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, dt_ps * unit.picoseconds)
        sim = make_sim(pdb.topology, system, integ, seed=1000 + rep)
        sim.context.setPositions(pdb.positions)
        sim.minimizeEnergy()
        sim.context.setVelocitiesToTemperature(300 * unit.kelvin, 2000 + rep)
        sim.step(int(100 / dt_ps))  # 100 ps equilibration (discarded)
        stride = int(save_ps / dt_ps)
        n = int(ns * 1000 / save_ps)
        xyz = np.zeros((n, pdb.topology.getNumAtoms(), 3), dtype=np.float32)
        t0 = time.time()
        for f in range(n):
            sim.step(stride)
            xyz[f] = sim.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        log(f"alanine replica {rep}: {ns} ns in {time.time() - t0:.0f} s")
        traj = md.Trajectory(xyz, md.Topology.from_openmm(pdb.topology))
        _, phi = md.compute_phi(traj)
        _, psi = md.compute_psi(traj)
        for f in range(n):
            rows.append((rep, round((f + 1) * save_ps, 3), round(float(np.degrees(phi[f, 0])), 2),
                         round(float(np.degrees(psi[f, 0])), 2)))
    TRAJ.mkdir(parents=True, exist_ok=True)
    with open(TRAJ / "alanine_dipeptide_phi_psi.csv", "w") as fh:
        fh.write("replica,time_ps,phi_deg,psi_deg\n")
        fh.writelines(f"{r},{t},{p},{s}\n" for r, t, p, s in rows)
    provenance("alanine_dipeptide_phi_psi.csv", {
        "system": "ACE-ALA-NME built with AmberTools tleap (scripts/build_alanine_dipeptide.leap)",
        "force_field": "amber14-all.xml (ff14SB) + implicit/obc2.xml (GB-OBC II)",
        "ensemble": "NVT, LangevinMiddleIntegrator 300 K, friction 1/ps, 2 fs, HBonds constraints, no cutoff",
        "replicas": n_rep, "ns_per_replica": ns, "save_interval_ps": save_ps,
        "equilibration_ps_discarded": 100, "seeds": "integrator 1000+rep, velocities 2000+rep"})


# --------------------------------------------------------------------------- chignolin
def prepare_chignolin():
    from pdbfixer import PDBFixer
    fixer = PDBFixer(str(STRUCT / "5AWL.pdb"))
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    ff = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    model = app.Modeller(fixer.topology, fixer.positions)
    model.addHydrogens(ff, pH=7.0)
    return ff, model


def equilibrate(sim, system, model, restrained, nvt_ps, npt_ps, dt_ps, seed):
    ridx = add_restraints(system, model.positions, restrained)
    sim.context.reinitialize(preserveState=False)
    sim.context.setPositions(model.positions)
    sim.minimizeEnergy(maxIterations=2000)
    sim.context.setVelocitiesToTemperature(300 * unit.kelvin, seed)
    baro = system.addForce(openmm.MonteCarloBarostat(1 * unit.bar, 300 * unit.kelvin))
    system.getForce(baro).setFrequency(0)  # off during NVT
    sim.context.reinitialize(preserveState=True)
    sim.step(int(nvt_ps / dt_ps))
    system.getForce(baro).setFrequency(25)
    sim.context.reinitialize(preserveState=True)
    sim.step(int(npt_ps / dt_ps))
    system.removeForce(ridx)
    sim.context.reinitialize(preserveState=True)


def chignolin(n_rep=3, ns=float(os.environ.get("CHIGNOLIN_NS", 0.25)), save_ps=2.0,
              eq_ps=float(os.environ.get("CHIGNOLIN_EQ_PS", 25.0))):
    import mdtraj as md
    ff, model = prepare_chignolin()
    protein_atoms = [a.index for a in model.topology.atoms()]
    app.PDBFile.writeFile(model.topology, model.positions, open(TRAJ / "chignolin_protein.pdb", "w"))
    model.addSolvent(ff, model="tip3p", padding=1.0 * unit.nanometer, ionicStrength=0.15 * unit.molar)
    log("chignolin solvated:", model.topology.getNumAtoms(), "atoms")
    heavy = [a.index for a in model.topology.atoms() if a.index in set(protein_atoms) and a.element.symbol != "H"]
    dt_ps = 0.004
    rows = []
    only = os.environ.get("CHIGNOLIN_REPLICA")
    for rep in ([int(only)] if only is not None else range(n_rep)):
        system = ff.createSystem(model.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer,
                                 constraints=app.HBonds, hydrogenMass=1.5 * unit.amu)
        integ = openmm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, dt_ps * unit.picoseconds)
        sim = make_sim(model.topology, system, integ, seed=3000 + rep)
        equilibrate(sim, system, model, heavy, nvt_ps=eq_ps, npt_ps=eq_ps, dt_ps=dt_ps, seed=4000 + rep)
        stride, n = int(save_ps / dt_ps), int(ns * 1000 / save_ps)
        xyz = np.zeros((n, len(protein_atoms), 3), dtype=np.float32)
        t0 = time.time()
        for f in range(n):
            sim.step(stride)
            st = sim.context.getState(getPositions=True, getEnergy=True, enforcePeriodicBox=False)
            pos = st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            xyz[f] = pos[protein_atoms]
            box = st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
            ke = st.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
            ndof = 3 * system.getNumParticles() - system.getNumConstraints() - 3
            temp = 2 * ke / (ndof * unit.MOLAR_GAS_CONSTANT_R.value_in_unit(unit.kilojoule_per_mole / unit.kelvin))
            rows.append((rep, round((f + 1) * save_ps, 1),
                         round(st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole), 1),
                         round(temp, 2), round(float(np.linalg.det(box)), 4)))
        log(f"chignolin replica {rep}: {ns} ns in {time.time() - t0:.0f} s")
        top = md.Topology.from_openmm(model.topology).subset(protein_atoms)
        traj = md.Trajectory(xyz, top, time=np.arange(1, n + 1) * save_ps)
        # make the peptide whole and centred: it is a single small molecule in a box much larger than itself
        traj.superpose(traj, 0)
        traj.save_dcd(str(TRAJ / f"chignolin_rep{rep}.dcd"))
    suffix = f"_rep{only}" if only is not None else ""
    with open(TRAJ / f"chignolin_thermo{suffix}.csv", "w") as fh:
        fh.write("replica,time_ps,potential_kJmol,temperature_K,volume_nm3\n")
        fh.writelines(",".join(map(str, r)) + "\n" for r in rows)
    if only is not None:
        return  # merge the per-replica logs with `python scripts/generate_trajectories.py merge_chignolin`
    provenance("chignolin", {
        "system": "CLN025 chignolin variant YYDPETGTWY, PDB 5AWL chain A, prepared with PDBFixer, H added at pH 7",
        "solvent": "TIP3P-FB water, 1.0 nm padding, cubic box, 0.15 M NaCl", "atoms": model.topology.getNumAtoms(),
        "force_field": "amber14-all.xml (ff14SB) + amber14/tip3pfb.xml",
        "protocol": f"minimize; {eq_ps:g} ps NVT + {eq_ps:g} ps NPT with 400 kJ/mol/nm2 heavy-atom restraints; production NPT",
        "integrator": "LangevinMiddleIntegrator 300 K, 1/ps, 4 fs, HBonds constraints, hydrogen mass 1.5 amu (HMR)",
        "barostat": "MonteCarloBarostat 1 bar, every 25 steps", "nonbonded": "PME, 1.0 nm cutoff",
        "replicas": n_rep, "ns_per_replica": ns, "save_interval_ps": save_ps,
        "files": "chignolin_protein.pdb (topology), chignolin_rep*.dcd (protein only, superposed on frame 0), chignolin_thermo.csv"})


def merge_chignolin():
    import glob
    parts = sorted(glob.glob(str(TRAJ / "chignolin_thermo_rep*.csv")))
    lines = ["replica,time_ps,potential_kJmol,temperature_K,volume_nm3\n"]
    for p in parts:
        lines += open(p).readlines()[1:]
        os.remove(p)
    (TRAJ / "chignolin_thermo.csv").write_text("".join(lines))


# --------------------------------------------------------------------------- EGFR-erlotinib
def erlotinib_from_1m17():
    """Erlotinib (AQ4) from PDB 1M17 with bond orders from SMILES and hydrogens added (RDKit)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    lines = [l for l in open(STRUCT / "1M17.pdb") if l.startswith("HETATM") and l[17:20] == "AQ4"]
    lig = Chem.MolFromPDBBlock("".join(lines), removeHs=False)
    lig = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(ERLOTINIB_SMILES), lig)
    return Chem.AddHs(lig, addCoords=True, addResidueInfo=True)


def ligand():
    from openff.toolkit import Molecule
    from openmmforcefields.generators import GAFFTemplateGenerator
    from rdkit import Chem
    mol = erlotinib_from_1m17()
    LIG.mkdir(parents=True, exist_ok=True)
    Chem.MolToMolFile(mol, str(LIG / "erlotinib_1M17_H.sdf"))
    Chem.MolToPDBFile(mol, str(LIG / "erlotinib_1M17_H.pdb"))
    off = Molecule.from_rdkit(mol, allow_undefined_stereo=True)
    off.name = "AQ4"
    t0 = time.time()
    gen = GAFFTemplateGenerator(molecules=off, forcefield="gaff-2.11")
    xml = gen.generate_residue_template(off)
    (LIG / "erlotinib_gaff2.xml").write_text(xml)
    log(f"erlotinib GAFF2/AM1-BCC parameters in {time.time() - t0:.0f} s")
    provenance("erlotinib_gaff2.xml", {
        "molecule": "erlotinib (PDB ligand AQ4), neutral, coordinates from PDB 1M17, H added by RDKit",
        "smiles": ERLOTINIB_SMILES, "force_field": "GAFF 2.11 via openmmforcefields GAFFTemplateGenerator",
        "charges": "AM1-BCC (AmberTools sqm via OpenFF toolkit)"})


def egfr(ps=float(os.environ.get("EGFR_PS", 60.0)), save_ps=1.0, eq_ps=float(os.environ.get("EGFR_EQ_PS", 20.0))):
    from pdbfixer import PDBFixer
    import mdtraj as md
    # Keep residues 672-964 (the whole kinase domain, 688-955 in this numbering). Residues 965-976 are not
    # resolved in the crystal; dropping the short C-terminal tail 977-995 avoids modelling that gap.
    keep = [l for l in open(STRUCT / "1M17.pdb")
            if not (l.startswith(("ATOM", "HETATM", "ANISOU")) and l[17:20] != "HOH" and int(l[22:26]) > 964)]
    tmp = TRAJ / "_1m17_trimmed.pdb"
    tmp.write_text("".join(keep))
    fixer = PDBFixer(str(tmp))
    tmp.unlink()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    fixer.missingResidues = {}  # do not build the disordered N-terminal tag
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    ff = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml", str(LIG / "erlotinib_gaff2.xml"))
    model = app.Modeller(fixer.topology, fixer.positions)
    model.addHydrogens(ff, pH=7.0)
    lig = app.PDBFile(str(LIG / "erlotinib_1M17_H.pdb"))
    model.add(lig.topology, lig.positions)
    solute = [a.index for a in model.topology.atoms()]
    app.PDBFile.writeFile(model.topology, model.positions, open(TRAJ / "egfr_erlotinib_complex.pdb", "w"))
    model.addSolvent(ff, model="tip3p", padding=1.0 * unit.nanometer, ionicStrength=0.15 * unit.molar,
                     boxShape="octahedron")
    log("EGFR-erlotinib solvated:", model.topology.getNumAtoms(), "atoms")
    system = ff.createSystem(model.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer,
                             constraints=app.HBonds, hydrogenMass=1.5 * unit.amu)
    dt_ps = 0.004
    integ = openmm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, dt_ps * unit.picoseconds)
    sim = make_sim(model.topology, system, integ, seed=5000)
    heavy = [i for i, a in enumerate(model.topology.atoms()) if i < len(solute) and a.element.symbol != "H"]
    t0 = time.time()
    equilibrate(sim, system, model, heavy, nvt_ps=eq_ps, npt_ps=eq_ps, dt_ps=dt_ps, seed=5001)
    log(f"EGFR equilibration in {time.time() - t0:.0f} s")
    stride, n = int(save_ps / dt_ps), int(ps / save_ps)
    xyz = np.zeros((n, len(solute), 3), dtype=np.float32)
    t0 = time.time()
    for f in range(n):
        sim.step(stride)
        xyz[f] = sim.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)[:len(solute)]
        if f % 10 == 9:
            log(f"  EGFR {(f + 1) * save_ps:.0f} ps")
    log(f"EGFR production {ps} ps in {time.time() - t0:.0f} s")
    top = md.Topology.from_openmm(model.topology).subset(solute)
    traj = md.Trajectory(xyz, top, time=np.arange(1, n + 1) * save_ps)

    traj.superpose(traj, 0, atom_indices=top.select("protein and name CA"))
    traj.save_dcd(str(TRAJ / "egfr_erlotinib.dcd"))
    provenance("egfr_erlotinib", {
        "system": "EGFR kinase domain (PDB 1M17 chain A residues 672-964, unresolved N-terminal tag not built) + erlotinib (AQ4)",
        "solvent": "TIP3P-FB water, 1.0 nm padding, truncated octahedron, 0.15 M NaCl", "atoms": model.topology.getNumAtoms(),
        "force_field": "amber14-all.xml (ff14SB) + amber14/tip3pfb.xml; erlotinib GAFF 2.11 / AM1-BCC (data/ligands/erlotinib_gaff2.xml)",
        "protocol": f"minimize; {eq_ps:g} ps NVT + {eq_ps:g} ps NPT with 400 kJ/mol/nm2 restraints on solute heavy atoms; production NPT",
        "integrator": "LangevinMiddleIntegrator 300 K, 1/ps, 4 fs, HBonds constraints, hydrogen mass 1.5 amu (HMR)",
        "nonbonded": "PME, 1.0 nm cutoff", "production_ps": ps, "save_interval_ps": save_ps,
        "files": "egfr_erlotinib_complex.pdb (topology, solute only), egfr_erlotinib.dcd (solute only, superposed on protein CA)"})


if __name__ == "__main__":
    TRAJ.mkdir(parents=True, exist_ok=True)
    steps = sys.argv[1:] or ["all"]
    if "all" in steps:
        steps = ["alanine", "chignolin", "ligand", "egfr"]
    for s in steps:
        log("==", s)
        globals()[s]()
