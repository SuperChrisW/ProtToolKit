#!/usr/bin/env python3
"""
Computational alanine scanning of an AF3-relaxed complex using PyRosetta.

Outputs CSV with:
resnum,chain,resname,side,ddg_mean,ddg_sd,dG_wt,dG_mut_mean,n_trials

Usage:
    python run_ddG_scan.py \
        --pdb complex_relaxed.pdb \
        --partners AB_C \
        --ntrials 3 \
        --shell 8.0 \
        --out ddg_scan.csv
"""

import argparse, csv, random, math
from statistics import mean, pstdev

import pyrosetta
from pyrosetta import rosetta, toolbox as Rtoolbox

def init_rosetta():
    pyrosetta.init(
        "-mute all "
        "-detect_disulf true "
        "-include_current true "
        "-ex1 -ex2aro "
        "-use_input_sc "
        "-ignore_unrecognized_res false "
        "-ignore_zero_occupancy false "
        "-linmem_ig 10 "
        "-relax:cartesian "
        "-score:weights ref2015"
    )

def get_partners_sets(pose, partners_str):
    # partners_str like "AB_C" (left side vs right side)
    left, right = partners_str.split("_")
    pdbinfo = pose.pdb_info()
    left_set, right_set = set(), set()
    for i in range(1, pose.size()+1):
        ch = pdbinfo.chain(i)
        if ch in left:
            left_set.add(i)
        if ch in right:
            right_set.add(i)
    return left_set, right_set

def select_interface_residues(pose, left_set, right_set, dist_cut=8.0):
    # simple heavy-atom distance based interface definition
    # mark residues that have any heavy-atom pair within dist_cut across partners
    xyz = pose
    interface = set()
    for i in left_set:
        ri = pose.residue(i)
        if not ri.is_protein(): 
            continue
        for j in right_set:
            rj = pose.residue(j)
            if not rj.is_protein():
                continue
            close = False
            for ai in range(1, ri.natoms()+1):
                if ri.atom_type(ai).is_hydrogen(): 
                    continue
                xi = ri.atom(ai).xyz()
                for aj in range(1, rj.natoms()+1):
                    if rj.atom_type(aj).is_hydrogen():
                        continue
                    if xi.distance(rj.atom(aj).xyz()) <= dist_cut:
                        close = True
                        break
                if close: break
            if close:
                #interface.add(i) # skip add receptor chain residues here
                interface.add(j)
    return sorted(interface)

def mutate_to_ala_with_repack(pose, resi, pack_radius=8.0, scorefxn=None, min_bb=False):
    """Clone pose, mutate position to Ala, repack neighbors, optional minimization. Return mutated pose."""
    if pose.residue(resi).name1() == 'A':
        return pose.clone()  # already Ala

    if scorefxn is None:
        scorefxn = rosetta.core.scoring.get_score_function()  # ref2015 default

    mut_pose = pose.clone()
    # Mutate & repack neighbors within pack_radius
    # mutant_aa: single-letter code
    Rtoolbox.mutants.mutate_residue(mut_pose, resi, 'A', pack_radius=pack_radius, pack_scorefxn=scorefxn)

    # Optional: local minimization around the mutated site (chi, and backbone if requested)
    mm = rosetta.core.kinematics.MoveMap()
    # Build a neighborhood mask by distance to mutated residue CA
    res_xyz = mut_pose.residue(resi).xyz("CA") if mut_pose.residue(resi).has("CA") else mut_pose.residue(resi).nbr_atom_xyz()
    for i in range(1, mut_pose.size()+1):
        if not mut_pose.residue(i).is_protein(): 
            continue
        # distance-based shell
        dmin = min(res_xyz.distance(mut_pose.residue(i).xyz(a)) for a in range(1, mut_pose.residue(i).natoms()+1))
        if dmin <= pack_radius:
            mm.set_chi(i, True)
            if min_bb:
                mm.set_bb(i, True)
    mm.set_jump(True)

    min_mover = rosetta.protocols.minimization_packing.MinMover()
    min_mover.score_function(scorefxn)
    min_mover.movemap(mm)
    min_mover.min_type("lbfgs_armijo_nonmonotone")
    min_mover.tolerance(1e-3)
    min_mover.apply(mut_pose)

    return mut_pose

def compute_dG_binding(pose, partners):
    iam = rosetta.protocols.analysis.InterfaceAnalyzerMover(
        partners, False, rosetta.core.scoring.ScoreFunctionFactory.create_score_function("ref2015")
    )
    iam.set_pack_separated(True)
    iam.set_calc_dSASA(True)
    iam.apply(pose)
    return iam.get_interface_dG()

def alanine_scan_with_toolbox(pose, partners, iface_residues, ntrials=3, pack_radius=8.0, min_bb=False, scorefxn=None):
    """Return list of rows with ΔΔG per residue using mutate_residue()."""
    if scorefxn is None:
        scorefxn = rosetta.core.scoring.ScoreFunctionFactory.create_score_function("ref2015")

    wt_pose = pose.clone()
    dG_wt = compute_dG_binding(wt_pose, partners)

    results = []
    pdbinfo = pose.pdb_info()
    for resi in iface_residues:
        r = pose.residue(resi)
        if not r.is_protein():
            continue
        if r.name1() == 'A':
            continue  # no-op

        trials = []
        for _ in range(ntrials):
            mut_pose = mutate_to_ala_with_repack(pose, resi, pack_radius=pack_radius, scorefxn=scorefxn, min_bb=min_bb)
            dG_mut = compute_dG_binding(mut_pose, partners)
            trials.append(dG_mut)

        ddg_vals = [d - dG_wt for d in trials]
        results.append({
            "resnum": pdbinfo.number(resi),
            "chain": pdbinfo.chain(resi),
            "resname": r.name3(),
            "side": "binder_or_receptor",  # fill from your left/right sets if you want
            "ddg_mean": sum(ddg_vals)/len(ddg_vals),
            "ddg_sd": (0.0 if len(ddg_vals)==1 else (sum((x - (sum(ddg_vals)/len(ddg_vals)))**2 for x in ddg_vals)/len(ddg_vals))**0.5),
            "dG_wt": dG_wt,
            "dG_mut_mean": sum(trials)/len(trials),
            "n_trials": ntrials
        })
    return results

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--partners", required=True, help='Partner chains like "AB_C" (receptor_binder)')
    ap.add_argument("--ntrials", type=int, default=3)
    ap.add_argument("--shell", type=float, default=8.0)
    ap.add_argument("--out", default="ddg_scan.csv")
    ap.add_argument("--min_bb", action="store_true", help="Allow backbone minimization in shell")
    args = ap.parse_args()

    init_rosetta()
    pose = rosetta.core.import_pose.pose_from_file(args.pdb)
    left, right = get_partners_sets(pose, args.partners)
    iface_residues = select_interface_residues(pose, left, right, dist_cut=args.shell)

    results = alanine_scan_with_toolbox(
        pose,
        partners=args.partners,
        iface_residues=iface_residues,
        ntrials=args.ntrials,
        pack_radius=args.shell,
        min_bb=args.min_bb
    )

    # Sort: strongest hotspots first (more positive ΔΔG indicates loss of binding upon Ala)
    results.sort(key=lambda x: (-x["ddg_mean"], x["side"], x["chain"], x["resnum"]))

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["resnum","chain","resname","side","ddg_mean","ddg_sd","dG_wt","dG_mut_mean","n_trials"])
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    # Quick console summary
    print(f"WT dG_bind (InterfaceAnalyzer dG_separated): {results[0]['dG_wt']:.3f} REU")
    print("Top 10 hotspots (ΔΔG_bind, REU):")
    for r in results[:10]:
        print(f"{r['chain']}{r['resnum']:>4} {r['resname']:>3} [{r['side'][:1]}]  +{r['ddg_mean']:.2f} ±{r['ddg_sd']:.2f}")

if __name__ == "__main__":
    main()
