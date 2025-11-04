# MD trajectory analysis of hydrogen bonds
import os
import math
import pickle
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt

# ================= Constants =================
DIST_CUTOFF_NM = 0.35  # nm
ANGLE_CUTOFF_DEG = 130
ANGLE_CUTOFF_RAD = ANGLE_CUTOFF_DEG * (math.pi / 180)


# ================= Utility Functions =================
def select_atom(topology, atom_info):
    """Select atom index from topology by chain, residue number, and atom name."""
    # atom_info: [[chain_id, res_seq, atom_name], ...] in shape of [N, 3]
    # return: list of mdtraj.Atom objects
    atoms = []
    for chain_id, res_seq, atom_name in atom_info:
        sel = topology.select(f"chainid {chain_id} and resSeq {res_seq} and name {atom_name}")
        if len(sel) == 0:
            #raise ValueError(f"Atom not found: chain {chain_id}, residue {res_seq}, atom {atom_name}")
            atoms.append(None)
        else:
            atoms.append(topology.atom(int(sel[0])))
    return atoms


def plot_histogram(data, cutoff, xlabel, ylabel, title, save_path, cutoff_label):
    """Plot and save histogram with cutoff line."""
    plt.figure(figsize=(5, 3))
    plt.hist(data, bins=50, color="steelblue", alpha=0.8)
    plt.axvline(x=cutoff, color='r', linestyle='--', label=cutoff_label)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()


# ================= Analysis Function =================
def analyze_hydrogen_bonds(
    id: str = "01",
    res_acceptor: list = [1, 47, "OE2"],  # chain 1, residue 47, atom OE2
    res_donor: list = [0, 195, "NZ"],     # chain 0, residue 195, atom NZ
    res_hydrogens: list = [[0, 195, "HZ1"], [0, 195, "HZ2"], [0, 195, "HZ3"]],
    save_path: str = "."
):
    """
    Analyze hydrogen bonds between specific donor/acceptor atoms in MD trajectory.
    Plots histograms for distance and angle distributions.
    """
    # Input trajectory file
    input_pdb = f"/home/lwang/models/gromacs-docker/TrajMPNN/AF_hsa{id}_relaxed_md/rep1/data2/traj.pdb"
    if not os.path.isfile(input_pdb):
        raise FileNotFoundError(f"Input file not found: {input_pdb}")
    traj = md.load(input_pdb)
    print(f"Loaded trajectory with {len(traj)} frames from {input_pdb}")

    # Remove equilibration frames
    traj = traj[10:]
    top = traj.topology

    # Atom selection
    donor = select_atom(top, res_donor)
    acceptor = select_atom(top, res_acceptor)
    hydrogens = [select_atom(top,atom_info) for atom_info in res_hydrogens]
    print(f"Selected donor atoms: {[d.index for d in donor]}")
    print(f"Selected acceptor atoms: {[a.index for a in acceptor]}")
    print(f"Selected hydrogen atoms: {[[h.index for h in h_list] for h_list in hydrogens]}")

    assert len(donor) == len(acceptor) == len(hydrogens), "Mismatch in number of donor, acceptor, or hydrogen pairs."
    dist_pairs, angle_pairs, angle_batch = [], [], []
    for i, (d,h,a) in enumerate(zip(donor, hydrogens, acceptor)):
        if not d is None and not a is None:
            dist_pairs.append([d.index, a.index])
            for j in range(len(h)):
                angle_pairs.append([d.index, h[j].index, a.index])
            angle_batch.extend([i]*len(h))
    dist_pairs = np.array(dist_pairs)
    angle_pairs = np.array(angle_pairs)
    angle_batch = np.array(angle_batch)

    # Compute distances and angles
    hb_dist = md.compute_distances(traj, dist_pairs, periodic=False)
    hb_dist = np.array(hb_dist).T
    hb_angles = md.compute_angles(traj, angle_pairs, periodic=False)
    print(f"Computed distances shape: {hb_dist.shape}, angles shape: {hb_angles.shape}")
    hb_angle = np.array([np.max(hb_angles[:,angle_batch == i], axis=1) for i in np.unique(angle_batch)])
    print(f"Computed distances shape: {hb_dist.shape}, angles shape: {hb_angle.shape}")

    # Apply hydrogen bond criteria, mask shape: [#hbonds, #time steps]
    dist_mask = hb_dist < DIST_CUTOFF_NM
    angle_mask = hb_angle > ANGLE_CUTOFF_RAD
    assert dist_mask.shape == angle_mask.shape
    combined_mask = dist_mask & angle_mask
    hbond_count = np.sum(combined_mask, axis=1)
    hbond_series = combined_mask.astype(int)

    print(f"Frames analyzed: {len(traj)}")
    print(f"H-bond count (dist < {DIST_CUTOFF_NM} nm & angle > {ANGLE_CUTOFF_DEG}°): {hbond_count}")

    # Ensure save directory exists
    os.makedirs(save_path, exist_ok=True)
    data_to_save = {
    "query_atoms": (donor, hydrogens, acceptor),
    "hbond_count": hbond_count,
    "hbond_series": hbond_series,
    "hb_dist": hb_dist,
    "hb_angle": hb_angle,
    "dist_pairs": dist_pairs,
    "angle_pairs": angle_pairs,
    "angle_batch": angle_batch,
    }
    with open(f"{save_path}/{id}_hbond_metadata", "wb") as f:
        pickle.dump(data_to_save, f)

    # Plot distributions
    '''plot_histogram(
        hb_dist,
        DIST_CUTOFF_NM,
        xlabel="Distance (nm)",
        ylabel="Frequency",
        title=f"H-bond Distance ({donor}-{acceptor})",
        save_path=f"{save_path}/{res_acceptor[1]}+{res_donor[1]}_dist_{id}.png",
        cutoff_label=f"Cutoff ({DIST_CUTOFF_NM} nm)",
    )
    plot_histogram(
        hb_angle,
        ANGLE_CUTOFF_RAD,
        xlabel="Angle (radians)",
        ylabel="Frequency",
        title=f"H-bond Angle ({donor}-{acceptor})",
        save_path=f"{save_path}/{res_acceptor[1]}+{res_donor[1]}_angle_{id}.png",
        cutoff_label=f"Cutoff ({ANGLE_CUTOFF_DEG}°)",
    )'''

    return hbond_count, hbond_series


# ================= Main =================
if __name__ == '__main__':
    hb_count = []
    hb_series_all = []
    res_acceptor = [
        [1, 47, "OE2"],
        [0, 317, "OE2"], #Glu305
        [0, 317, "OE1"], #Glu305
        [0, 317, "OE2"], #Glu305
        [0, 317, "OE1"], #Glu305
        [0, 314, "OD1"], #ASN302
        #[0, 314, "O"], #ASN302/ absent in hsa01
        [0, 226, "OE1"], #214
        [0, 226, "OE2"], #214
        [0, 263, "OD1"], #251
        ]
    res_donor = [
        [0, 208, "NZ"],
        [1, 24, "N"], 
        [1, 24, "N"], 
        [1, 25, "N"], 
        [1, 25, "N"], 
        [1, 25, "OG"], 
        #[1, 28, "OH"],
        [1, 37, "N"],
        [1, 37, "N"],
        [1, 37, "OG1"],
        ]   
    res_hydrogens = [
        [[0, 208, "HZ1"], [0, 208, "HZ2"], [0, 208, "HZ3"]],
        [[1, 24, "H"]],
        [[1, 24, "H"]],
        [[1, 25, "H"]],
        [[1, 25, "H"]],
        [[1, 25, "HG"]],
        #[[1, 28, "HH"]],
        [[1, 37, "H"]],
        [[1, 37, "H"]],
        [[1, 37, "HG1"]],
    ]

    for i in range(1, 11):
        try:
            print(f"Analyzing system {i:02d}...")
            count, series = analyze_hydrogen_bonds(
                    id=f"{i:02d}",
                    res_acceptor=res_acceptor,
                    res_donor=res_donor,
                    res_hydrogens=res_hydrogens,
                    save_path="/home/lwang/models/gromacs-docker/TrajMPNN/1tf0_fixed_md/binder_hb_fig"
                    )
            hb_count.append(count)
            hb_series_all.append(series)
        except Exception as e:
            hb_count.append(None)
            hb_series_all.append(None)
            continue
    print("Hydrogen bond counts for all systems:", hb_count)
