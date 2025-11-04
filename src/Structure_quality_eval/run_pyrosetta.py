import os
import sys
import contextlib
import pandas as pd
from pyrosetta import init, pose_from_pdb, get_fa_scorefxn
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

# --- Silence PyRosetta output ---
@contextlib.contextmanager
def suppress_stdout_stderr():
    # Redirect stdout and stderr at the OS level (not Python level)
    with open(os.devnull, 'w') as devnull:
        old_stdout_fd = os.dup(1)
        old_stderr_fd = os.dup(2)

        os.dup2(devnull.fileno(), 1)
        os.dup2(devnull.fileno(), 2)
        try:
            yield
        finally:
            os.dup2(old_stdout_fd, 1)
            os.dup2(old_stderr_fd, 2)

# --- Initial binding energy before relaxation ---
def compute_binding_energy(pose: Pose, chains="A_B", scorefxn=None):
    ia = InterfaceAnalyzerMover(chains, False, scorefxn)  # <- FIXED here
    with suppress_stdout_stderr():
        ia.apply(pose)
    return ia.get_interface_dG()

def run_pyr(pdb_path, chains="A_B"):
    # --- Initialization ---
    init("-relax:fast")
    # --- Load complex ---
    pose = pose_from_pdb(pdb_path)

    # --- Score function ---
    scorefxn = get_fa_scorefxn()

    initial_dG = compute_binding_energy(pose, chains, scorefxn)
    print(f"Initial ΔG_binding: {initial_dG:.2f} REU")

    # --- Relax the complex ---
    relax = FastRelax()
    relax.set_scorefxn(scorefxn)
    relax.apply(pose)

    # --- Binding energy after relaxation ---
    relaxed_dG = compute_binding_energy(pose)
    print(f"Relaxed ΔG_binding: {relaxed_dG:.2f} REU")

    # --- Save the relaxed structure ---
    pose.dump_pdb(pdb_path.replace('.pdb', '_relaxed.pdb'))

    return initial_dG, relaxed_dG

if __name__ == "__main__":
    ### it will take ± 17 mins for one EGFR model
    work_dir = "/home/lwang/models/Prot_data_prepare/dataset/Andrejs_protMPNN_AF"
    model_list = [file for file in os.listdir(work_dir) if file.endswith('.pdb') and (file.startswith('fold_hsa_') or file.startswith('fold_egfr_'))]
    
    custom_model_list = [
        'fold_egfr_bindermut02_2_model_0',
        'fold_egfr_bindermut04_2_model_2',
        'fold_egfr_bindermut06_1_model_0',
        'fold_egfr_bindermut09_2_model_1',
        'fold_egfr_bindermut10_1_model_0',
        'fold_hsa_gamut03_2_model_0',
        'fold_hsa_gamut04_1_model_0',
        'fold_hsa_gamut05_1_model_4',
        'fold_hsa_gamut08_1_model_3',
        'fold_hsa_gamut10_1_model_3',
    ]
    model_list = [f"{model}.pdb" for model in custom_model_list]
    print(f"Found {len(model_list)} models to process.")
    dG_results = {}
    for i, model in enumerate(model_list):
        if os.path.exists(os.path.join(work_dir, model.replace('.pdb', '_relaxed.pdb'))):
            print(f"{model} already relaxed, skip.")
            continue
        print(f"Processing {model}...")
        model_fpath = os.path.join(work_dir, model)
        init_dG, relax_dG = run_pyr(model_fpath)
        dG_results[model] = [init_dG, relax_dG]
        print(model, init_dG, relax_dG)

    #dG_df = pd.DataFrame(dG_results).T
    #dG_df.columns = ['init_ddG', 'relax_ddG']
    #dG_df.to_csv(os.path.join(work_dir, 'dG_results.csv'), index=True)