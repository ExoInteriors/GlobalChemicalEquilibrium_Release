import os
import shutil
from src.constants import repo_root


def copy_inputs(version='Sulfur_Version', input_dir=None):
    input_dir = input_dir or os.path.join(repo_root, f'input_Folder_{version}')
    subfolders = [f.path for f in os.scandir(input_dir) if f.is_dir()]

    # Explicit sulfur-version sources
    solver_src = os.path.join(repo_root, version, 'solver')
    param_src = os.path.join(repo_root, version, 'param.dat')
    Gibbs_src = os.path.join(repo_root, version, 'Gibbs.dat')

    with open(param_src, encoding="utf-8") as pf:
        param_lines = pf.readlines()
    gibbs_root_script = os.path.join(repo_root, 'Gibbs.py') 
    gibbs_root_exists = os.path.isfile(gibbs_root_script)

    for s in sorted(subfolders):
        print(s)
        shutil.copy(solver_src, os.path.join(s, 'solver'))
        target_param = os.path.join(s, 'param.dat')
        with open(target_param, 'w', encoding="utf-8") as out:
            for line in param_lines:
                if gibbs_root_exists and line.strip().startswith("Gibbs energy file =") and line.strip().endswith(".py"):
                    rel_path = os.path.relpath(gibbs_root_script, start=s)
                    out.write(f"Gibbs energy file = {rel_path}\n")
                else:
                    out.write(line)
        shutil.copyfile(Gibbs_src, os.path.join(s, 'Gibbs.dat'))


if __name__ == "__main__":
	copy_inputs()