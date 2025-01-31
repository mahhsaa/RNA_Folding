import os
import math
from collections import defaultdict

base_pairs = ['AA', 'AU', 'AC', 'AG', 'UU', 'UC', 'UG', 'CC', 'CG', 'GG']

def parse_pdb(pdb_path):
    c3_atoms = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                atom_name = line[12:16].strip()
                if atom_name == "C3'":
                    res_name = line[17:20].strip()
                    chain = line[21]
                    res_id = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    c3_atoms.append({
                        'chain': chain,
                        'res_id': res_id,
                        'res_name': res_name,
                        'coords': (x, y, z)
                    })
    return c3_atoms

def process_pdb(pdb_path, observed, reference):
    c3_atoms = parse_pdb(pdb_path)
    chains = defaultdict(list)
    for atom in c3_atoms:
        chains[atom['chain']].append(atom)
    for chain_id, atoms in chains.items():
        sorted_atoms = sorted(atoms, key=lambda x: x['res_id'])
        for i in range(len(sorted_atoms)):
            for j in range(i + 4, len(sorted_atoms)):
                atom_i = sorted_atoms[i]
                atom_j = sorted_atoms[j]
                dx = atom_i['coords'][0] - atom_j['coords'][0]
                dy = atom_i['coords'][1] - atom_j['coords'][1]
                dz = atom_i['coords'][2] - atom_j['coords'][2]
                distance = math.sqrt(dx**2 + dy**2 + dz**2)
                if distance >= 20:
                    continue
                bin_idx = int(distance)
                pair = ''.join(sorted([atom_i['res_name'], atom_j['res_name']]))
                if pair in base_pairs:
                    observed[pair][bin_idx] += 1
                reference[bin_idx] += 1

def main(pdb_dir, output_dir):
    observed = {pair: [0]*20 for pair in base_pairs}
    reference = [0]*20
    for filename in os.listdir(pdb_dir):
        if filename.endswith('.pdb'):
            pdb_path = os.path.join(pdb_dir, filename)
            process_pdb(pdb_path, observed, reference)
    total_ref = sum(reference)
    ref_freq = [count / total_ref if total_ref != 0 else 0 for count in reference]
    for pair in base_pairs:
        total_obs = sum(observed[pair])
        scores = []
        for bin_idx in range(20):
            if total_obs == 0 or observed[pair][bin_idx] == 0:
                scores.append(10.0)
            else:
                obs_freq = observed[pair][bin_idx] / total_obs
                ratio = obs_freq / (ref_freq[bin_idx] if ref_freq[bin_idx] != 0 else 1e-9)
                score = -math.log(ratio)
                scores.append(min(score, 10.0))
        with open(os.path.join(output_dir, f'{pair}.txt'), 'w') as f:
            for score in scores:
                f.write(f"{score}\n")

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_dir', help='Directory containing PDB files')
    parser.add_argument('output_dir', help='Directory to save score files')
    args = parser.parse_args()
    main(args.pdb_dir, args.output_dir)
