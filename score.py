import os
import math

def parse_pdb(pdb_path):
    c3_atoms = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                atom_name = line[12:16].strip()
                if atom_name == "C3'":
                    res_name = line[17:20].strip()
                    res_id = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    c3_atoms.append({
                        'res_id': res_id,
                        'res_name': res_name,
                        'coords': (x, y, z)
                    })
    return c3_atoms

def load_scores(score_dir):
    scores = {}
    for pair in ['AA', 'AU', 'AC', 'AG', 'UU', 'UC', 'UG', 'CC', 'CG', 'GG']:
        with open(os.path.join(score_dir, f'{pair}.txt'), 'r') as f:
            scores[pair] = [float(line.strip()) for line in f]
    return scores

def compute_energy(pdb_path, scores):
    c3_atoms = parse_pdb(pdb_path)
    c3_atoms.sort(key=lambda x: x['res_id'])
    total = 0.0
    for i in range(len(c3_atoms)):
        for j in range(i + 4, len(c3_atoms)):
            atom_i = c3_atoms[i]
            atom_j = c3_atoms[j]
            dx = atom_i['coords'][0] - atom_j['coords'][0]
            dy = atom_i['coords'][1] - atom_j['coords'][1]
            dz = atom_i['coords'][2] - atom_j['coords'][2]
            distance = math.sqrt(dx**2 + dy**2 + dz**2)
            if distance >= 20:
                continue
            bin_idx = int(distance)
            frac = distance - bin_idx
            pair = ''.join(sorted([atom_i['res_name'], atom_j['res_name']]))
            if pair not in scores:
                continue
            if bin_idx >= 19:
                score = scores[pair][19]
            else:
                score = (1 - frac) * scores[pair][bin_idx] + frac * scores[pair][bin_idx + 1]
            total += score
    return total

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_file', help='PDB file to score')
    parser.add_argument('score_dir', help='Directory with score files')
    args = parser.parse_args()
    scores = load_scores(args.score_dir)
    energy = compute_energy(args.pdb_file, scores)
    print(f"Estimated Gibbs Free Energy: {energy:.2f}")
