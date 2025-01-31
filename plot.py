import os

base_pairs = ['AA', 'AU', 'AC', 'AG', 'UU', 'UC', 'UG', 'CC', 'CG', 'GG']

def plot_profiles(score_dir):
    for pair in base_pairs:
        with open(os.path.join(score_dir, f'{pair}.txt'), 'r') as f:
            scores = [float(line.strip()) for line in f]
        print(f"\n=== Scoring Profile: {pair} ===")
        for bin in range(20):
            print(f"Bin {bin:2d} ({bin}-{bin+1} Å): {'*' * int(scores[bin])} ({scores[bin]:.2f})")

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('score_dir', help='Directory with score files')
    args = parser.parse_args()
    plot_profiles(args.score_dir)
