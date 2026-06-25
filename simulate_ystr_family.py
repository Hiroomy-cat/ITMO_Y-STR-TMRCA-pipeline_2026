"""
Y-STR Family Simulation
=======================
Simulates a Y-STR haplotype tree over N_GENERATIONS generations,
each individual producing SONS_PER_PARENT sons.

Outputs:
  simulation_output/genotypes/   — one TSV file per individual
  simulation_output/family_tree.png  — genealogical tree image
  simulation_output/mrca_pairs.csv   — MRCA distances for all pairs in last gen
"""

import os
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from itertools import combinations

# ── Configuration ────────────────────────────────────────────────────────────
RANDOM_SEED      = 42
BED_FILE         = "output_hg19_full_with_rates.bed"
OUTPUT_DIR       = "simulation_output"
GENO_DIR         = os.path.join(OUTPUT_DIR, "genotypes")
N_GENERATIONS    = 8
SONS_PER_PARENT  = 2
ANCESTOR_REPEAT  = 10          # starting repeat count for every locus
MIN_REPEATS      = 1           # floor for repeat count
# ─────────────────────────────────────────────────────────────────────────────

np.random.seed(RANDOM_SEED)

os.makedirs(GENO_DIR, exist_ok=True)


# ── 1. Load loci with known mutation rates ────────────────────────────────────
loci = []
with open(BED_FILE) as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 7:
            continue
        rate_str = parts[6]
        if rate_str == "NA":
            continue
        loci.append({"name": parts[5], "rate": float(rate_str)})

# Collapse duplicate locus names (same name, same rate → keep one entry)
seen = {}
for loc in loci:
    if loc["name"] not in seen:
        seen[loc["name"]] = loc
loci = list(seen.values())

print(f"Loci with mutation rates : {len(loci)}")


# ── 2. Ancestral haplotype ────────────────────────────────────────────────────
ancestor_hap = {loc["name"]: ANCESTOR_REPEAT for loc in loci}


# ── 3. Mutation model ─────────────────────────────────────────────────────────
def mutate(parent_hap: dict) -> dict:
    """Return a child haplotype derived from parent_hap.

    Each locus mutates independently with probability = its rate.
    Mutations are single-step ±1 (equal probability each direction).
    Repeat count is floored at MIN_REPEATS.
    """
    child = parent_hap.copy()
    for loc in loci:
        if np.random.random() < loc["rate"]:
            direction = np.random.choice([-1, 1])
            child[loc["name"]] = max(MIN_REPEATS, child[loc["name"]] + direction)
    return child


# ── 4. Build family tree ──────────────────────────────────────────────────────
# Node keys are lineage path strings, e.g. "1_2_1" for gen-3 individual.
# The ancestor key is "ancestor".

def path_to_name(path: list) -> str:
    return "ancestor" if not path else "S_" + "_".join(str(x) for x in path)


tree = {}
tree["ancestor"] = {
    "name":       "ancestor",
    "path":       [],
    "haplotype":  ancestor_hap,
    "parent":     None,
    "children":   [],
    "generation": 0,
}

current_gen = ["ancestor"]

for gen in range(1, N_GENERATIONS + 1):
    next_gen = []
    for parent_id in current_gen:
        parent = tree[parent_id]
        for son_idx in range(1, SONS_PER_PARENT + 1):
            child_path = parent["path"] + [son_idx]
            child_id   = path_to_name(child_path)
            child_hap  = mutate(parent["haplotype"])
            tree[child_id] = {
                "name":       child_id,
                "path":       child_path,
                "haplotype":  child_hap,
                "parent":     parent_id,
                "children":   [],
                "generation": gen,
            }
            parent["children"].append(child_id)
            next_gen.append(child_id)
    current_gen = next_gen

print(f"Total individuals simulated: {len(tree)}")


# ── 5. Write genotype files ───────────────────────────────────────────────────
for node_id, node in tree.items():
    filepath = os.path.join(GENO_DIR, f"{node_id}.tsv")
    with open(filepath, "w", newline="") as fh:
        fh.write(f"# Sample     : {node_id}\n")
        fh.write(f"# Generation : {node['generation']}\n")
        fh.write(f"# Parent     : {node['parent']}\n")
        fh.write(f"# Children   : {', '.join(node['children']) if node['children'] else 'none'}\n")
        fh.write(f"# RandomSeed : {RANDOM_SEED}\n")
        fh.write("#\n")
        fh.write("Locus\tRepeats\n")
        for loc in loci:
            fh.write(f"{loc['name']}\t{node['haplotype'][loc['name']]}\n")

print(f"Genotype files written to : {GENO_DIR}/")


# ── 6. MRCA for all pairs in the last generation ──────────────────────────────
last_gen = sorted(
    [nid for nid, n in tree.items() if n["generation"] == N_GENERATIONS]
)

def find_mrca(id1: str, id2: str) -> tuple:
    """Return (mrca_node_id, generations_back) relative to last generation."""
    anc1 = {}
    cur = id1
    while cur is not None:
        anc1[cur] = tree[cur]["generation"]
        cur = tree[cur]["parent"]

    cur = id2
    while cur is not None:
        if cur in anc1:
            gens_back = N_GENERATIONS - tree[cur]["generation"]
            return cur, gens_back
        cur = tree[cur]["parent"]
    return None, None


mrca_rows = []
for a, b in combinations(last_gen, 2):
    mrca_id, gens_back = find_mrca(a, b)
    mrca_rows.append({
        "sample_1":      a,
        "sample_2":      b,
        "mrca_node":     mrca_id,
        "mrca_gen":      tree[mrca_id]["generation"],
        "generations_back": gens_back,
    })

mrca_path = os.path.join(OUTPUT_DIR, "mrca_pairs.csv")
with open(mrca_path, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=mrca_rows[0].keys())
    writer.writeheader()
    writer.writerows(mrca_rows)

print(f"MRCA pairs written to     : {mrca_path}")
print(f"  Total pairs : {len(mrca_rows)}")
gens_dist = {}
for row in mrca_rows:
    gb = row["generations_back"]
    gens_dist[gb] = gens_dist.get(gb, 0) + 1
for gb in sorted(gens_dist):
    print(f"  {gb} generation(s) back : {gens_dist[gb]} pair(s)")


# ── 7. Genealogical tree visualization ───────────────────────────────────────
def leaf_count(node_id: str) -> int:
    node = tree[node_id]
    if not node["children"]:
        return 1
    return sum(leaf_count(c) for c in node["children"])


def assign_x(node_id: str, x0: float, x1: float, depth: int, pos: dict):
    xmid = (x0 + x1) / 2.0
    pos[node_id] = (xmid, -depth)
    if not tree[node_id]["children"]:
        return
    total = leaf_count(node_id)
    cur = x0
    for child_id in tree[node_id]["children"]:
        frac = leaf_count(child_id) / total
        assign_x(child_id, cur, cur + (x1 - x0) * frac, depth + 1, pos)
        cur += (x1 - x0) * frac


positions = {}
assign_x("ancestor", 0.0, 1.0, 0, positions)

# Scale figure width with the number of leaves (256 for 8 gens)
n_leaves  = 2 ** N_GENERATIONS
fig_w     = max(22, n_leaves * 0.18)
fig_h     = max(10, N_GENERATIONS * 1.4)
NODE_R    = min(0.006, 0.6 / n_leaves)   # shrink nodes for large trees

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
ax.set_aspect("auto")
ax.axis("off")

palette = plt.cm.tab10(np.linspace(0, 0.7, N_GENERATIONS + 1))

# Draw edges first
for node_id, node in tree.items():
    if node["parent"] is not None:
        px, py = positions[node["parent"]]
        cx, cy = positions[node_id]
        ax.plot([px, cx], [py, cy], color="#aaaaaa", linewidth=0.5, zorder=1)

# Draw nodes and labels
for node_id, node in tree.items():
    x, y = positions[node_id]
    gen  = node["generation"]
    circ = plt.Circle((x, y), NODE_R, color=palette[gen], zorder=3,
                       linewidth=0.3, edgecolor="white")
    ax.add_patch(circ)

    label = node_id.replace("S_", "").replace("_", ".")
    if gen == 0:
        ax.text(x, y + NODE_R * 3, "Ancestor", ha="center", va="bottom",
                fontsize=8, fontweight="bold", color="#333333")
    elif gen == N_GENERATIONS and n_leaves <= 64:
        # only draw leaf labels when the tree is small enough to read
        ax.text(x, y - NODE_R * 3, label, ha="center", va="top",
                fontsize=4, rotation=60, color="#444444")
    elif gen <= 2:
        ax.text(x + NODE_R * 2, y, label, ha="left", va="center",
                fontsize=6, color="#333333")

# Generation axis labels on the left
for gen in range(N_GENERATIONS + 1):
    sample_id = next(nid for nid, n in tree.items() if n["generation"] == gen)
    _, y = positions[sample_id]
    ax.text(-0.025, y, f"Gen {gen}", ha="right", va="center",
            fontsize=8, color="#555555")

# Legend
handles = [mpatches.Patch(color=palette[g], label=f"Generation {g}")
           for g in range(N_GENERATIONS + 1)]
ax.legend(handles=handles, loc="lower right", fontsize=7,
          framealpha=0.85, ncol=3)

ax.set_title(
    f"Y-STR Family Simulation  |  seed={RANDOM_SEED}  |  "
    f"{len(loci)} loci  |  {N_GENERATIONS} generations  |  "
    f"{len(last_gen)} individuals in last gen",
    fontsize=11, pad=12
)
ax.set_xlim(-0.07, 1.05)
ax.set_ylim(-N_GENERATIONS - 0.5, 0.3)

plt.tight_layout()
tree_path = os.path.join(OUTPUT_DIR, "family_tree.png")
plt.savefig(tree_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"Family tree image saved   : {tree_path}")


# ── 8. Summary ────────────────────────────────────────────────────────────────
print("\n── Simulation summary ──────────────────────────────────────────────")
print(f"  Random seed        : {RANDOM_SEED}")
print(f"  Loci simulated     : {len(loci)}")
print(f"  Generations        : {N_GENERATIONS}")
print(f"  Sons per parent    : {SONS_PER_PARENT}")
print(f"  Total individuals  : {len(tree)}")
print(f"  Last gen size      : {len(last_gen)}")
print(f"  MRCA pairs written : {len(mrca_rows)}")
print("────────────────────────────────────────────────────────────────────")
