import argparse
import os
import sys

import numpy as np
import pandas as pd
import stylia

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

pathogens_file = os.path.join(root, "..", "src", "pathogens.csv")
proteome_dir = os.path.join(root, "..", "output", "01_uniprot_proteome")
output_dir = os.path.join(root, "..", "output", "01_uniprot_proteome")
os.makedirs(output_dir, exist_ok=True)

# Format: print | Style: article — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("article")

# Pathogens not yet resolved in src/pathogens.csv (no verified taxonomy ID / reference
# proteome accession, no real download). Used only with --demo, to preview the full
# 15-row layout by reusing A. baumannii's real numbers as a stand-in. Never written to
# src/pathogens.csv or output/ — these names exist only in memory for this run.
DEMO_ONLY_PATHOGENS = [
    ("calbicans", "Candida albicans"),
    ("campylobacter", "Campylobacter spp."),
    ("ecoli", "Escherichia coli"),
    ("efaecium", "Enterococcus faecium"),
    ("enterobacter", "Enterobacter spp."),
    ("hpylori", "Helicobacter pylori"),
    ("kpneumoniae", "Klebsiella pneumoniae"),
    ("mtuberculosis", "Mycobacterium tuberculosis"),
    ("ngonorrhoeae", "Neisseria gonorrhoeae"),
    ("paeruginosa", "Pseudomonas aeruginosa"),
    ("pfalciparum", "Plasmodium falciparum"),
    ("saureus", "Staphylococcus aureus"),
    ("smansoni", "Schistosoma mansoni"),
    ("spneumoniae", "Streptococcus pneumoniae"),
]


def compute_stats(df):
    reviewed = df["reviewed"] == "reviewed"
    return {
        "n_proteins": len(df),
        "prop_reviewed": reviewed.mean(),
        "prop_unreviewed": 1 - reviewed.mean(),
        "prop_pdb": df["xref_pdb"].notna().mean(),
        "prop_pocket": (df["active_site"].notna() | df["binding_site"].notna()).mean(),
        "annotation_scores": df["annotation_score"].dropna().values,
    }


def load_summaries(pathogens, demo=False):
    labels = []
    n_proteins = []
    prop_reviewed = []
    prop_unreviewed = []
    prop_pdb = []
    prop_pocket = []
    annotation_scores = []

    reference_stats = None
    for _, row in pathogens.iterrows():
        code = row["pathogen_code"]
        proteins_csv = os.path.join(proteome_dir, code, "{}_proteins.csv".format(code))
        if os.path.exists(proteins_csv):
            stats = compute_stats(pd.read_csv(proteins_csv))
            reference_stats = reference_stats or stats
        elif demo and reference_stats is not None:
            print("{}: no real download yet, reusing reference numbers for the demo figure".format(code))
            stats = reference_stats
        else:
            print("Skipping {}: {} not found".format(code, proteins_csv))
            continue

        labels.append(row["organism_name"])
        n_proteins.append(stats["n_proteins"])
        prop_reviewed.append(stats["prop_reviewed"])
        prop_unreviewed.append(stats["prop_unreviewed"])
        prop_pdb.append(stats["prop_pdb"])
        prop_pocket.append(stats["prop_pocket"])
        annotation_scores.append(stats["annotation_scores"])

    return {
        "labels": labels,
        "n_proteins": n_proteins,
        "prop_reviewed": prop_reviewed,
        "prop_unreviewed": prop_unreviewed,
        "prop_pdb": prop_pdb,
        "prop_pocket": prop_pocket,
        "annotation_scores": annotation_scores,
    }


def style_yaxis(ax, n, labels, show_labels):
    positions = np.arange(n)
    ax.set_yticks(positions)
    ax.set_yticklabels(labels if show_labels else [""] * n)
    ax.invert_yaxis()


def plot_protein_counts(ax, summary, show_labels):
    nc = stylia.NamedColors()
    n = len(summary["labels"])
    positions = np.arange(n)
    ax.hlines(positions, 0, summary["n_proteins"], color=nc.cobalt)
    ax.plot(summary["n_proteins"], positions, "o", color=nc.cobalt)
    style_yaxis(ax, n, summary["labels"], show_labels)
    stylia.label(ax, xlabel="Number of proteins", ylabel="", title="Protein sequences", abc="A")


def plot_reviewed_stack(ax, summary, show_labels):
    nc = stylia.NamedColors()
    n = len(summary["labels"])
    positions = np.arange(n)
    ax.barh(positions, summary["prop_reviewed"], color=nc.crimson, label="Reviewed")
    ax.barh(
        positions,
        summary["prop_unreviewed"],
        left=summary["prop_reviewed"],
        color=nc.silver,
        label="Unreviewed",
    )
    style_yaxis(ax, n, summary["labels"], show_labels)
    ax.set_xlim(0, 1)
    ax.legend()
    stylia.label(ax, xlabel="Proportion of proteins", ylabel="", title="Reviewed status", abc="B")


def plot_pdb_proportion(ax, summary, show_labels):
    nc = stylia.NamedColors()
    n = len(summary["labels"])
    positions = np.arange(n)
    ax.barh(positions, summary["prop_pdb"], color=nc.turquoise)
    style_yaxis(ax, n, summary["labels"], show_labels)
    ax.set_xlim(0, 1)
    stylia.label(ax, xlabel="Proportion with PDB structure", ylabel="", title="Structural coverage", abc="C")


def plot_annotation_score(ax, summary, show_labels):
    nc = stylia.NamedColors()
    n = len(summary["labels"])
    positions = np.arange(n)
    ax.boxplot(
        summary["annotation_scores"],
        positions=positions,
        vert=False,
        widths=0.6,
        patch_artist=True,
        boxprops=dict(facecolor=nc.get("cobalt", lighten=0.5), color=nc.cobalt),
        medianprops=dict(color=nc.cobalt),
        whiskerprops=dict(color=nc.cobalt),
        capprops=dict(color=nc.cobalt),
        flierprops=dict(markeredgecolor=nc.cobalt),
    )
    style_yaxis(ax, n, summary["labels"], show_labels)
    stylia.label(ax, xlabel="Annotation score", ylabel="", title="Annotation quality", abc="D")


def plot_pocket_proportion(ax, summary, show_labels):
    nc = stylia.NamedColors()
    n = len(summary["labels"])
    positions = np.arange(n)
    ax.barh(positions, summary["prop_pocket"], color=nc.tangerine)
    style_yaxis(ax, n, summary["labels"], show_labels)
    ax.set_xlim(0, 1)
    stylia.label(ax, xlabel="Proportion with pocket residues", ylabel="", title="Pocket annotation", abc="E")


def build_pathogens(demo):
    pathogens = pd.read_csv(pathogens_file)
    if demo:
        demo_rows = pd.DataFrame(DEMO_ONLY_PATHOGENS, columns=["pathogen_code", "organism_name"])
        pathogens = pd.concat([pathogens, demo_rows], ignore_index=True)
    return pathogens


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--demo",
        action="store_true",
        help=(
            "Preview the full multi-pathogen layout by filling in pathogens with no real "
            "download yet using A. baumannii's real numbers as a stand-in. In-memory only "
            "— never written to src/pathogens.csv or output/."
        ),
    )
    args = parser.parse_args()

    pathogens = build_pathogens(args.demo)
    summary = load_summaries(pathogens, demo=args.demo)

    fig, axs = stylia.create_figure(1, 5)
    plot_protein_counts(axs.next(), summary, show_labels=True)
    plot_reviewed_stack(axs.next(), summary, show_labels=False)
    plot_pdb_proportion(axs.next(), summary, show_labels=False)
    plot_annotation_score(axs.next(), summary, show_labels=False)
    plot_pocket_proportion(axs.next(), summary, show_labels=False)

    filename = "01b_uniprot_overview_DEMO.png" if args.demo else "01b_uniprot_overview.png"
    output_path = os.path.join(output_dir, filename)
    stylia.save_figure(output_path)
    print("Saved figure to {}".format(output_path))


if __name__ == "__main__":
    main()
