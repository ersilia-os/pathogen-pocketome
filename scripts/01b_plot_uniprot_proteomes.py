import argparse
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd
import stylia
from matplotlib.lines import Line2D

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

pathogens_file = os.path.join(root, "..", "src", "pathogens.csv")
proteome_dir = os.path.join(root, "..", "output", "01_uniprot_proteome")
go_terms_path = os.path.join(root, "..", "data", "raw", "go", "go_terms.csv")
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


ANNOTATION_SCORES = [1, 2, 3, 4, 5]

# GO aspect + number of terms used for the "top functions" panel (E). Fixed to the
# aspect that literally means "function" (as opposed to biological_process or
# cellular_component) — see docs/2026-09-24_uniprot-proteome-pipeline.md.
TOP_GO_ASPECT = "molecular_function"
TOP_N_GO_TERMS = 5


def compute_stats(df):
    reviewed = df["reviewed"] == "reviewed"
    return {
        "n_proteins": len(df),
        "prop_reviewed": reviewed.mean(),
        "prop_unreviewed": 1 - reviewed.mean(),
        "n_pdb": df["xref_pdb"].notna().sum(),
        "annotation_score_props": [
            (df["annotation_score"] == score).mean() for score in ANNOTATION_SCORES
        ],
    }


def determine_top_go_ids(df, go_terms, aspect=TOP_GO_ASPECT, top_n=TOP_N_GO_TERMS):
    counter = Counter()
    for value in df["go_ids"].dropna():
        counter.update(x.strip() for x in value.split(";") if x.strip())
    candidates = [
        (go_id, count) for go_id, count in counter.items()
        if go_id in go_terms.index and go_terms.loc[go_id, "aspect"] == aspect
    ]
    candidates.sort(key=lambda item: item[1], reverse=True)
    return [go_id for go_id, _ in candidates[:top_n]]


def compute_go_term_props(df, go_ids):
    filled = df["go_ids"].fillna("")
    return [filled.str.contains(go_id, regex=False).mean() for go_id in go_ids]


def load_summaries(pathogens, demo=False):
    go_terms = pd.read_csv(go_terms_path).set_index("go_id")

    labels = []
    n_proteins = []
    prop_reviewed = []
    prop_unreviewed = []
    n_pdb = []
    annotation_score_props = []
    go_term_props = []

    reference_stats = None
    reference_df = None
    top_go_ids = None
    top_go_names = None
    for _, row in pathogens.iterrows():
        code = row["pathogen_code"]
        proteins_csv = os.path.join(proteome_dir, code, "{}_proteins.csv".format(code))
        if os.path.exists(proteins_csv):
            df = pd.read_csv(proteins_csv)
            stats = compute_stats(df)
            if reference_stats is None:
                reference_stats = stats
                reference_df = df
                top_go_ids = determine_top_go_ids(df, go_terms)
                top_go_names = [go_terms.loc[go_id, "name"] for go_id in top_go_ids]
        elif demo and reference_stats is not None:
            print("{}: no real download yet, reusing reference numbers for the demo figure".format(code))
            df = reference_df
            stats = reference_stats
        else:
            print("Skipping {}: {} not found".format(code, proteins_csv))
            continue

        labels.append(row["organism_name"])
        n_proteins.append(stats["n_proteins"])
        prop_reviewed.append(stats["prop_reviewed"])
        prop_unreviewed.append(stats["prop_unreviewed"])
        n_pdb.append(stats["n_pdb"])
        annotation_score_props.append(stats["annotation_score_props"])
        go_term_props.append(compute_go_term_props(df, top_go_ids))

    return {
        "labels": labels,
        "n_proteins": n_proteins,
        "prop_reviewed": prop_reviewed,
        "prop_unreviewed": prop_unreviewed,
        "n_pdb": n_pdb,
        "annotation_score_props": annotation_score_props,
        "top_go_names": top_go_names,
        "go_term_props": go_term_props,
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


def plot_pdb_count(ax, summary, show_labels):
    nc = stylia.NamedColors()
    n = len(summary["labels"])
    positions = np.arange(n)
    ax.hlines(positions, 0, summary["n_pdb"], color=nc.turquoise)
    ax.plot(summary["n_pdb"], positions, "o", color=nc.turquoise)
    style_yaxis(ax, n, summary["labels"], show_labels)
    stylia.label(ax, xlabel="Number of proteins with PDB structure", ylabel="", title="Structural coverage", abc="C")


def plot_annotation_score_proportion(ax, summary, show_labels):
    n = len(summary["labels"])
    positions = np.arange(n)
    cm = stylia.FadingColormap("cobalt")
    cm.fit(ANNOTATION_SCORES)
    colors = cm.transform(ANNOTATION_SCORES)
    left = np.zeros(n)
    for i, (score, color) in enumerate(zip(ANNOTATION_SCORES, colors)):
        props = np.array([p[i] for p in summary["annotation_score_props"]])
        ax.barh(positions, props, left=left, color=color, label=str(score))
        left += props
    style_yaxis(ax, n, summary["labels"], show_labels)
    ax.set_xlim(0, 1)
    ax.legend(title="Score")
    stylia.label(ax, xlabel="Proportion of proteins", ylabel="", title="Annotation quality", abc="D")


def plot_top_go_dotplot(ax, summary, show_labels):
    nc = stylia.NamedColors()
    n = len(summary["labels"])
    positions = np.arange(n)
    names = summary["top_go_names"]
    x_positions = np.arange(len(names))

    xs, ys, sizes = [], [], []
    for i, props in enumerate(summary["go_term_props"]):
        for j, prop in enumerate(props):
            xs.append(x_positions[j])
            ys.append(positions[i])
            sizes.append(prop)

    size_scale = 3000
    legend_props = [0.01, 0.10]
    ax.scatter(xs, ys, s=np.array(sizes) * size_scale, color=nc.tangerine, alpha=0.7, edgecolor=nc.tangerine)
    legend_handles = [
        Line2D(
            [0], [0], marker="o", linestyle="", color=nc.tangerine, markeredgecolor=nc.tangerine,
            alpha=0.7, markersize=np.sqrt(p * size_scale), label="{:.0%}".format(p),
        )
        for p in legend_props
    ]
    ax.legend(handles=legend_handles, title="Proportion", loc="upper right")

    ax.set_xticks(x_positions)
    ax.set_xticklabels([str(i + 1) for i in range(len(names))])
    style_yaxis(ax, n, summary["labels"], show_labels)
    stylia.label(ax, xlabel="", ylabel="", title="Top functions", abc="E")


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
    plot_pdb_count(axs.next(), summary, show_labels=False)
    plot_annotation_score_proportion(axs.next(), summary, show_labels=False)
    plot_top_go_dotplot(axs.next(), summary, show_labels=False)

    filename = "01b_uniprot_overview_DEMO.png" if args.demo else "01b_uniprot_overview.png"
    output_path = os.path.join(output_dir, filename)
    stylia.save_figure(output_path)
    print("Saved figure to {}".format(output_path))


if __name__ == "__main__":
    main()
