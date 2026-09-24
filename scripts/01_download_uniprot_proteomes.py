import os
import sys
import zipfile

import pandas as pd
import requests

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import UNIPROT_REST_BASE, UNIPROT_PROTEIN_FIELDS

pathogens_file = os.path.join(root, "..", "src", "pathogens.csv")
raw_dir = os.path.join(root, "..", "data", "raw", "uniprot")
output_dir = os.path.join(root, "..", "output", "01_uniprot_proteome")
tmp_dir = os.path.join(root, "..", "tmp", "01_uniprot_proteome")
os.makedirs(raw_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)
os.makedirs(tmp_dir, exist_ok=True)

COLUMN_RENAME = {
    "Entry": "uniprot_ac",
    "Entry Name": "entry_name",
    "Protein names": "protein_name",
    "Gene Names": "gene_names",
    "Organism (ID)": "organism_id",
    "Length": "length",
    "Mass": "mass",
    "Reviewed": "reviewed",
    "Protein existence": "protein_existence",
    "Annotation": "annotation_score",
    "PDB": "xref_pdb",
    "AlphaFoldDB": "xref_alphafolddb",
    "InterPro": "xref_interpro",
    "Gene Ontology IDs": "go_ids",
    "EC number": "ec_number",
    "Active site": "active_site",
    "Binding site": "binding_site",
    "Subcellular location [CC]": "subcellular_location",
    "Transmembrane": "transmembrane",
}


def download_proteome_tsv(proteome_id, dest_path):
    # /uniprotkb/search paginates to 25 results by default; /stream returns the full set.
    url = "{}/uniprotkb/stream".format(UNIPROT_REST_BASE)
    params = {
        "query": "proteome:{}".format(proteome_id),
        "format": "tsv",
        "fields": ",".join(UNIPROT_PROTEIN_FIELDS),
    }
    response = requests.get(url, params=params, timeout=120)
    response.raise_for_status()
    with open(dest_path, "w") as f:
        f.write(response.text)


def download_proteome_fasta(proteome_id, dest_path):
    url = "{}/uniprotkb/stream".format(UNIPROT_REST_BASE)
    params = {
        "query": "proteome:{}".format(proteome_id),
        "format": "fasta",
    }
    response = requests.get(url, params=params, timeout=120)
    response.raise_for_status()
    with open(dest_path, "w") as f:
        f.write(response.text)


def split_fasta(fasta_path, dest_dir):
    with open(fasta_path) as f:
        content = f.read()

    records = ["\n>" + block if i > 0 else block for i, block in enumerate(content.split("\n>"))]
    ac_paths = []
    for record in records:
        record = record.strip()
        if not record:
            continue
        header = record.splitlines()[0]
        accession = header.split("|")[1]
        record_path = os.path.join(dest_dir, "{}.fasta".format(accession))
        with open(record_path, "w") as f:
            f.write(record + "\n")
        ac_paths.append(record_path)
    return ac_paths


def zip_fasta_files(fasta_paths, dest_zip):
    with zipfile.ZipFile(dest_zip, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in fasta_paths:
            zf.write(path, arcname=os.path.basename(path))


def main():
    pathogens = pd.read_csv(pathogens_file)

    for _, row in pathogens.iterrows():
        pathogen_code = row["pathogen_code"]
        proteome_id = row["reference_proteome_id"]

        pathogen_raw_dir = os.path.join(raw_dir, pathogen_code)
        pathogen_output_dir = os.path.join(output_dir, pathogen_code)
        pathogen_tmp_dir = os.path.join(tmp_dir, pathogen_code)
        os.makedirs(pathogen_raw_dir, exist_ok=True)
        os.makedirs(pathogen_output_dir, exist_ok=True)
        os.makedirs(pathogen_tmp_dir, exist_ok=True)

        raw_tsv_path = os.path.join(pathogen_raw_dir, "{}_{}.tsv".format(pathogen_code, proteome_id))
        raw_fasta_path = os.path.join(pathogen_raw_dir, "{}_{}.fasta".format(pathogen_code, proteome_id))
        download_proteome_tsv(proteome_id, raw_tsv_path)
        download_proteome_fasta(proteome_id, raw_fasta_path)

        df = pd.read_csv(raw_tsv_path, sep="\t")
        df = df.rename(columns=COLUMN_RENAME)
        df.insert(0, "pathogen_code", pathogen_code)
        proteins_csv_path = os.path.join(pathogen_output_dir, "{}_proteins.csv".format(pathogen_code))
        df.to_csv(proteins_csv_path, index=False)

        fasta_paths = split_fasta(raw_fasta_path, pathogen_tmp_dir)
        sequences_zip_path = os.path.join(pathogen_output_dir, "{}_sequences.zip".format(pathogen_code))
        zip_fasta_files(fasta_paths, sequences_zip_path)
        for path in fasta_paths:
            os.remove(path)
        os.rmdir(pathogen_tmp_dir)

        print(
            "{}: {} proteins written ({} sequences zipped)".format(
                pathogen_code, len(df), len(fasta_paths)
            )
        )


if __name__ == "__main__":
    main()
