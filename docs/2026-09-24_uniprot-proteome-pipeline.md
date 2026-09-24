# UniProt proteome pipeline — design decisions

Decision log for `scripts/01_download_uniprot_proteomes.py`, the first step of the pocketome pipeline (per-pathogen protein tables + sequences pulled from UniProt).

## Reference proteome only, not every strain

Most of the 15 pathogens have many sequenced strains in UniProt, each with its own proteome. We use only the designated **reference proteome** per pathogen, to get one non-redundant protein set per organism rather than strain-level duplication.

For *A. baumannii*, querying `https://rest.uniprot.org/proteomes/search?query=organism_id:470` (470 = the species-level NCBI taxonomy ID) returns every proteome registered under that species-level taxid — strain-specific proteomes (e.g. AB0057, ATCC 19606) are registered under their own strain-level taxids and don't appear. Of the proteomes under taxid 470, exactly one is flagged `"Reference proteome"`: `UP000032746` (strain AB5075-UW, 3839 proteins). So "species-level taxonomy ID + reference-proteome flag" resolves deterministically to one proteome, with no manual strain selection.

This was verified live against the UniProt REST API, not assumed — see the query above. The same two-step lookup (resolve taxonomy ID → filter proteomes by that ID for the reference flag) will need to be repeated for each of the other 14 pathogens before they're added to `src/pathogens.csv`.

## Reviewed + unreviewed entries

Swiss-Prot (reviewed) curation is sparse outside model organisms. Restricting to reviewed-only would drop the large majority of most of these pathogens' proteomes. We keep both reviewed and unreviewed (TrEMBL) entries, and record the `reviewed` flag as a column so downstream steps can filter or weight by curation confidence without us having dropped anything upfront.

## Field selection

Fields pulled per protein, grouped by purpose:

- **Core identifiers**: UniProt AC, entry name, gene name(s), protein name, organism ID, length, mass.
- **Structural evidence** (feeds the pocketome's structural-annotation stage): PDB cross-references, AlphaFoldDB cross-reference.
- **Functional annotation** (matches the README's "Interpro, UniProt" functional layer): InterPro domain cross-references, GO term IDs, EC number.
- **Binding-site relevant curation**: UniProt-curated active-site and binding-site residue positions (populated for a subset of well-studied proteins, e.g. OXA-23 beta-lactamase), subcellular location, transmembrane region annotations.
- **Quality flags**: reviewed status, protein existence level, annotation score — needed since we're keeping unreviewed entries and may want to triage confidence downstream.

## Sequence storage

FASTA sequences are downloaded once per pathogen as a single multi-sequence file (one UniProtKB API call), then split into individual per-protein FASTA files and bundled into a single zip per pathogen. This keeps the metadata CSV free of long sequence strings while still giving downstream steps (e.g. per-protein structure prediction or docking) easy access to individual sequences without re-parsing a multi-FASTA.

## Mapping file scope

`src/pathogens.csv` stores both the NCBI taxonomy ID and the pre-resolved reference-proteome accession per pathogen, rather than resolving the proteome ID at run time. This locks in a specific `UPxxxxxxx` for reproducibility; if UniProt ever re-designates a pathogen's reference proteome, the row needs a manual update (and a note of the change).
