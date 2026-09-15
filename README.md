# Pathogen Pocketome

[![status](https://img.shields.io/badge/status-pending-red)](https://github.com/)

The pathogen pocketome.

## Background

Starting from canonical protein sequences for 15 pathogens, the idea is to build up a structural annotation (PDB structures and predicted models from AF2, AF3, and homology modelling) and a functional annotation (Interpro, UniProt), then move into binding site identification, drawing on ligand evidence, AlphaFill, etc., to arrive at the pathogenic pocketome itself. From there, the plan would involve clustering billion-sized libraries down to something more tractable, maybe 10–50k compounds, and docking them against all the identified pockets, which will need to be GPU-accelerated. On top of that, surrogate modelling is envisioned (potentially with LQ, though faster inference may be needed), coupled with generative models in an iterative training-learning loop. The project is meant to be complementary to Ersilia's existing ChEMBL models, which are phenotypic, and also to target-based ChEMBL models, which are still to be defined.

This is still very preliminary.

## Project tracking

* No meetings logged yet.

## Progress

Not started yet.

## Repository structure

See [CLAUDE.md](CLAUDE.md) for the repository structure and conventions.

## About the Ersilia Open Source Initiative

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling sustainable research in the Global South. Ersilia's main asset is the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia), an open-source repository of AI/ML models for antimicrobial drug discovery.

![Ersilia Logo](assets/Ersilia_Brand.png)
