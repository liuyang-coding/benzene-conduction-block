# Reference-heart single-cell expression
# Author: Yang Liu | Date: 2026-10-27 | Python: 3.13.11
# Input: AnnData with full-library raw-count X, gene-symbol var_names and
# obs donor/cell_type/scrublet_score; candidates DataFrame with unique symbol.
# Dependencies: numpy, pandas, scipy, matplotlib; anndata to load inputs.

import numpy as np
import pandas as pd
from scipy import sparse


# Doublet filtering and donor-balanced expression ----
def single_cell_summary(atlas, candidates, *, counts_are_full_library, chunk_size=10000):
    """Use full-library CP10k/log1p, then equal-weight donors within cell type."""
    if counts_are_full_library is not True:
        raise ValueError("Full-library raw counts are required for the CP10k denominator")
    if atlas.obs.donor.isna().any():
        raise ValueError("Donor metadata cannot be missing")
    mapping = single_cell_gene_mapping(atlas, candidates)
    positions = mapping.position.to_numpy()
    found = positions >= 0
    donors = atlas.obs.donor.astype(str).to_numpy()
    cell_types = atlas.obs.cell_type.astype(str).to_numpy()
    scores = pd.to_numeric(atlas.obs.scrublet_score, errors="coerce").to_numpy(dtype=float)
    eligible = np.isfinite(scores) & (scores < 0.25)
    eligible &= ~atlas.obs.cell_type.isna().to_numpy()
    eligible &= ~np.isin(cell_types, ["doublets", "NotAssigned", "nan", ""])
    eligible &= np.array([bool(value.strip()) and value != "nan" for value in donors])
    if not eligible.any():
        raise ValueError("No eligible single cells remain")
    groups = sorted(set(zip(donors[eligible], cell_types[eligible])))
    lookup = {group: index for index, group in enumerate(groups)}
    sums = np.zeros((len(groups), len(mapping)))
    detections = np.zeros_like(sums)
    counts = np.zeros(len(groups), dtype=int)
    for start in range(0, atlas.n_obs, chunk_size):
        end = min(start + chunk_size, atlas.n_obs)
        matrix = atlas.X[start:end, :]
        raw = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
        if not np.isfinite(raw).all() or np.any(raw < 0) or not np.allclose(raw, np.round(raw)):
            raise ValueError("X must contain raw nonnegative integer counts")
        library = np.asarray(matrix.sum(axis=1)).ravel()
        keep = eligible[start:end] & (library > 0)
        selected = matrix[:, positions[found]]
        selected = selected.toarray() if sparse.issparse(selected) else np.asarray(selected)
        normalized = np.log1p(selected * 10000.0 / np.maximum(library[:, None], 1.0))
        local_donors, local_types = donors[start:end], cell_types[start:end]
        for group in set(zip(local_donors[keep], local_types[keep])):
            rows = keep & (local_donors == group[0]) & (local_types == group[1])
            index = lookup[group]
            sums[index, found] += normalized[rows].sum(axis=0)
            detections[index, found] += (selected[rows] > 0).sum(axis=0)
            counts[index] += int(rows.sum())
    if not counts.sum():
        raise ValueError("No eligible cells have positive library size")
    records = []
    for index, (donor, cell_type) in enumerate(groups):
        if not counts[index]:
            continue
        for gene, symbol in enumerate(mapping.symbol):
            records.append({"donor": donor, "cell_type": cell_type, "symbol": symbol,
                "cells": counts[index],
                "mean_log1p_cp10k": sums[index, gene] / counts[index] if found[gene] else np.nan,
                "percent_expressing": detections[index, gene] / counts[index] * 100 if found[gene] else np.nan})
    by_donor = pd.DataFrame(records)
    summary = by_donor.groupby(["symbol", "cell_type"], sort=False).agg(
        mean_log1p_cp10k=("mean_log1p_cp10k", "mean"),
        percent_expressing=("percent_expressing", "mean"), donors=("donor", "nunique")).reset_index()
    types = sorted(summary.cell_type.unique())
    expression = summary.pivot(index="symbol", columns="cell_type", values="mean_log1p_cp10k").reindex(
        index=mapping.symbol, columns=types)
    detection = summary.pivot(index="symbol", columns="cell_type", values="percent_expressing").reindex(
        index=mapping.symbol, columns=types)
    return {"expression": expression, "detection": detection, "by_donor": by_donor,
            "summary": summary}
# result = single_cell_summary(atlas, candidates, counts_are_full_library=True)


# Expression heatmap ----
def single_cell_heatmap(result, *, standardize=True):
    """Return a Figure; within-gene Z scores are display-only and NA stays gray."""
    import matplotlib.pyplot as plt

    values = result["expression"].copy()
    if standardize:
        deviation = values.std(axis=1, ddof=0).replace(0, 1)
        values = values.sub(values.mean(axis=1), axis=0).div(deviation, axis=0)
    cmap = plt.get_cmap("RdBu_r" if standardize else "viridis").copy()
    cmap.set_bad("#DDDDDD")
    fig, ax = plt.subplots(figsize=(max(6, len(values.columns) * 0.7), max(3, len(values) * 0.3)))
    image = ax.imshow(np.ma.masked_invalid(values.to_numpy(dtype=float)), aspect="auto", cmap=cmap,
                      vmin=-2 if standardize else 0, vmax=2 if standardize else None)
    ax.set_xticks(range(len(values.columns)), values.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(values.index)), values.index)
    ax.set_title("Reference heart cell expression (equal-weight donors)")
    fig.colorbar(image, ax=ax, label="Within-gene Z score" if standardize else "Mean log1p(CP10k)")
    fig.tight_layout()
    return fig
# fig = single_cell_heatmap(result)
# fig.savefig("expression.tiff", dpi=300, pil_kwargs={"compression": "tiff_lzw"})

# Gene mapping ----
def single_cell_gene_mapping(atlas, candidates):
    """Map exact symbols/known single-gene aliases; never split composite assays."""
    if candidates.empty or candidates.symbol.isna().any():
        raise ValueError("Provide nonmissing candidate symbols")
    symbols = candidates.symbol.astype(str).str.strip().str.upper()
    if symbols.eq("").any() or symbols.duplicated().any():
        raise ValueError("Candidate symbols must be unique and nonempty")
    names = pd.Index(atlas.var_names.astype(str))
    upper = names.str.upper()
    aliases = {"HLA_E": "HLA-E", "NTPROBNP": "NPPB", "CCN1": "CYR61", "SEPTIN8": "SEPT8"}
    records = []
    for symbol in symbols:
        position, matched = -1, None
        for candidate in (symbol, aliases.get(symbol)):
            if candidate is None:
                continue
            hits = np.flatnonzero(upper == candidate.upper())
            if len(hits) > 1:
                raise ValueError(f"Ambiguous atlas symbol: {symbol}")
            if len(hits):
                position, matched = int(hits[0]), str(names[hits[0]])
                break
        records.append({"symbol": symbol, "atlas_symbol": matched, "position": position})
    return pd.DataFrame(records)
