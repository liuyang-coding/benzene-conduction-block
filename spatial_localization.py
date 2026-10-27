# Reference-heart spatial localization
# Author: Yang Liu | Date: 2026-10-27 | Python: 3.13.11
# Inputs: {section_id: AnnData}, candidates DataFrame with unique symbol.
# Each atlas: author log1p(CP10k) X, uns['log1p'], gene-symbol var_names,
# obs donor/region/annotation_final/total_counts, obsm['X_spatial'] coordinates.
# Dependencies: numpy, pandas, scipy, statsmodels, matplotlib; anndata for inputs.

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import ttest_1samp
from statsmodels.stats.multitest import multipletests


# Section and donor summaries ----
def spatial_summaries(sections, candidates, *, author_log1p_cp10k, chunk_size=5000):
    """Preserve author expression scale; average sections within donor, then donors."""
    if author_log1p_cp10k is not True:
        raise ValueError("Author log1p(CP10k) values are required and must not be renormalized")
    records, all_spots = [], []
    for section_id, atlas in sections.items():
        if "log1p" not in atlas.uns:
            raise ValueError("Author log1p marker is required")
        if atlas.obs[["donor", "region", "annotation_final"]].isna().any().any():
            raise ValueError(f"{section_id}: donor, region and structure cannot be missing")
        coordinates = np.asarray(atlas.obsm["X_spatial"])
        if coordinates.shape != (atlas.n_obs, 2):
            raise ValueError(f"{section_id}: expected n_spots by 2 coordinates")
        mapping = spatial_gene_mapping(atlas, candidates)
        positions = mapping.position.to_numpy()
        found = positions >= 0
        values = np.full((atlas.n_obs, len(mapping)), np.nan)
        for start in range(0, atlas.n_obs, chunk_size):
            end = min(start + chunk_size, atlas.n_obs)
            matrix = atlas.X[start:end, :]
            raw = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
            if not np.isfinite(raw).all() or np.any(raw < 0):
                raise ValueError(f"{section_id}: expected finite nonnegative log1p(CP10k)")
            selected = matrix[:, positions[found]]
            values[start:end, found] = selected.toarray() if sparse.issparse(selected) else np.asarray(selected)
        obs = atlas.obs[["donor", "region", "annotation_final", "total_counts"]].copy()
        obs = obs.rename(columns={"annotation_final": "structure"})
        for column in ["donor", "region", "structure"]:
            obs[column] = obs[column].astype(str)
            if obs[column].str.strip().eq("").any():
                raise ValueError(f"{section_id}: empty {column}")
        if obs.donor.nunique() != 1 or obs.region.nunique() != 1:
            raise ValueError(f"{section_id}: each section must have exactly one donor and region")
        obs.loc[(obs.region == "AVN") & (obs.structure == "node"), "structure"] = "AV_node"
        obs.loc[(obs.region == "SAN") & (obs.structure == "node"), "structure"] = "SA_node"
        obs.loc[obs.structure.str.lower().isin(["purkinje cell", "purkinje"]), "structure"] = "Purkinje"
        counts = pd.to_numeric(obs.total_counts, errors="coerce").to_numpy(dtype=float)
        valid = np.isfinite(coordinates).all(axis=1) & np.isfinite(counts) & (counts > 0)
        obs["section"], obs["spot"] = section_id, obs.index.astype(str)
        obs["x"], obs["y"] = coordinates[:, 0], coordinates[:, 1]
        obs, values = obs.loc[valid].reset_index(drop=True), values[valid]
        if obs.empty:
            raise ValueError(f"{section_id}: no eligible spots")
        symbols = mapping.symbol.tolist()
        all_spots.append(pd.concat([obs, pd.DataFrame(values, columns=symbols)], axis=1))
        for structure, rows in list(obs.groupby("structure", observed=True).groups.items()) + [("All", obs.index)]:
            array = values[np.asarray(rows)]
            for index, symbol in enumerate(symbols):
                finite = np.isfinite(array[:, index])
                records.append({"section": section_id, "donor": obs.donor.iloc[0],
                    "region": obs.region.iloc[0], "structure": structure, "symbol": symbol,
                    "spots": len(rows),
                    "mean_log1p_cp10k": array[finite, index].mean() if finite.any() else np.nan,
                    "percent_expressing": (array[finite, index] > 0).mean() * 100 if finite.any() else np.nan})
    by_section = pd.DataFrame(records)
    by_donor = by_section.groupby(["donor", "region", "structure", "symbol"], sort=False).agg(
        mean_log1p_cp10k=("mean_log1p_cp10k", "mean"), percent_expressing=("percent_expressing", "mean"),
        sections=("section", "nunique"), spots=("spots", "sum")).reset_index()
    summary = by_donor.groupby(["region", "structure", "symbol"], sort=False).agg(
        mean_log1p_cp10k=("mean_log1p_cp10k", "mean"), percent_expressing=("percent_expressing", "mean"),
        donors=("donor", "nunique")).reset_index()
    return {"by_section": by_section, "by_donor": by_donor, "summary": summary,
            "spots": pd.concat(all_spots, ignore_index=True)}
# result = spatial_summaries(sections, candidates, author_log1p_cp10k=True)


# Conduction structures versus paired myocardium ----
def spatial_contrasts(result):
    """>=10 spots/side in one section; >=3 donors and nonzero variance for t tests."""
    sections = result["by_section"]
    genes = sections.symbol.drop_duplicates().tolist()
    contrasts = {"AV_node": "myocardium_atrial", "AV_bundle": "myocardium_ventricular",
                 "Purkinje": "myocardium_ventricular"}
    paired_rows = []
    for target, reference in contrasts.items():
        left = sections.loc[(sections.structure == target) & (sections.spots >= 10)]
        right = sections.loc[(sections.structure == reference) & (sections.spots >= 10)]
        pairs = left.merge(right, on=["section", "donor", "region", "symbol"],
                           suffixes=("_target", "_reference"), validate="one_to_one")
        for row in pairs.itertuples():
            paired_rows.append({"contrast": target, "section": row.section, "donor": row.donor,
                "symbol": row.symbol, "difference": row.mean_log1p_cp10k_target - row.mean_log1p_cp10k_reference})
    paired = pd.DataFrame(paired_rows, columns=["contrast", "section", "donor", "symbol", "difference"])
    donor = paired.groupby(["contrast", "donor", "symbol"]).difference.mean().reset_index()
    rows = []
    for contrast in contrasts:
        for symbol in genes:
            values = donor.loc[(donor.contrast == contrast) & (donor.symbol == symbol), "difference"].dropna().to_numpy(dtype=float)
            pvalue = float(ttest_1samp(values, 0).pvalue) if len(values) >= 3 and np.std(values) > 0 else np.nan
            rows.append({"contrast": contrast, "symbol": symbol, "donors": len(values),
                "mean_difference": values.mean() if len(values) else np.nan, "p_value": pvalue,
                "status": "tested" if np.isfinite(pvalue) else "insufficient donors or zero variance"})
    tests = pd.DataFrame(rows)
    tests["BH_FDR"] = np.nan
    finite = np.isfinite(tests.p_value)
    if finite.any():
        # One family across all finite gene-by-three-structure tests.
        tests.loc[finite, "BH_FDR"] = multipletests(tests.loc[finite, "p_value"], method="fdr_bh")[1]
    return {"paired_sections": paired, "paired_donors": donor, "tests": tests}
# comparisons = spatial_contrasts(result)


# Single-section expression plot ----
def spatial_expression_plot(result, section_id, symbol, *, image=None, image_scale=1.0, outline_structure=None):
    """Return a Figure; optional author histology/scale and display-only pooled P99 cap."""
    import matplotlib.pyplot as plt

    symbol = str(symbol).strip().upper()
    spots = result["spots"]
    selected = spots.loc[spots.section == section_id]
    if selected.empty or not np.isfinite(image_scale) or image_scale <= 0:
        raise ValueError("Unknown section or invalid image_scale")
    pooled = spots[symbol].to_numpy(dtype=float)
    pooled = pooled[np.isfinite(pooled)]
    vmax = max(float(np.quantile(pooled, 0.99)), 1e-6) if len(pooled) else 1.0
    fig, ax = plt.subplots(figsize=(6, 6))
    if image is not None:
        ax.imshow(image)
    detected = np.isfinite(selected[symbol]) & (selected[symbol] > 0)
    artist = ax.scatter(selected.loc[detected, "x"] * image_scale,
                        selected.loc[detected, "y"] * image_scale,
                        c=selected.loc[detected, symbol], s=7, cmap="magma", vmin=0, vmax=vmax, alpha=0.8)
    if not np.isfinite(selected[symbol]).any():
        ax.text(0.5, 0.5, "Gene unavailable", transform=ax.transAxes, ha="center")
    if outline_structure is not None:
        target = selected.structure == outline_structure
        ax.scatter(selected.loc[target, "x"] * image_scale, selected.loc[target, "y"] * image_scale,
                   s=12, facecolors="none", edgecolors="#00BFC4", linewidths=0.5)
    if image is None:
        ax.invert_yaxis()
    else:
        ax.set_xlim(0, image.shape[1])
        ax.set_ylim(image.shape[0], 0)
    ax.set_aspect("equal")
    ax.set_title(f"{symbol}: {section_id}")
    ax.axis("off")
    fig.colorbar(artist, ax=ax, label="Author log1p(CP10k); display capped at pooled P99", shrink=0.7)
    fig.tight_layout()
    return fig
# fig = spatial_expression_plot(result, "section_1", "NPPB", outline_structure="AV_node")
# fig.savefig("spatial.tiff", dpi=300, pil_kwargs={"compression": "tiff_lzw"})

# Gene mapping ----
def spatial_gene_mapping(atlas, candidates):
    """Preserve unmatched assays as NA; composites are never split into genes."""
    if candidates.empty or candidates.symbol.isna().any():
        raise ValueError("candidates requires a nonempty, nonmissing symbol column")
    symbols = candidates.symbol.astype(str).str.strip().str.upper()
    if symbols.eq("").any() or symbols.duplicated().any():
        raise ValueError("Candidate symbols must be unique and nonempty")
    names = pd.Index(atlas.var_names.astype(str))
    aliases = {"HLA_E": "HLA-E", "NTPROBNP": "NPPB", "CCN1": "CYR61", "SEPTIN8": "SEPT8"}
    records = []
    for symbol in symbols:
        position, matched = -1, None
        for candidate in (symbol, aliases.get(symbol)):
            if candidate is None:
                continue
            hits = np.flatnonzero(names.str.upper() == candidate.upper())
            if len(hits) > 1:
                raise ValueError(f"Ambiguous atlas symbol: {symbol}")
            if len(hits):
                position, matched = int(hits[0]), str(names[hits[0]])
                break
        records.append({"symbol": symbol, "atlas_symbol": matched, "position": position})
    return pd.DataFrame(records)
