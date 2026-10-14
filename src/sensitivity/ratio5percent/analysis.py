"""Michel template interface, ratio measurement and pure-pion exclusion calculations."""

from pathlib import Path

import numpy as np
import uproot


def normalize(values):
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError("Histogram must contain finite, nonnegative bin counts")
    if values.sum() <= 0:
        raise ValueError("Histogram integral must be positive")
    return values / values.sum()


def mix_templates(pion, kaon, ratios):
    """Return normalized shapes (pi + r*K)/(1+r), not fluctuated toy data."""
    pion, kaon = normalize(pion), normalize(kaon)
    ratios = np.atleast_1d(np.asarray(ratios, dtype=float))
    if pion.shape != kaon.shape:
        raise ValueError("Pion and kaon templates must have the same bin count")
    if not np.all(np.isfinite(ratios)) or np.any(ratios < 0):
        raise ValueError("K/pi ratios must be finite and nonnegative")
    return (pion[None, :] + ratios[:, None] * kaon[None, :]) / (1 + ratios[:, None])


def weighted_chi2(observed_shape, expected_shapes, events):
    """ROOT WW CHI2 for equal totals N and bin variances equal to counts.

    With both histograms normalized to N and errors sqrt(count), ROOT's
    weighted-histogram homogeneity statistic reduces to sum((n-mu)^2/(n+mu)).
    This is not a one-sample Neyman chi-square or a profile likelihood.
    """
    if not np.isfinite(events) or events <= 0:
        raise ValueError("Event count must be finite and positive")
    observed = normalize(observed_shape) * events
    expected = np.asarray([normalize(row) for row in expected_shapes]) * events
    variance = observed[None, :] + expected
    contributions = np.divide((observed[None, :] - expected) ** 2, variance,
                              out=np.zeros_like(expected), where=variance > 0)
    return contributions.sum(axis=1)


def load_templates(template_path, pion_name, kaon_name):
    """Interface: nonnegative TH1 cos(theta) templates with identical [-1,1] bins.

    Template counts are arbitrary; normalize each species independently downstream.
    Experimental event counts must never be inferred from template integrals.
    """
    with uproot.open(template_path) as source:
        pion, edges = source[pion_name].to_numpy()
        kaon, kaon_edges = source[kaon_name].to_numpy()
    if not np.array_equal(edges, kaon_edges):
        raise ValueError("Pion and kaon bin edges differ")
    if edges[0] != -1 or edges[-1] != 1 or np.any(np.diff(edges) <= 0):
        raise ValueError("Michel cos(theta) bin edges must span [-1,1] and increase")
    return normalize(pion), normalize(kaon), edges


def generate_distributions(template_path, output_path, true_ratio=0.05,
                           scan_start=0.01, scan_end=0.10, points=1000,
                           pion_name="piMichelDis", kaon_name="kaMichelDis"):
    if points < 4 or not 0 <= scan_start < scan_end or not scan_start <= true_ratio <= scan_end:
        raise ValueError("Need >=4 scan points and 0 <= start <= true ratio <= end")
    pion, kaon, edges = load_templates(template_path, pion_name, kaon_name)
    ratios = np.linspace(scan_start, scan_end, points)
    shapes = mix_templates(pion, kaon, ratios)
    observed = mix_templates(pion, kaon, [true_ratio])[0]
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(output_path) as target:
        target["simulated_data"] = (observed, edges)
        for index, shape in enumerate(shapes):
            target[f"simulated_data_{index}"] = (shape, edges)
        target.mktree("scan_config", {"physical_ratio": ratios})
    return ratios


def scan_distributions(input_path, output_path, events=10_000_000, ratio_axis="notebook"):
    if ratio_axis not in ("notebook", "physical"):
        raise ValueError("ratio_axis must be notebook or physical")
    with uproot.open(input_path) as source:
        observed, edges = source["simulated_data"].to_numpy()
        ratios = source["scan_config"]["physical_ratio"].array(library="np")
        shapes = []
        for index in range(len(ratios)):
            shape, candidate_edges = source[f"simulated_data_{index}"].to_numpy()
            if not np.array_equal(edges, candidate_edges):
                raise ValueError("Scan histogram bin edges differ")
            shapes.append(shape)
    values = weighted_chi2(observed, shapes, events)
    # Reproduce the historical (index-1)*(end-start)/points mapping explicitly.
    axis = ratios if ratio_axis == "physical" else (
        ratios[0] + (np.arange(len(ratios)) - 1) * (ratios[-1] - ratios[0]) / len(ratios)
    )
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(output_path) as target:
        target.mktree("tree", {"ratio": axis, "result": values, "physical_ratio": ratios})
        target["analysis_config"] = f"events={events}; ratio_axis={ratio_axis}; statistic=ROOT WW CHI2 (equal totals and Poisson bin errors)"
    return axis, values


DEFAULT_EXCLUSION_RATIOS = (0.01, 0.02, 0.03, 0.05, 0.07, 0.10, 0.15, 0.17, 0.20)


def required_events(pion, kaon, ratios, threshold=9.0):
    """Solve N=threshold/chi2(N=1), using WW chi-square's linear scaling.

    This reproduces the legacy threshold rule, not a calibrated coverage claim.
    """
    ratios = np.asarray(ratios, dtype=float)
    if (ratios.ndim != 1 or len(ratios) < 2 or not np.all(np.isfinite(ratios))
            or np.any(ratios <= 0) or np.any(np.diff(ratios) <= 0)):
        raise ValueError("Exclusion ratios must be finite, positive, strictly increasing")
    if not np.isfinite(threshold) or threshold <= 0:
        raise ValueError("Exclusion chi-square threshold must be positive")
    shapes = mix_templates(pion, kaon, ratios)
    coefficients = weighted_chi2(pion, shapes, 1)
    if np.any(coefficients <= 0):
        raise ValueError("Templates cannot distinguish pure pion from the mixtures")
    events = threshold / coefficients
    if not np.all(np.isfinite(events)):
        raise ValueError("Required event count is not finite")
    return events


def scan_exclusion(template_path, output_path, ratios=DEFAULT_EXCLUSION_RATIOS,
                   pion_name="piAntiMuon", kaon_name="KAntiMuon", threshold=9.0):
    pion, kaon, _ = load_templates(template_path, pion_name, kaon_name)
    events = required_events(pion, kaon, ratios, threshold)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(output_path) as output:
        output.mktree("sensitivityResult", {
            "ratio": np.asarray(ratios, dtype=float), "targetStatistics": events,
        })
        output["analysis_config"] = (
            f"input={Path(template_path).resolve()}; pion={pion_name}; kaon={kaon_name}; "
            f"threshold={threshold}; statistic=ROOT WW CHI2; solve=linear scaling"
        )
    return events
