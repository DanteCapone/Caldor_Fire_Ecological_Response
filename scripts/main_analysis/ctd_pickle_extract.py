"""Extract compact, analysis-ready CTD products from the LTP/MLTP pickle files.

The source pickles store one depth x cast matrix per variable.  This script
keeps the pickle files as the authoritative input while avoiding a very large
long-form intermediate table.  It writes profile summaries for plotting,
1 m Mid-lake temperature profiles for stability calculations, and LTP surface
temperature means for Figure 2 v2.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


FOCAL_START = pd.Timestamp("2021-06-01", tz="UTC")
FOCAL_END_EXCLUSIVE = pd.Timestamp("2022-05-01", tz="UTC")
MONTH_ORDER = [6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4]
MIN_CLIM_CASTS = 3

VARIABLES = {
    "T": "Temperature",
    "rho_gsw_p0": "Potential_Density",
    "cond": "Conductivity",
    "do": "Dissolved_Oxygen",
    "chla": "Chl_Fluorescence",
    "turbidity": "Turbidity",
    "par": "PAR",
}

STATIONS = {
    "ltp": {"label": "Index", "depth_cap": 160},
    "mltp": {"label": "Mid-lake", "depth_cap": 500},
}


def bin_to_one_metre(depth: np.ndarray, values: np.ndarray, cap: int) -> np.ndarray:
    """Average the native 0.2 m grid into 1 m bins, preserving cast columns."""
    rounded = np.rint(depth).astype(int)
    output = np.full((cap + 1, values.shape[1]), np.nan, dtype=float)
    with np.errstate(invalid="ignore"):
        for z in range(cap + 1):
            rows = rounded == z
            if rows.any():
                block = values[rows, :]
                count = np.isfinite(block).sum(axis=0)
                total = np.nansum(block, axis=0)
                output[z, :] = np.divide(
                    total,
                    count,
                    out=np.full(values.shape[1], np.nan, dtype=float),
                    where=count > 0,
                )
    return output


def smooth_profiles(values: np.ndarray) -> np.ndarray:
    """Apply the established centered 5 m moving median to each cast."""
    return (
        pd.DataFrame(values)
        .rolling(window=5, center=True, min_periods=1)
        .median()
        .to_numpy()
    )


def profile_rows(
    values: np.ndarray,
    dates: pd.Series,
    station: str,
    variable: str,
) -> list[dict]:
    rows: list[dict] = []
    years = dates.dt.year.to_numpy()
    months = dates.dt.month.to_numpy()
    focal = (dates >= FOCAL_START) & (dates < FOCAL_END_EXCLUSIVE)

    for month in MONTH_ORDER:
        target_year = 2021 if month >= 6 else 2022
        obs_mask = focal.to_numpy() & (months == month)
        clim_mask = (months == month) & (years != target_year)

        obs = values[:, obs_mask]
        clim = values[:, clim_mask]
        for depth_m in range(values.shape[0]):
            obs_values = obs[depth_m, :] if obs.shape[1] else np.array([])
            clim_values = clim[depth_m, :] if clim.shape[1] else np.array([])
            obs_values = obs_values[np.isfinite(obs_values)]
            clim_values = clim_values[np.isfinite(clim_values)]

            obs_n = int(obs_values.size)
            clim_n = int(clim_values.size)
            if obs_n == 0 and clim_n < MIN_CLIM_CASTS:
                continue

            rows.append(
                {
                    "Station_ID": station,
                    "variable": variable,
                    "month": month,
                    "month_label": pd.Timestamp(2001, month, 1).strftime("%b"),
                    "depth_m": depth_m,
                    "obs_mean": np.mean(obs_values) if obs_n else np.nan,
                    "obs_sd": np.std(obs_values, ddof=1) if obs_n > 1 else np.nan,
                    "obs_n_casts": obs_n,
                    "clim_mean": np.mean(clim_values)
                    if clim_n >= MIN_CLIM_CASTS
                    else np.nan,
                    "clim_sd": np.std(clim_values, ddof=1)
                    if clim_n >= MIN_CLIM_CASTS and clim_n > 1
                    else np.nan,
                    "clim_n_casts": clim_n,
                }
            )
    return rows


def main(project_root: Path) -> None:
    source_dir = (
        project_root
        / "data"
        / "lake_environmental_data"
        / "ctd"
        / "01_ctd"
    )
    output_dir = project_root / "data" / "processed" / "ctd"
    output_dir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict] = []
    availability_rows: list[dict] = []
    mltp_temperature: pd.DataFrame | None = None
    ltp_surface_rows: list[dict] = []

    for station_key, station_config in STATIONS.items():
        source = source_dir / f"{station_key}_ctd.pkl"
        if not source.exists():
            raise FileNotFoundError(f"Missing CTD pickle: {source}")

        print(f"Reading {source.name} ...", flush=True)
        payload = pd.read_pickle(source)
        dates = pd.to_datetime(payload["dateutc"], utc=True)
        depth = np.asarray(payload["depth"], dtype=float)
        cap = int(station_config["depth_cap"])

        for source_name, display_name in VARIABLES.items():
            matrix = np.asarray(payload[source_name], dtype=float)
            one_metre = smooth_profiles(bin_to_one_metre(depth, matrix, cap))
            summary_rows.extend(
                profile_rows(
                    one_metre,
                    dates,
                    str(station_config["label"]),
                    display_name,
                )
            )

            years = dates.dt.year.to_numpy()
            months = dates.dt.month.to_numpy()
            focal = ((dates >= FOCAL_START) & (dates < FOCAL_END_EXCLUSIVE)).to_numpy()
            for month in MONTH_ORDER:
                target_year = 2021 if month >= 6 else 2022
                obs_mask = focal & (months == month)
                clim_mask = (months == month) & (years != target_year)
                availability_rows.append(
                    {
                        "Station_ID": station_config["label"],
                        "variable": display_name,
                        "month": month,
                        "month_label": pd.Timestamp(2001, month, 1).strftime("%b"),
                        "focal_casts": int(np.isfinite(one_metre[:, obs_mask]).any(axis=0).sum()),
                        "climatology_casts": int(
                            np.isfinite(one_metre[:, clim_mask]).any(axis=0).sum()
                        ),
                    }
                )

            if station_key == "mltp" and source_name == "T":
                mltp_temperature = pd.DataFrame(one_metre)
                mltp_temperature.columns = [f"cast_{i + 1:04d}" for i in range(len(dates))]
                mltp_temperature.insert(0, "depth_m", np.arange(cap + 1))
                mltp_temperature = mltp_temperature.melt(
                    id_vars="depth_m", var_name="cast_id", value_name="temperature_c"
                )
                cast_lookup = pd.DataFrame(
                    {
                        "cast_id": [f"cast_{i + 1:04d}" for i in range(len(dates))],
                        "dateutc": dates.dt.strftime("%Y-%m-%dT%H:%M:%SZ"),
                        "date": dates.dt.strftime("%Y-%m-%d"),
                    }
                )
                mltp_temperature = mltp_temperature.merge(cast_lookup, on="cast_id")
                mltp_temperature = mltp_temperature[
                    np.isfinite(mltp_temperature["temperature_c"])
                ][["cast_id", "dateutc", "date", "depth_m", "temperature_c"]]

            if station_key == "ltp" and source_name == "T":
                surface = one_metre[:11, :]
                counts = np.isfinite(surface).sum(axis=0)
                means = np.divide(
                    np.nansum(surface, axis=0),
                    counts,
                    out=np.full(surface.shape[1], np.nan),
                    where=counts > 0,
                )
                ltp_surface_rows = [
                    {
                        "cast_id": f"cast_{i + 1:04d}",
                        "dateutc": dates.iloc[i].strftime("%Y-%m-%dT%H:%M:%SZ"),
                        "date": dates.iloc[i].strftime("%Y-%m-%d"),
                        "surface_temperature_c": means[i],
                        "n_depth_bins_0_10m": int(counts[i]),
                    }
                    for i in range(len(dates))
                    if np.isfinite(means[i])
                ]

    summary = pd.DataFrame(summary_rows)
    summary["obs_lo"] = summary["obs_mean"] - summary["obs_sd"].fillna(0)
    summary["obs_hi"] = summary["obs_mean"] + summary["obs_sd"].fillna(0)
    summary["clim_lo"] = summary["clim_mean"] - summary["clim_sd"].fillna(0)
    summary["clim_hi"] = summary["clim_mean"] + summary["clim_sd"].fillna(0)
    summary.to_csv(
        output_dir / "ctd_profile_summary_jun2021_apr2022_from_pkl.csv", index=False
    )
    pd.DataFrame(availability_rows).to_csv(
        output_dir / "ctd_profile_availability_jun2021_apr2022_from_pkl.csv",
        index=False,
    )
    if mltp_temperature is None:
        raise RuntimeError("MLTP temperature matrix was not extracted")
    mltp_temperature.to_csv(
        output_dir / "mltp_temperature_profiles_1m_from_pkl.csv", index=False
    )
    pd.DataFrame(ltp_surface_rows).to_csv(
        output_dir / "ltp_surface_temperature_0_10m_from_pkl.csv", index=False
    )
    print(f"Wrote CTD products to {output_dir}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    args = parser.parse_args()
    main(args.project_root.resolve())
