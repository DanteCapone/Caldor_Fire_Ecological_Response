"""Calculate MLTP Schmidt stability and Lake Number with Lake-Tools.

The two metric equations are imported from the vendored ``lake_params``
package. CTD potential density comes directly from ``mltp_ctd.pkl`` and TB4
wind supplies the wind forcing. Outputs are written to ``data/processed/ctd``.
"""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

from lake_tools_lake_params import lakenumber, schmidtstability, u_star


MIN_PROFILE_DEPTH_M = 400.0
WIND_HEIGHT_M = 10.0
WIND_DETECTION_FLOOR_M_S = 0.1
THERMOCLINE_SEARCH_MAX_M = 150.0


def trapz(y: np.ndarray, x: np.ndarray) -> float:
    """NumPy-version-independent trapezoidal integral."""
    fn = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    return float(fn(y, x))


def wind_for_cast(
    cast_time: pd.Timestamp,
    wind_time_ns: np.ndarray,
    wind_speed: np.ndarray,
    wind_direction: np.ndarray,
) -> dict[str, float | int | str | pd.Timestamp]:
    """Return the nearest raw TB4 wind observation without averaging."""
    if wind_time_ns.size == 0:
        return {
            "n_wind_observations": 0,
            "wind_window": "unavailable",
            "wind_observation_time": pd.NaT,
            "wind_time_offset_minutes": np.nan,
            "wind_speed_raw_m_s": np.nan,
            "wind_speed_forcing_m_s": np.nan,
            "wind_dir_deg": np.nan,
        }
    center = cast_time.value
    right = int(np.searchsorted(wind_time_ns, center, side="left"))
    candidates = [i for i in (right - 1, right) if 0 <= i < wind_time_ns.size]
    selected = min(candidates, key=lambda i: abs(int(wind_time_ns[i]) - center))
    selected_time = pd.to_datetime(int(wind_time_ns[selected]), utc=True)
    offset_minutes = abs((selected_time - cast_time).total_seconds()) / 60.0
    raw_speed = float(wind_speed[selected])
    forcing_speed = max(raw_speed, WIND_DETECTION_FLOOR_M_S)
    return {
        "n_wind_observations": 1,
        "wind_window": "nearest_raw_observation",
        "wind_observation_time": selected_time,
        "wind_time_offset_minutes": float(offset_minutes),
        "wind_speed_raw_m_s": raw_speed,
        "wind_speed_forcing_m_s": forcing_speed,
        "wind_dir_deg": float(wind_direction[selected]),
    }


def prepare_profile(
    source_depth: np.ndarray,
    source_density: np.ndarray,
    target_depth: np.ndarray,
) -> tuple[np.ndarray | None, float, int]:
    """Interpolate a cast to 1 m and extend endpoint densities to 0/500 m."""
    valid = np.isfinite(source_depth) & np.isfinite(source_density)
    if valid.sum() < 20:
        return None, np.nan, int(valid.sum())
    z = source_depth[valid]
    rho = source_density[valid]
    order = np.argsort(z)
    z = z[order]
    rho = rho[order]
    unique_z, unique_index = np.unique(z, return_index=True)
    rho = rho[unique_index]
    maximum_observed_depth = float(unique_z.max())
    if maximum_observed_depth < MIN_PROFILE_DEPTH_M:
        return None, maximum_observed_depth, int(unique_z.size)
    interpolated = np.interp(target_depth, unique_z, rho)
    return interpolated, maximum_observed_depth, int(unique_z.size)


def thermocline_depth(depth: np.ndarray, density: np.ndarray) -> tuple[float, float]:
    """Locate the maximum positive density gradient after 5 m smoothing."""
    smooth = (
        pd.Series(density)
        .rolling(window=5, center=True, min_periods=1)
        .mean()
        .to_numpy()
    )
    gradient = np.gradient(smooth, depth)
    eligible = (
        (depth >= 1.0)
        & (depth <= THERMOCLINE_SEARCH_MAX_M)
        & np.isfinite(gradient)
    )
    if not eligible.any():
        return np.nan, np.nan
    eligible_indices = np.flatnonzero(eligible)
    index = eligible_indices[int(np.argmax(gradient[eligible]))]
    return float(depth[index]), float(gradient[index])


def calculate_metrics(project_root: Path) -> pd.DataFrame:
    # Lake-Tools was authored against NumPy versions that exposed ``trapz``.
    # NumPy 2.4 removed that alias, so restore it without changing the equation.
    if not hasattr(np, "trapz"):
        np.trapz = np.trapezoid  # type: ignore[attr-defined]
    pkl_file = (
        project_root
        / "data/lake_environmental_data/ctd/01_ctd/mltp_ctd.pkl"
    )
    bathy_file = project_root / "data/processed/ctd/lake_tahoe_hypsography_usgs.csv"
    wind_file = project_root / "data/lake_environmental_data/00_met/TB4.csv"

    with pkl_file.open("rb") as handle:
        ctd = pickle.load(handle)
    source_depth = np.asarray(ctd["depth"], dtype=float)
    density_matrix = np.asarray(ctd["rho_gsw_p0"], dtype=float)
    cast_times = pd.to_datetime(ctd["dateutc"], utc=True)

    bathy = pd.read_csv(bathy_file, usecols=["depth_m", "area_m2"]).dropna()
    bathy = bathy.sort_values("depth_m").drop_duplicates("depth_m")
    target_depth = np.arange(0.0, 501.0, 1.0)
    area = np.interp(target_depth, bathy["depth_m"], bathy["area_m2"])
    surface_area = float(area[0])
    maximum_depth = float(target_depth[-1])
    volume = trapz(area, target_depth)
    center_volume_depth = trapz(target_depth * area, target_depth) / volume

    wind = pd.read_csv(
        wind_file,
        usecols=["dateutc", "WindSpeed", "WindDir"],
        parse_dates=["dateutc"],
    ).dropna(subset=["dateutc", "WindSpeed", "WindDir"])
    wind = wind.loc[
        wind["WindSpeed"].ge(0)
        & wind["WindDir"].between(0, 360, inclusive="both")
    ].sort_values("dateutc")
    wind_time_ns = pd.to_datetime(wind["dateutc"], utc=True).astype("int64").to_numpy()
    wind_speed = wind["WindSpeed"].to_numpy(dtype=float)
    wind_direction = wind["WindDir"].to_numpy(dtype=float)

    rows: list[dict[str, object]] = []
    for cast_index, cast_time in enumerate(cast_times):
        density, max_observed, n_observed = prepare_profile(
            source_depth, density_matrix[:, cast_index], target_depth
        )
        row: dict[str, object] = {
            "cast_id": cast_index + 1,
            "dateutc": cast_time,
            "date": cast_time.date(),
            "maximum_observed_depth_m": max_observed,
            "n_observed_depths": n_observed,
            "profile_endpoint_extension": "constant density to 0 and 500 m",
        }
        if density is None:
            row.update(
                schmidt_stability_j_m2=np.nan,
                thermocline_depth_m=np.nan,
                maximum_density_gradient_kg_m4=np.nan,
                epilimnion_density_kg_m3=np.nan,
                n_wind_observations=0,
                wind_window="not_calculated_profile_incomplete",
                wind_observation_time=pd.NaT,
                wind_time_offset_minutes=np.nan,
                wind_speed_raw_m_s=np.nan,
                wind_speed_forcing_m_s=np.nan,
                wind_dir_deg=np.nan,
                water_friction_velocity_m_s=np.nan,
                lake_number=np.nan,
            )
            rows.append(row)
            continue

        stability = float(
            schmidtstability(
                target_depth, density, area, surface_area, False, 0.0
            )
        )
        thermo_depth, max_gradient = thermocline_depth(target_depth, density)
        epi_mask = target_depth <= thermo_depth
        epi_volume = trapz(area[epi_mask], target_depth[epi_mask])
        epi_density = (
            trapz(density[epi_mask] * area[epi_mask], target_depth[epi_mask])
            / epi_volume
        )
        wind_values = wind_for_cast(
            cast_time, wind_time_ns, wind_speed, wind_direction
        )
        forcing_speed = float(wind_values["wind_speed_forcing_m_s"])
        u_star_value = (
            float(u_star(np.asarray([forcing_speed]), WIND_HEIGHT_M, epi_density)[0])
            if np.isfinite(forcing_speed)
            else np.nan
        )
        lake_number = (
            float(
                lakenumber(
                    stability,
                    np.asarray([u_star_value]),
                    surface_area,
                    epi_density,
                    maximum_depth,
                    thermo_depth,
                    center_volume_depth,
                )[0]
            )
            if np.isfinite(u_star_value) and u_star_value > 0 and stability > 0
            else np.nan
        )
        row.update(
            schmidt_stability_j_m2=stability,
            thermocline_depth_m=thermo_depth,
            maximum_density_gradient_kg_m4=max_gradient,
            epilimnion_density_kg_m3=epi_density,
            **wind_values,
            water_friction_velocity_m_s=u_star_value,
            lake_number=lake_number,
        )
        rows.append(row)

    return pd.DataFrame(rows)


def calculate_wind_resolution_diagnostic(
    metrics: pd.DataFrame, project_root: Path
) -> pd.DataFrame:
    """Estimate Lake Number at raw TB4 resolution for diagnostic plotting.

    Hydrographic terms are linearly interpolated between MLTP casts (with
    endpoint hold outside the first/last cast); wind is never averaged. This
    product is explicitly diagnostic and does not represent additional CTD
    observations.
    """
    focal = metrics.loc[
        pd.to_datetime(metrics["dateutc"], utc=True).dt.year.eq(2021)
        & metrics[
            [
                "schmidt_stability_j_m2",
                "thermocline_depth_m",
                "epilimnion_density_kg_m3",
            ]
        ].notna().all(axis=1)
    ].copy()
    focal = focal.sort_values("dateutc")
    cast_time = pd.to_datetime(focal["dateutc"], utc=True).astype("int64").to_numpy()

    wind = pd.read_csv(
        project_root / "data/lake_environmental_data/00_met/TB4.csv",
        usecols=["dateutc", "WindSpeed", "WindDir"],
        parse_dates=["dateutc"],
    ).dropna(subset=["dateutc", "WindSpeed", "WindDir"])
    wind["dateutc"] = pd.to_datetime(wind["dateutc"], utc=True)
    wind = wind.loc[
        wind["dateutc"].dt.year.eq(2021)
        & wind["WindSpeed"].ge(0)
        & wind["WindDir"].between(0, 360, inclusive="both")
    ].sort_values("dateutc")
    wind_time = wind["dateutc"].astype("int64").to_numpy()

    bathy = pd.read_csv(
        project_root / "data/processed/ctd/lake_tahoe_hypsography_usgs.csv",
        usecols=["depth_m", "area_m2"],
    ).dropna().sort_values("depth_m").drop_duplicates("depth_m")
    depth = np.arange(0.0, 501.0, 1.0)
    area = np.interp(depth, bathy["depth_m"], bathy["area_m2"])
    volume = trapz(area, depth)
    center_volume_depth = trapz(depth * area, depth) / volume
    surface_area = float(area[0])

    stability = np.interp(
        wind_time, cast_time, focal["schmidt_stability_j_m2"].to_numpy()
    )
    thermo = np.interp(
        wind_time, cast_time, focal["thermocline_depth_m"].to_numpy()
    )
    epi_density = np.interp(
        wind_time, cast_time, focal["epilimnion_density_kg_m3"].to_numpy()
    )
    wind_speed = wind["WindSpeed"].to_numpy(dtype=float)
    forcing_speed = np.maximum(wind_speed, WIND_DETECTION_FLOOR_M_S)
    friction = u_star(forcing_speed, WIND_HEIGHT_M, epi_density)
    lake = lakenumber(
        stability, friction, surface_area, epi_density, 500.0, thermo,
        center_volume_depth,
    )
    return pd.DataFrame(
        {
            "dateutc": wind["dateutc"].to_numpy(),
            "wind_speed_m_s": wind_speed,
            "wind_speed_forcing_m_s": forcing_speed,
            "wind_dir_deg": wind["WindDir"].to_numpy(dtype=float),
            "interpolated_schmidt_stability_j_m2": stability,
            "interpolated_thermocline_depth_m": thermo,
            "interpolated_epilimnion_density_kg_m3": epi_density,
            "water_friction_velocity_m_s": friction,
            "lake_number_diagnostic": lake,
            "hydrographic_temporal_method": (
                "linear interpolation between MLTP casts; endpoint hold"
            ),
        }
    )


def write_outputs(metrics: pd.DataFrame, project_root: Path) -> None:
    output_dir = project_root / "data/processed/ctd"
    output_dir.mkdir(parents=True, exist_ok=True)
    metrics_file = output_dir / "mltp_lake_tools_stability_metrics.csv"
    audit_file = output_dir / "mltp_lake_tools_stability_2021_audit.csv"
    methods_file = output_dir / "mltp_lake_tools_stability_methods.csv"
    wind_resolution_file = output_dir / "mltp_lake_number_raw_wind_2021.csv"
    metrics.to_csv(metrics_file, index=False)
    # Keep legacy filenames synchronized so downstream scripts cannot silently
    # fall back to the superseded rLakeAnalyzer calculation.
    metrics.to_csv(output_dir / "mltp_stability_metrics_from_pkl_tb4.csv", index=False)
    metric_dates = pd.to_datetime(metrics["date"])
    compatibility_schmidt = pd.DataFrame(
        {
            "CTD_ID": [f"lake_tools_cast_{cast_id:04d}" for cast_id in metrics["cast_id"]],
            "Event_ID": "",
            "date": metrics["date"],
            "maximum_depth_m": metrics["maximum_observed_depth_m"],
            "n_depths": metrics["n_observed_depths"],
            "schmidt_stability_j_m2": metrics["schmidt_stability_j_m2"],
            "deep_enough": metrics["maximum_observed_depth_m"].ge(MIN_PROFILE_DEPTH_M),
        }
    ).loc[metric_dates.dt.year.between(2005, 2025)]
    compatibility_schmidt.to_csv(
        output_dir / "mltp_schmidt_stability_2005_2025.csv", index=False
    )

    focal = metrics.loc[pd.to_datetime(metrics["date"]).dt.year.eq(2021)].copy()
    focal["schmidt_complete"] = np.isfinite(focal["schmidt_stability_j_m2"])
    focal["lake_number_complete"] = np.isfinite(focal["lake_number"])
    focal.to_csv(audit_file, index=False)
    calculate_wind_resolution_diagnostic(metrics, project_root).to_csv(
        wind_resolution_file, index=False
    )
    if focal.empty:
        raise RuntimeError("No MLTP casts were found in 2021.")
    if not focal["schmidt_complete"].all():
        missing = focal.loc[
            ~focal["schmidt_complete"],
            ["cast_id", "date", "schmidt_complete", "lake_number_complete"],
        ]
        raise RuntimeError(f"Incomplete 2021 Schmidt stability metrics:\n{missing}")
    if not focal["lake_number_complete"].all():
        missing = focal.loc[
            ~focal["lake_number_complete"],
            ["cast_id", "date", "wind_speed_raw_m_s", "wind_speed_forcing_m_s"],
        ]
        raise RuntimeError(f"Incomplete 2021 Lake Number metrics:\n{missing}")

    methods = pd.DataFrame(
        {
            "item": [
                "Metric implementation",
                "Source version",
                "CTD density",
                "Depth grid",
                "Profile inclusion",
                "Endpoint handling",
                "Thermocline",
                "Representative density",
                "Wind forcing",
                "Calm-wind handling",
                "Water friction velocity",
                "Lake Number completeness",
            ],
            "value": [
                "Vendored Lake-Tools schmidtstability, lakenumber, and u_star",
                "savalbuena/Lake-Tools commit f9cb02b516275b7f9ac2dc0687485c66d6fc56fc",
                "rho_gsw_p0 potential density from MLTP pickle",
                "1 m from 0 to 500 m with USGS Lake Tahoe hypsography",
                "Observed CTD depth must reach at least 400 m",
                "Linear interpolation within cast; nearest endpoint density extended to 0 and 500 m",
                "Maximum positive density gradient from a centered 5 m rolling-mean profile, searched 1-150 m",
                "Hypsography-weighted mean epilimnion density from 0 m to thermocline",
                "Nearest valid raw TB4 WindSpeed and WindDir observation; no temporal averaging",
                "Exact zero TB4 readings retained as observations but evaluated at the 0.1 m/s logger resolution floor; no temporal averaging",
                "Current Lake-Tools u_star at 10 m (Cd=0.001 below 5 m/s; otherwise 0.0013)",
                "Pipeline requires complete Schmidt stability and Lake Number for every 2021 cast; the documented 0.1 m/s calm-wind floor prevents undefined zero-wind values",
            ],
        }
    )
    methods.to_csv(methods_file, index=False)
    methods.to_csv(output_dir / "mltp_stability_methods.csv", index=False)
    print(f"Saved {metrics_file}")
    print(f"Saved {audit_file}")
    print(f"Saved {wind_resolution_file}")
    print(f"2021 casts: {len(focal)}; Lake Number missing: {focal['lake_number'].isna().sum()}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    root = args.project_root.resolve()
    write_outputs(calculate_metrics(root), root)


if __name__ == "__main__":
    main()
