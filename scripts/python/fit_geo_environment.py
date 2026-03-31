import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def clamp(x, lo, hi):
    return max(lo, min(hi, x))


def rmse(y_true, y_pred):
    if len(y_true) == 0 or len(y_true) != len(y_pred):
        return 0.0
    d = np.asarray(y_true) - np.asarray(y_pred)
    return float(np.sqrt(np.mean(d * d)))


def parse_floats_from_line(line):
    return [float(x) for x in re.findall(r"[+-]?\d+\.\d+E[+-]\d+", line)]


def parse_geo_dat(path):
    lines = Path(path).read_text(encoding="utf-8", errors="ignore").splitlines()
    envs = []

    i = 0
    while i < len(lines):
        m = re.match(r"\s*No\.\s*(\d+)\s+(.*)", lines[i])
        if not m:
            i += 1
            continue

        env_id = int(m.group(1))
        env_name = m.group(2).strip()

        electron_energy = None
        electron_flux = None

        j = i + 1
        while j < len(lines):
            if re.match(r"\s*No\.\s*\d+\s+", lines[j]):
                break

            if "Electron energy(eV):" in lines[j]:
                electron_energy = parse_floats_from_line(lines[j])
            elif "Flux(1/m^2.s.sr.eV):" in lines[j] and electron_energy is not None and electron_flux is None:
                electron_flux = parse_floats_from_line(lines[j])

            j += 1

        if electron_energy is not None and electron_flux is not None:
            n = min(len(electron_energy), len(electron_flux))
            envs.append(
                {
                    "id": env_id,
                    "name": env_name,
                    "energy": np.array(electron_energy[:n], dtype=float),
                    "flux": np.array(electron_flux[:n], dtype=float),
                }
            )

        i = j

    return envs


def fit_three_models(energy_ev, intensity):
    observed_energy_ev = np.asarray(energy_ev, dtype=float)
    observed_intensity = np.asarray(intensity, dtype=float)
    n = observed_energy_ev.size

    if n < 5:
        return {
            "ok": False,
            "message": "观测能谱数据不足，保留当前分布",
        }

    model_double = np.zeros(n, dtype=float)
    model_kappa = np.zeros(n, dtype=float)
    model_power = np.zeros(n, dtype=float)

    sum_w = 0.0
    mean_e = 0.0
    for i in range(n):
        w = max(0.0, float(observed_intensity[i]))
        sum_w += w
        mean_e += float(observed_energy_ev[i]) * w
    mean_e = mean_e / sum_w if sum_w > 0.0 else 10.0

    e_split = mean_e
    cold_sum = 0.0
    hot_sum = 0.0
    cold_w = 0.0
    hot_w = 0.0

    for i in range(n):
        e = max(1e-9, float(observed_energy_ev[i]))
        w = max(0.0, float(observed_intensity[i]))
        if e <= e_split:
            cold_sum += e * w
            cold_w += w
        else:
            hot_sum += e * w
            hot_w += w

    t1 = max(1e-6, cold_sum / cold_w) if cold_w > 0.0 else max(1e-6, mean_e * 0.5)
    t2 = max(t1, hot_sum / hot_w) if hot_w > 0.0 else max(t1, mean_e * 1.5)
    hot_fraction = clamp(hot_w / sum_w, 0.02, 0.98) if sum_w > 0.0 else 0.2

    sx = 0.0
    sy = 0.0
    sxx = 0.0
    sxy = 0.0
    m = 0.0

    for i in range(n // 2, n):
        x = np.log(max(1e-9, float(observed_energy_ev[i])))
        y = np.log(max(1e-12, float(observed_intensity[i])))
        sx += x
        sy += y
        sxx += x * x
        sxy += x * y
        m += 1.0

    slope = -2.5
    if m > 2.0:
        denom = m * sxx - sx * sx
        if abs(denom) > 1e-12:
            slope = (m * sxy - sx * sy) / denom

    fitted_power = clamp(-slope, 1.1, 8.0)
    fitted_kappa = clamp(1.5 + 0.8 * fitted_power, 1.6, 12.0)

    for i in range(n):
        e = max(1e-9, float(observed_energy_ev[i]))
        model_double[i] = (1.0 - hot_fraction) * np.exp(-e / t1) + hot_fraction * np.exp(-e / t2)

        theta = max(1e-6, mean_e / max(2.0, fitted_kappa))
        model_kappa[i] = np.power(1.0 + e / (fitted_kappa * theta), -(fitted_kappa + 1.0))

        model_power[i] = np.power(e + 1e-3, -fitted_power)

    norm_obs = float(np.sum(observed_intensity))
    norm_obs = norm_obs if norm_obs > 0.0 else 1.0

    norm_d = float(np.sum(model_double))
    norm_k = float(np.sum(model_kappa))
    norm_p = float(np.sum(model_power))

    model_double = model_double * norm_obs / max(1e-12, norm_d)
    model_kappa = model_kappa * norm_obs / max(1e-12, norm_k)
    model_power = model_power * norm_obs / max(1e-12, norm_p)

    score_double = rmse(observed_intensity, model_double)
    score_kappa = rmse(observed_intensity, model_kappa)
    score_power = rmse(observed_intensity, model_power)

    selected = "DOUBLE_MAXWELL"
    best = score_double
    if score_kappa < best:
        best = score_kappa
        selected = "KAPPA"
    if score_power < best:
        selected = "POWER_LAW"

    return {
        "ok": True,
        "mean_e": mean_e,
        "t1": t1,
        "t2": t2,
        "hot_fraction": hot_fraction,
        "fitted_power": fitted_power,
        "fitted_kappa": fitted_kappa,
        "score_double": score_double,
        "score_kappa": score_kappa,
        "score_power": score_power,
        "selected": selected,
        "model_double": model_double,
        "model_kappa": model_kappa,
        "model_power": model_power,
    }


def main():
    root = Path(__file__).resolve().parents[2]
    dat_path = root / "demos" / "GEO_Environment.dat"
    out_dir = root / "results" / "geo_env_fit"
    out_dir.mkdir(parents=True, exist_ok=True)

    envs = parse_geo_dat(dat_path)
    if len(envs) != 5:
        print(f"[WARN] 识别到环境数量: {len(envs)}，预期 5")

    summary_rows = []

    for env in envs:
        fit = fit_three_models(env["energy"], env["flux"])
        if not fit["ok"]:
            print(f"[SKIP] Env {env['id']} 数据不足")
            continue

        summary_rows.append(
            {
                "env_id": env["id"],
                "env_name": env["name"],
                "rmse_double_maxwell": fit["score_double"],
                "rmse_kappa": fit["score_kappa"],
                "rmse_power_law": fit["score_power"],
                "selected": fit["selected"],
                "mean_e": fit["mean_e"],
                "t1": fit["t1"],
                "t2": fit["t2"],
                "hot_fraction": fit["hot_fraction"],
                "kappa": fit["fitted_kappa"],
                "power_index": fit["fitted_power"],
            }
        )

        x = env["energy"]
        y = np.asarray(env["flux"], dtype=float)

        fig, ax = plt.subplots(figsize=(9, 6))
        ax.plot(x, y, "o", ms=4, label="Observed", color="#1f77b4")
        ax.plot(x, fit["model_double"], "-", lw=2, label=f"Double Maxwell (RMSE={fit['score_double']:.3e})", color="#d62728")
        ax.plot(x, fit["model_kappa"], "-", lw=2, label=f"Kappa (RMSE={fit['score_kappa']:.3e})", color="#2ca02c")
        ax.plot(x, fit["model_power"], "-", lw=2, label=f"Power Law (RMSE={fit['score_power']:.3e})", color="#ff7f0e")

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Energy (eV)")
        ax.set_ylabel("Flux (1/m^2.s.sr.eV)")
        ax.set_title(f"Env {env['id']} - {env['name']}\nBest: {fit['selected']}")
        ax.grid(True, which="both", ls="--", alpha=0.35)
        ax.legend(fontsize=8)

        out_png = out_dir / f"env_{env['id']}_fit.png"
        fig.tight_layout()
        fig.savefig(out_png, dpi=160)
        plt.close(fig)
        print(f"[OK] 图已生成: {out_png}")

    summary_path = out_dir / "fit_summary.csv"
    with summary_path.open("w", encoding="utf-8") as f:
        f.write(
            "env_id,env_name,rmse_double_maxwell,rmse_kappa,rmse_power_law,selected,mean_e,t1,t2,hot_fraction,kappa,power_index\n"
        )
        for r in summary_rows:
            f.write(
                f"{r['env_id']},{r['env_name']},{r['rmse_double_maxwell']:.10e},{r['rmse_kappa']:.10e},{r['rmse_power_law']:.10e},{r['selected']},{r['mean_e']:.10e},{r['t1']:.10e},{r['t2']:.10e},{r['hot_fraction']:.10e},{r['kappa']:.10e},{r['power_index']:.10e}\n"
            )

    print(f"[DONE] 汇总: {summary_path}")


if __name__ == "__main__":
    main()
