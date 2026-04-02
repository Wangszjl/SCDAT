#!/usr/bin/env python3
import argparse
import math
import re
from pathlib import Path


def parse_numeric_tokens(text):
    return [float(token) for token in re.findall(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?", text)]


def smoothline3(x, y):
    n = min(len(x), len(y))
    if n <= 2:
        return y[:n]
    ys = [0.0] * n
    for i in range(n):
        if i == 0:
            first = 0
        elif i >= n - 1:
            first = n - 3
        else:
            first = i - 1
        idx = [first, first + 1, first + 2]
        t0, t1, t2 = x[idx[0]], x[idx[1]], x[idx[2]]
        p0, p1, p2 = y[idx[0]], y[idx[1]], y[idx[2]]
        det = (t0 - t1) * (t0 - t2) * (t1 - t2)
        if abs(det) < 1.0e-30:
            ys[i] = y[i]
            continue
        a0 = (p0 * t1 * t2 * (t1 - t2) + p1 * t2 * t0 * (t2 - t0) + p2 * t0 * t1 * (t0 - t1)) / det
        a1 = (p0 * (t2 * t2 - t1 * t1) + p1 * (t0 * t0 - t2 * t2) + p2 * (t1 * t1 - t0 * t0)) / det
        a2 = (p0 * (t1 - t2) + p1 * (t2 - t0) + p2 * (t0 - t1)) / det
        ys[i] = a0 + a1 * x[i] + a2 * x[i] * x[i]
    return ys


def insert(x, y, target):
    n = min(len(x), len(y))
    if n == 0:
        return 0.0
    if n == 1:
        return y[0]
    if n == 2:
        if x[0] == x[1]:
            return y[0]
        return (y[0] * (target - x[1]) - y[1] * (target - x[0])) / (x[0] - x[1])
    if target <= x[1]:
        k, m = 0, 2
    elif target >= x[n - 2]:
        k, m = n - 3, n - 1
    else:
        k, m = 1, n
        while m - k != 1:
            i = (k + m) // 2
            if target < x[i - 1]:
                m = i
            else:
                k = i
        k -= 1
        m -= 1
        if abs(target - x[k]) < abs(target - x[m]):
            k -= 1
        else:
            m += 1
    value = 0.0
    for i in range(k, m + 1):
        weight = 1.0
        for j in range(k, m + 1):
            if j != i:
                weight *= (target - x[j]) / (x[i] - x[j])
        value += weight * y[i]
    return value


def parse_input_file(path):
    values = [float(line.strip()) for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(values) < 12:
        raise ValueError("input_LEO_Charging.dat does not contain enough values")
    return {
        "jsun": int(round(values[0])),
        "jram": int(round(values[1])),
        "h_km": values[2],
        "flow_angle_deg": values[3],
        "jenv": int(round(values[4])),
    }


def parse_environment_case(path, case_id):
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    case_header = f"No. {case_id}"
    for idx, line in enumerate(lines):
        if line.strip().startswith(case_header):
            thermal = parse_numeric_tokens(lines[idx + 4])
            electron_energy = parse_numeric_tokens(lines[idx + 6])
            electron_flux = parse_numeric_tokens(lines[idx + 7])
            ion_energy = parse_numeric_tokens(lines[idx + 8])
            ion_flux = parse_numeric_tokens(lines[idx + 9])
            return thermal, electron_energy, electron_flux, ion_energy, ion_flux
    raise ValueError(f"Environment case {case_id} not found in {path}")


def reconstruct_body_je0(thermal, electron_energy, electron_flux, ion_energy):
    ne = [thermal[0], thermal[2], thermal[4]]
    te = [thermal[1], thermal[3], thermal[5]]
    ti = [thermal[7], thermal[10], thermal[13]]
    ee = [value for value in electron_energy if value > 0.0]
    ei = [value for value in ion_energy if value > 0.0]
    fluxe = [value for value in electron_flux[: len(ee)]]
    fluxes = smoothline3(ee, fluxe)

    emin = min(ee[0], ei[0])
    emax = max(ee[-1], ei[-1])
    for k in range(3):
        if te[k] != 0.0:
            emin = min(emin, 0.05 * min(te[k], ti[k]))
            emax = max(emax, 10.0 * max(te[k], ti[k]))

    num = 1000
    nmin = math.log10(emin)
    nmax = math.log10(emax)
    dn = (nmax - nmin) / num
    xe = 0.0

    for i in range(num):
        energy = 10 ** (nmin + i * dn)
        d_energy = 10 ** (nmin + (i + 1) * dn) - 10 ** (nmin + i * dn)
        if energy < ee[0]:
            fluxe0 = sum(
                5.325e10 * ne[k] * energy * (te[k] ** -1.5) * math.exp(-energy / te[k])
                if te[k] != 0.0
                else 0.0
                for k in range(3)
            )
        elif energy <= ee[-1]:
            fluxe0 = max(0.0, insert(ee, fluxes, energy))
        else:
            fluxe0 = sum(
                5.325e10 * ne[k] * energy * (te[k] ** -1.5) * math.exp(-energy / te[k])
                if te[k] != 0.0
                else 0.0
                for k in range(3)
            )
        xe += fluxe0 * d_energy

    je_na_per_m2 = -1.0e9 * math.pi * 1.60217733e-19 * xe
    return je_na_per_m2


def main():
    parser = argparse.ArgumentParser(
        description="Reconstruct the legacy LEO body Je(0) from input/environment files."
    )
    parser.add_argument(
        "--input",
        default="ref/LEO/input_LEO_Charging.dat",
        help="Path to input_LEO_Charging.dat",
    )
    parser.add_argument(
        "--environment",
        default="ref/LEO/LEO_Environment.dat",
        help="Path to LEO_Environment.dat",
    )
    args = parser.parse_args()

    input_config = parse_input_file(args.input)
    thermal, electron_energy, electron_flux, ion_energy, _ = parse_environment_case(
        args.environment, input_config["jenv"]
    )
    je0_na_per_m2 = reconstruct_body_je0(thermal, electron_energy, electron_flux, ion_energy)
    print(f"input={Path(args.input).resolve()}")
    print(f"environment={Path(args.environment).resolve()}")
    print(f"jenv={input_config['jenv']}")
    print(f"legacy_leo_body_je0_na_per_m2={je0_na_per_m2:.12e}")
    print(f"legacy_leo_body_je0_a_per_m2={je0_na_per_m2 * 1.0e-9:.12e}")


if __name__ == "__main__":
    main()
