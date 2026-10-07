"""Bind manuscript claims to results, compile, and enforce page/layout checks."""

import argparse
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

from pypdf import PdfReader


def sci(x):
    a, b = f"{x:.3e}".split("e")
    return rf"{a}\times10^{{{int(b)}}}"


def numbers(root):
    data = json.loads((root / "examples/results/metrics.json").read_text())
    methods = data["pulse"]["methods"]
    reduction = methods["cn"]["final_time_error"] / methods["rannacher"]["final_time_error"]
    commands = {
        "CNMinimum": f"{methods['cn']['first_min']:.6f}",
        "CNEnergy": f"{methods['cn']['first_energy_ratio']:.6f}",
        "StepSize": sci(data["pulse"]["dt"]),
        "ExpmError": sci(data["matrix_exponential_max_error"]),
        "MaxDrift": sci(max(m["relative_heat_drift"] for m in methods.values())),
        "Reduction": f"{reduction:.1f}",
        "CNOrder": f"{data['temporal_last_orders']['cn']:.5f}",
        "BEOrder": f"{data['temporal_last_orders']['be']:.5f}",
        "RanOrder": f"{data['temporal_last_orders']['rannacher']:.5f}",
    }
    commands["PulseRows"] = "\n".join(
        f"{m.upper() if m != 'rannacher' else 'Startup + CN'} & "
        f"${sci(v['first_min'])}$ & {v['first_energy_ratio']:.5f} & "
        f"${sci(v['final_time_error'])}$ " + r"\\"
        for m, v in methods.items()
    )
    commands["SpaceRows"] = "\n".join(
        f"{r['cells']} & ${sci(r['cell_average_l2_error'])}$ & "
        + ("--" if i == 0 else f"{data['spatial_orders'][i - 1]:.5f}")
        + r"\\"
        for i, r in enumerate(data["spatial"])
    )
    return "\n".join(rf"\newcommand{{\{k}}}{{{v}}}" for k, v in commands.items()) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", help="path to tectonic or pdflatex")
    args = parser.parse_args()
    folder = Path(__file__).resolve().parent
    root = folder.parent
    (folder / "numbers.tex").write_text(numbers(root), encoding="utf-8", newline="\n")
    engine = args.engine or shutil.which("tectonic") or shutil.which("pdflatex")
    if not engine:
        raise SystemExit("Install Tectonic or TeX Live and pass --engine if necessary")
    if Path(engine).exists():
        engine = str(Path(engine).resolve())
    name = Path(engine).stem.lower()
    build = folder / "build"
    build.mkdir(exist_ok=True)
    if name == "tectonic":
        command = [
            engine,
            "--keep-logs",
            "--keep-intermediates",
            "--outdir",
            str(build),
            "paper.tex",
        ]
        passes = 1
    elif name == "pdflatex":
        command = [
            engine,
            "-no-shell-escape",
            "-halt-on-error",
            "-interaction=nonstopmode",
            "-output-directory",
            str(build),
            "paper.tex",
        ]
        passes = 2
    else:
        raise SystemExit("Supported engines: tectonic, pdflatex")
    for _ in range(passes):
        result = subprocess.run(
            command, cwd=folder, capture_output=True, text=True, errors="replace", timeout=180
        )
        if result.returncode:
            raise SystemExit(result.stdout + result.stderr)
    log = (build / "paper.log").read_text(errors="replace")
    issues = re.findall(
        r"Overfull \\[hv]box[^\n]*|[^\n]*(?:undefined references|"
        r"Citation .* undefined|Reference .* undefined)[^\n]*",
        log,
    )
    if issues:
        raise SystemExit("Layout/reference errors:\n" + "\n".join(issues))
    pdf = build / "paper.pdf"
    pages = PdfReader(pdf).pages
    if not 1 <= len(pages) <= 5 or any(not (p.extract_text() or "").strip() for p in pages):
        raise SystemExit(f"Page count/empty-page check failed: {len(pages)}")
    inputs = [
        folder / "paper.tex",
        folder / "numbers.tex",
        folder / "build.py",
        root / "pyproject.toml",
        root / "README.md",
        root / "MANIFEST.in",
        root / ".github/workflows/validate.yml",
    ]
    for directory, glob in [
        ("src", "*.py"),
        ("tests", "*.py"),
        ("scripts", "*.py"),
        ("docs", "*.md"),
        ("examples/results", "*"),
    ]:
        inputs += sorted((root / directory).rglob(glob))
    report = {
        "pages": len(pages),
        "page_limit": 5,
        "engine": name,
        "overfull_boxes": 0,
        "undefined_references": 0,
        "inputs_sha256": {
            p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in inputs
            if p.is_file()
        },
        "pdf_sha256": hashlib.sha256(pdf.read_bytes()).hexdigest(),
    }
    shutil.copyfile(pdf, folder / "paper.pdf")
    (folder / "validation.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n"
    )
    print(f"Compiled {len(pages)} pages; no undefined references or overfull boxes.")


if __name__ == "__main__":
    main()
