"""dashboard.py -- turn a DashboardSpec into something the user sees.

Two halves: spec -> text (render_text, stdlib only) and spec ->
figure (render_matplotlib, Agg backend, headless). Both consume the
dashboard spec dict from results.build_dashboard, which is JSON-
serialisable so an HTTP API could hand it to a browser renderer.
"""
from __future__ import annotations

import os
import json
import logging

log = logging.getLogger(__name__)


def fmt_num(v, ndig=4):
    try:
        f = float(v)
    except Exception:
        return str(v)
    if f != 0 and (abs(f) < 1e-3 or abs(f) >= 1e4):
        return f"{f:.{ndig}g}"
    return f"{f:.{ndig}f}"


def series_stats(s):
    raw = s.get("value") or s.get("y") or []
    vals = [float(v) for v in raw if isinstance(v, (int, float))]
    if not vals:
        return 0.0, 0.0
    peak = max(abs(v) for v in vals)
    mean = sum(vals) / len(vals)
    return peak, mean


def render_text(spec, top_lines=8):
    out = []
    out.append("=" * 64)
    out.append("  AeroSolved run dashboard")
    out.append("  case        : " + str(spec.get("case", "?")))
    out.append("  postproc    : " + str(spec.get("postproc_path", "(none)")))
    out.append("=" * 64)

    scalars = spec.get("scalars", [])
    if scalars:
        out.append("")
        out.append("  Key scalars")
        out.append("    " + "-" * 54)
        for s in scalars:
            out.append("    %-30s %s %s" % (s["name"], fmt_num(s.get("value")), s.get("unit", "")))

    plots = spec.get("plots", [])
    if plots:
        out.append("")
        out.append("  Plots discovered: %d" % len(plots))
        for p in plots:
            out.append("    * " + p.title if False else "    * " + p["title"])
            series = p.get("series", [])
            for s in series[:4]:
                peak, mean = series_stats(s)
                out.append("      - %-40s  peak=%.3g  mean=%.3g" % (s.get("label",""), peak, mean))
            if len(series) > 4:
                out.append("      ... +%d more series" % (len(series) - 4))

    files = spec.get("files", [])
    if files:
        out.append("")
        out.append("  Data files    : %d" % len(files))
        for f in files[:top_lines]:
            out.append("    %-68s %s" % (f.get("path","?"), f.get("shape", [])))
        if len(files) > top_lines:
            out.append("    ... +%d more" % (len(files) - top_lines))

    warns = spec.get("warnings", [])
    if warns:
        out.append("")
        out.append("  Warnings")
        out.append("    " + "-" * 54)
        for w in warns:
            out.append("    ! %s" % w)
    return "\n".join(out)


def render_matplotlib(spec, out_dir, fmt="pdf"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(out_dir, exist_ok=True)
    saved = {}
    for i, p in enumerate(spec.get("plots", [])):
        title = p.get("title", "plot_%d" % i)
        fig, ax = plt.subplots()
        xscale = p.get("xscale")
        yscale = p.get("yscale")
        for s in p.get("series", []):
            label = s.get("label", "")
            if "diameter" in s and "value" in s:
                x, y = s.get("diameter", []), s.get("value", [])
                if xscale == "log" and yscale == "log":
                    ax.loglog(x, y, marker="o", ms=3, linewidth=1, label=label)
                elif xscale == "log":
                    ax.semilogx(x, y, marker="o", ms=3, label=label)
                else:
                    ax.plot(x, y, marker="o", ms=3, label=label)
            elif "category" in s and "value" in s:
                cats = s.get("category", [])
                left = cats.index(s["label"]) if s["label"] in cats else 0
                width = 0.8 / max(1, len(cats))
                ax.bar(left, s.get("value", []), width, label=label)
                ax.set_xticks(range(len(cats)))
                ax.set_xticklabels(cats)
            elif "x" in s and "y" in s:
                ax.plot(s["x"], s["y"], marker="o", ms=3, label=label)
        ax.set_xlabel(p.get("x_axis", "x"))
        ax.set_ylabel(p.get("y_axis", "y"))
        ax.set_title(title)
        if any("label" in s for s in p.get("series", [])):
            ax.legend(loc="best", fontsize=8)
        ax.grid(True, which="both" if xscale == "log" else "major", alpha=0.4, linewidth=0.5)
        fig.tight_layout()
        safe = "".join(c if (c.isalnum() or c in "-_., ") else "_" for c in title).strip()
        path = os.path.join(out_dir, "%s.%s" % (safe or ("plot_%d" % i), fmt))
        fig.savefig(path, dpi=120, bbox_inches="tight")
        saved[title] = path
        plt.close(fig)
    return saved


def spec_to_png(spec, out_dir, fmt="png"):
    return render_matplotlib(spec, out_dir, fmt=fmt)


def spec_json(spec, out_path):
    with open(out_path, "w") as f:
        json.dump(spec, f, indent=2, default=str)
    return out_path


def render(spec, kind="text", out_dir="."):
    kind = kind.lower()
    if kind == "text":
        return render_text(spec)
    if kind in ("matplotlib", "pdf"):
        return "\n".join(render_matplotlib(spec, out_dir, fmt=kind).values())
    if kind == "png":
        return "\n".join(spec_to_png(spec, out_dir).values())
    if kind == "json":
        return spec_json(spec, os.path.join(out_dir, "dashboard.json"))
    raise ValueError("unknown render kind: %r" % kind)
