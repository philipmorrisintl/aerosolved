"""case_config.py - discover cases, parse Allrun VARS blocks.

Given an aerosolved repo root, this module enumerates parameterized cases
(cases that have an `Allrun` script), and for each one extracts:

 * the solver name (aerosolEulerFoam, aerosolBuoyantEulerFoam, ...) read from
   the case's `system/controlDict` / `system/controlDict.m4`;
 * the **VARS defaults** - the dict of `-DVARNAME=$SHELLVAR` pairs that the
   Allrun feeds to `m4`, where each `SHELLVAR` has a default value the user
   can override;
 * the **positional CLI args** the Allrun takes (e.g. <mesh> and
   <sectional|moment>);
 * the set of **model choices** the Allrun selects between.

All of this is derived by static analysis of the Allrun + case sources, so no
case needs a hand-written entry to be UI-enabled; a new case that follows the
AeroSolved conventions just works.
"""

from __future__ import annotations

import os
import re
import dataclasses
import logging
from typing import Optional

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Regexes
# ---------------------------------------------------------------------------

# The `VARS="..."` shell variable that holds the list of -D flags. It is a
# multi-line, backslash-continued assignment, so it needs re.DOTALL and
# `^...$` in MULTILINE context.
VARS_BLOCK_RE = re.compile(
    r'^VARS=(?P<q>[\'"])(?P<body>.*?)(?P=q)\s*$',
    re.MULTILINE | re.DOTALL,
)

# Individual `-DNAME=$VAR` tokens inside the VARS block body.
DNAMES_IN_VARS = re.compile(r'-D(?P<name>\w+)=\$(?P<var>\w+)')

# `X1=-5` style shell assignments (the "sane defaults" the UI will prefill).
SHELL_ASSIGN_RE = re.compile(r'^\s*(?P<name>[A-Za-z_]\w*)\s*=\s*(?P<val>[^\\\n]+)')

# Positional args the Allrun takes: `$1`, `$2`, ...
POSARG_RE = re.compile(r'\s*\$([1-9][0-9]?)\b')

# The solver, read from `application <name>;` in a controlDict.
APPLICATION_RE = re.compile(r'^\s*application\s+(\w+)\s*;', re.MULTILINE)


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclasses.dataclass
class ShellVar:
    """A single `$VAR` the UI exposes as a knob."""
    name: str                       # shell var name (e.g. 'X1', 'Y3')
    m4_name: str                    # -D name (e.g. 'VARX1')
    default: str                    # default expression, as a string
    help: str = ""                  # optional human-readable hint
    kind: str = "number"            # "number" | "string" | "enum"
    choices: Optional[list[str]] = None

    def __post_init__(self):
        if self.kind == "number":
            try:
                float(self.default)
            except (ValueError, TypeError):
                self.kind = "raw"


@dataclasses.dataclass
class CaseDef:
    """Configuration of a single UI-enabled case."""
    name: str
    path: str
    allrun_path: str
    plot_py: Optional[str]
    model_var_name: str = "MODEL"
    model_choices: list[str] = dataclasses.field(default_factory=list)
    positional_args: list[str] = dataclasses.field(default_factory=list)
    application: Optional[str] = None
    vars: list[ShellVar] = dataclasses.field(default_factory=list)
    vars_by_m4_name: dict[str, ShellVar] = dataclasses.field(default_factory=dict)
    raw_vars_block: str = ""
    result_shape: str = "unknown"

    def __str__(self):
        return f"<CaseDef {self.name!r} vars={len(self.vars)} args={self.positional_args}>"

    def knob(self, m4_name: str) -> Optional[ShellVar]:
        return self.vars_by_m4_name.get(m4_name)

    def to_vars_dict(self) -> dict[str, str]:
        return {v.name: v.default for v in self.vars}


# ---------------------------------------------------------------------------

def find_repo_root(start: str) -> Optional[str]:
    """Walk up the filesystem to find the aerosolved repo root
    (a dir that has `Allwmake` and a `cases/` subdir)."""
    cur = os.path.abspath(start)
    for _ in range(10):
        if os.path.isfile(os.path.join(cur, "Allwmake")) and os.path.isdir(os.path.join(cur, "cases")):
            return cur
        parent = os.path.dirname(cur)
        if parent == cur:
            break
        cur = parent
    return None


def dedupe_keep_order(seq: list[str]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for s in seq:
        if s not in seen:
            out.append(s)
            seen.add(s)
    return out


def _scan_allrun(allrun_text: str) -> tuple[list[ShellVar], list[str], str]:
    """Static analysis of an Allrun source.

    Returns (shell_vars, positional_args, raw_vars_block).
     """
    raw_vars_block = ""
    m = VARS_BLOCK_RE.search(allrun_text)
    if m:
       raw_vars_block = m.group("body")

    # 1. Collect all simple shell assignments (the knob defaults).
    shell_assigns: dict[str, str] = {}
    for line in allrun_text.splitlines():
        mm = SHELL_ASSIGN_RE.match(line)
        if mm:
            vname = mm.group("name")
            val = mm.group("val").strip().rstrip("\\").strip()
            shell_assigns[vname] = val

    # 2. Map each `-DNAME=$VAR` in the VARS block to its default.
    shell_vars: list[ShellVar] = []
    seen_m4: set[str] = set()
    if raw_vars_block:
        for mm in DNAMES_IN_VARS.finditer(raw_vars_block):
            m4_name = mm.group("name")
            var_name = mm.group("var")
            if m4_name in seen_m4:
                continue
            seen_m4.add(m4_name)
            default = shell_assigns.get(var_name, "")
            shell_vars.append(ShellVar(name=var_name, m4_name=m4_name, default=default))

    # 3. Positional args. Heuristically classify each discovered `$N`.
    positional_args: list[str] = []
    seen_args: set[int] = set()
    for line in allrun_text.splitlines():
        for am in POSARG_RE.finditer(line):
            n = int(am.group(1))
            if n in seen_args:
                continue
            seen_args.add(n)
            kind = "arg%d" % n
            if re.search(r'case\s+\$\{%d\}' % n, allrun_text):
                kind = "model"
            elif re.search(r'case\s+\$\{%d\}' % n, allrun_text) or \
                 re.search(r'Invalid\s+mesh', allrun_text, re.I):
                kind = "mesh"
            positional_args.append(kind)
    positional_args = dedupe_keep_order(positional_args)

    return shell_vars, positional_args, raw_vars_block


def _infer_model_choices(allrun_text: str) -> list[str]:
    """Find the enum values a model-type arg takes.

    Detects `case $N in  moment) ... sectional) ...  esac` blocks.
    """
    choices: set[str] = set()
    for block in re.finditer(r'case\s*\$(\d+)\s+in\n(.*?)esac', allrun_text, re.DOTALL):
        body = block.group(2)
        for case_m in re.finditer(r'(\w+)\)\s+(.*?)(?=\n\s*\w+\)\s+|\n\s*esac)',
                                 body, re.DOTALL):
            name = case_m.group(1)
            if re.search(r'INVALID|invalid', case_m.group(2), re.I):
                continue
            if name not in ("*") and len(name) > 1:
                choices.add(name)
    return sorted(choices)


def _read_application(case_path: str) -> Optional[str]:
    """Read the solver name from the case's controlDict (or controlDict.m4)."""
    for cand in ("system/controlDict.m4", "system/controlDict"):
        p = os.path.join(case_path, cand)
        if os.path.isfile(p):
            m = APPLICATION_RE.search(open(p).read())
            if m:
                return m.group(1)
    return None


# ---------------------------------------------------------------------------

def load_case(repo_root: str, case_name: str) -> CaseDef:
    """Load and analyse one case by name from a repo root."""
    case_path = os.path.join(repo_root, "cases", case_name)
    allrun_path = os.path.join(case_path, "Allrun")
    if not os.path.isfile(allrun_path):
        raise FileNotFoundError(f"No Allrun script under {case_path!r}")
    allrun_text = open(allrun_path).read()

    shell_vars, positional_args, raw_vars_block = _scan_allrun(allrun_text)
    model_choices = _infer_model_choices(allrun_text)
    application = _read_application(case_path)

    plot_py = os.path.join(case_path, "plot.py")
    plot_py = plot_py if os.path.isfile(plot_py) else None

    return CaseDef(
        name=case_name,
        path=case_path,
        allrun_path=allrun_path,
        plot_py=plot_py,
        model_choices=model_choices,
        positional_args=positional_args,
        application=application or "aerosolEulerFoam",
        vars=shell_vars,
        vars_by_m4_name={v.m4_name: v for v in shell_vars},
        raw_vars_block=raw_vars_block,
        result_shape=case_name,
    )


def discover_cases(repo_root: str) -> list[CaseDef]:
    """Enumerate every parameterized case (one with an `Allrun`) under cases/."""
    cases_dir = os.path.join(repo_root, "cases")
    if not os.path.isdir(cases_dir):
        return []
    out: list[CaseDef] = []
    for name in sorted(os.listdir(cases_dir)):
        allrun = os.path.join(cases_dir, name, "Allrun")
        if not os.path.isfile(allrun):
            continue
        try:
            out.append(load_case(repo_root, name))
        except Exception as e:
            log.warning("case %r skipped: %s", name, e)
    return out
