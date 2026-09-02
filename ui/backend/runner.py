"""runner.py -- orchestrate an AeroSolved case run and stream events.

A "run" mirrors the command-line workflow, with knob overrides applied
automatically:

 1. Copy the source case to an isolated run directory.
 2. Patch the copy's `Allrun`: knob shell-var defaults are rewritten to the
    values the user set in the UI.
 3. Execute `./Allrun <mesh> <model>`, streaming stdout+stderr line by line.
 4. Track phases (mesh -> fields -> solve -> postProcess) from log markers and
    emit structured events the frontend displays.
 5. On success, locate the postProcessing directory for the results layer.

The runner is display-independent: it yields JSON-serialisable events, so it may
be driven by a Qt worker thread or a headless test.
"""

from __future__ import annotations

import os
import re
import glob
import shutil
import time
import dataclasses
import logging
import subprocess
from typing import Iterator, Optional

from backend import m4_expand

log = logging.getLogger(__name__)


@dataclasses.dataclass
class RunEvent:
    type: str                      # log | phase | progress | error | done
    phase: str = ""
    message: str = ""
    detail: str = ""
    ts: float = dataclasses.field(default_factory=time.time)

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)


_PHASE_MARKERS = [
    ("mesh",        [re.compile(r"blockMesh", re.I)]),
    ("fields",      [re.compile(r"setFields|mapFields", re.I)]),
    ("decompose",   [re.compile(r"decomposePar", re.I)]),
    ("solve",       [re.compile(r"Courant Number|deltaT|Iteration", re.I)]),
    ("postProcess", [re.compile(r"postProcess|sample|Writing.*plot", re.I)]),
]

PHASE_ORDER = ["mesh", "fields", "decompose", "solve", "postProcess"]
_PHASE_PROGRESS = {"mesh": 0.1, "fields": 0.3, "decompose": 0.4, "solve": 0.6, "postProcess": 0.9}


def detect_phase(line: str) -> Optional[str]:
    for phase, pats in _PHASE_MARKERS:
        for p in pats:
            if p.search(line):
                return phase
    return None


_KNOB_ASSIGN_RE = re.compile(
    r"^(?P<indent>\s*)(?P<name>[A-Za-z_]\w*)\s*=\s*(?P<val>[^\\\n]+)"
)


def patch_allrun_knobs(allrun_text: str, knobs: dict[str, str]) -> str:
    """Rewrite shell-var defaults in `allrun_text` to `knobs` values.

    Only lines that are simple `NAME = value` assignments are touched, so loop
    bodies, conditionals, etc. are left intact (conservative).
    """
    if not knobs:
        return allrun_text
    out: list[str] = []
    for line in allrun_text.splitlines(keepends=True):
        m = _KNOB_ASSIGN_RE.match(line)
        if m and m.group("name") in knobs:
            trail = line[m.end():]
            out.append("%s%s=%s%s" % (m.group("indent"), m.group("name"),
                                      knobs[m.group("name")], trail))
        else:
            out.append(line)
    return "".join(out)


@dataclasses.dataclass
class RunConfig:
    case_path: str
    run_dir: str
    model: str = "sectional"
    mesh_cells: Optional[str] = None
    knobs: dict[str, str] = dataclasses.field(default_factory=dict)
    timeout: Optional[int] = None
    env: dict[str, str] = dataclasses.field(default_factory=dict)
    raw_vars: Optional[str] = None



def _kill(proc: subprocess.Popen) -> None:
    """Kill the process group (Allrun spawns blockMesh/solver children)."""
    try:
        os.killpg(os.getpgid(proc.pid), 9)
    except Exception:
        try:
            proc.kill()
        except Exception:
            pass


class RunManager:
    PHASES = PHASE_ORDER

    def __init__(self, cfg: RunConfig):
        self.cfg = cfg

    def prepare(self) -> None:
        dst = self.cfg.run_dir
        allrun_dst = os.path.join(dst, "Allrun")
        if os.path.isfile(allrun_dst):
            for pat in ("Time", "processor0", "log.*"):
                for p in glob.glob(os.path.join(dst, pat)):
                    shutil.rmtree(p, ignore_errors=True) if os.path.isdir(p) else os.remove(p)
            text = open(allrun_dst).read()
        else:
            shutil.copytree(self.cfg.case_path, dst, symlinks=True,
                            ignore=shutil.ignore_patterns(".git", "Time"))
            text = open(allrun_dst).read()

        if self.cfg.knobs:
            text = patch_allrun_knobs(text, self.cfg.knobs)
            with open(allrun_dst, "w") as fh:
                fh.write(text)
        # Pre-expand *.m4 templates (parameterized cases) with knob values so
        # the solver stage sees resolved VARS.  Idempotent: MacroExpander only
        # touches *.m4 files that still exist, and Allrun's own setMacros is a
        # no-op afterwards.
        m4_expand.MacroExpander(
            self.cfg.run_dir,
            knobs=self.cfg.knobs or {},
            raw_vvars_override=self.cfg.raw_vars,
        ).pre_expand_inplace()

    def command(self) -> list[str]:
        allrun = os.path.join(self.cfg.run_dir, "Allrun")
        args = [allrun]
        if self.cfg.mesh_cells is not None:
            args += [str(self.cfg.mesh_cells), self.cfg.model]
        elif self.cfg.model:
            args += [self.cfg.model]
        return args

    def stream(self) -> Iterator[RunEvent]:
        self.prepare()
        cmd = self.command()
        yield RunEvent(type="phase", phase="mesh", message="preparing run directory")

        env = os.environ.copy()
        for k, v in self.cfg.env.items():
            env[str(k)] = str(v)

        seen: set[str] = set()
        try:
            proc = subprocess.Popen(
                cmd, cwd=self.cfg.run_dir, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, text=True, bufsize=1, env=env,
                start_new_session=True,
            )
        except Exception as e:
            yield RunEvent(type="error", message=str(e), detail="spawn")
            yield RunEvent(type="done", phase="error", message=str(e), detail="spawn")
            return

        start = time.time()
        assert proc.stdout is not None
        for line in proc.stdout:
            if self.cfg.timeout and (time.time() - start) > self.cfg.timeout:
                _kill(proc)
                yield RunEvent(type="error",
                               message="run timed out after %ds" % self.cfg.timeout,
                               detail="timeout")
                yield RunEvent(type="done", phase="error", message="timeout", detail="timeout")
                return
            yield RunEvent(type="log", message=line.rstrip("\n"))
            ph = detect_phase(line)
            if ph and ph not in seen:
                seen.add(ph)
                yield RunEvent(type="phase", phase=ph, message="ph-> %s" % ph)
                yield RunEvent(type="progress", phase=ph,
                               detail=str(_PHASE_PROGRESS.get(ph, 0.5)))

        proc.wait()
        rc = proc.returncode
        if rc == 0:
            yield RunEvent(type="done", phase="success", message="run complete", detail="0")
        else:
            yield RunEvent(type="done", phase="error",
                           message="run exited non-zero", detail=str(rc))

    def results_dir(self) -> Optional[str]:
        base = self.cfg.run_dir
        for cand in (os.path.join(base, "postProcessing"),
                     os.path.join(base, "processor0", "postProcessing")):
            if os.path.isdir(cand):
                return cand
        return None
