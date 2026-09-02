"""m4_expand.py -- expand AeroSolved's `*.m4` case templates.

Mirrors `setMacros` in `scripts/AeroSolvedRunFunctions`:

    setMacros() {
        find . -name '*.m4' | while read IN; do
            OUT=$(echo $IN | rev | cut -c 4- | rev)
            m4 $1 $IN > $OUT
        done
    }

i.e. for every `*.m4` in the case tree, run `m4 -DNAME=VALUE ... $IN` and write
the result to the same path minus the trailing `.m4`. This is how OpenFOAM
dictionaries that reference `$VARNx` get resolved at run time.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import logging

log = logging.getLogger(__name__)


def m4_available() -> bool:
    return shutil.which("m4") is not None


def shell_quote(s: str) -> str:
    """Single-quote for the shell; escapes inner quotes."""
    return "'" + s.replace("'", r"'\''") + "'"


class MacroExpander:
    """Given `knobs` (dict of shell var name -> string value), expand the
    `*.m4` templates in a case directory tree."""

    def __init__(self, case_path: str, knobs=None,
                 raw_vvars_override: str | None = None,
                 dry_run: bool = False):
        self.case_path = case_path
        self.knobs = knobs or {}
        self.raw_vvars_override = raw_vvars_override
        self.dry_run = dry_run

    def m4_args(self) -> str:
        """Return the `-DVAR=VALUE` argument string that `m4` receives.

        Two sources of `-D` flags:
          1) from `self.knobs` (UI-driven; knob values take precedence);
          2) from `self.raw_vvars_override` (advanced mode: a pasted
              `VARS="..."` string; only names NOT in the knobs are added,
           so an explicit knob always wins over the raw block).
        """
        parts: list[str] = []

        # 1. Knob-driven -D flags.
        for name, val in self.knobs.items():
            parts.append("-D%s=%s" % (name, str(val)))

        # 2. Raw VARS override: parse -DNAME=$VAR; add names absent from the
        #    knob set (knobs win on collision).
        if self.raw_vvars_override:
            knob_names = set(self.knobs)
            for mm in re.finditer(r"-D(\w+)=\$(\w+)", self.raw_vvars_override):
                name, var = mm.group(1), mm.group(2)
                if name in knob_names:
                    continue
                val = os.environ.get(var, "")
                parts.append("-D%s=%s" % (name, val))

        return " ".join(parts)

    def _find_m4_templates(self) -> list[str]:
        out: list[str] = []
        for root, dirs, files in os.walk(self.case_path):
            dirs[:] = [d for d in dirs if d not in (".git", "Time")]
            for fn in files:
                if fn.endswith(".m4"):
                    out.append(os.path.join(root, fn))
        return sorted(out)

    def _out_path(self, m4_path: str) -> str:
        assert m4_path.endswith(".m4"), m4_path
        return m4_path[:-3]

    def preflight_check(self) -> None:
        """Raise if `m4` is unavailable but there are templates to expand."""
        if not self._find_m4_templates():
            return
        if not m4_available():
            raise RuntimeError(
                "m4 is required to expand the *.m4 templates in this case but is "
                "not in PATH. Install `m4` (e.g. `sudo apt-get install m4`) or use "
                "a pre-expanded case." % self.case_path)

    def expand_into(self, out_root: str) -> list[str]:
        """Write each resolved template into `out_root/<rel path>` (minus `.m4`).
        Returns the list of output paths created."""
        templates = self._find_m4_templates()
        if not templates:
            return []
        if self.dry_run:
            return [self._out_path(t) for t in templates]
        if not m4_available():
            raise RuntimeError("`m4` not found; cannot expand %d template(s)" % len(templates))

        m4_args = self.m4_args()
        out_paths: list[str] = []
        for tp in templates:
            rel_out = os.path.relpath(self._out_path(tp), self.case_path)
            out_full = os.path.join(out_root, rel_out)
            os.makedirs(os.path.dirname(out_full) or ".", exist_ok=True)
            cmd = "m4 %s %s > %s" % (
                shell_quote(m4_args) if m4_args else "",
                shell_quote(tp),
                shell_quote(out_full),
            )
            r = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
            if r.returncode != 0:
                raise RuntimeError(
                    "m4 expansion of %s failed:\n  stdout: %s\n  stderr: %s"
                    % (tp, r.stdout[-500:], r.stderr[-500:]))
            out_paths.append(out_full)
        return out_paths

    def expand_inplace(self) -> list[str]:
        """Expand templates in the case directory itself (writes next to each .m4)."""
        return self.expand_into(self.case_path)

    def pre_expand_inplace(self) -> list[str]:
        """Expand in-place only for cases that still need it.

        Skips expansion when `m4` is unavailable.  Returns the list of
        output files created (empty when there is nothing to do).
        """
        if not self._find_m4_templates():
            return []
        if not m4_available():
            n = len(self._find_m4_templates())
            log.warning("m4 missing -- skipping %d *.m4 template(s) in %s"
                            % (n, self.case_path))
            return []
        return self.expand_inplace()
