"""
Release hygiene script -- a handful of checks/actions that used to be
hand-cranked (and hand-remembered) before every release: regenerating
requirements.txt, confirming pyproject.toml's version and CHANGELOG.md's
top entry agree, doing a clean build + twine check + full test run, and
flagging untracked files sitting in the repo that might get missed.

This performs those checks/actions directly, but deliberately does NOT
perform the release itself (no twine upload, no git tag/push, no git
commit) -- those are irreversible/external actions that stay a human
decision. If every check passes, it prints the exact commands for the
remaining steps instead of running them.

None of this belongs in the shipped public API -- it lives here (rather
than tests/ or a top-level scripts/ dir) because it's maintainer tooling
about the package itself, not usage of it.

Run from the project root with the project's own venv, e.g.:
    skip/Scripts/python.exe -m skipalignments.util.release_check

Exits 0 if every check passes, 1 otherwise -- suitable for a pre-release
step, not (yet) wired into CI.
"""
import os
import re
import subprocess
import sys
from pathlib import Path

# src/skipalignments/util/release_check.py -> repo root is three levels up
PROJECT_ROOT = Path(__file__).resolve().parents[3]

# Mirrors what `source skip/Scripts/activate` puts on PATH: the venv's own
# Scripts dir first (where console-script shims and any project-local
# executable, e.g. a hardlinked ebi.exe, live). Without this, a bare
# `skip/Scripts/python.exe -m ...` invocation (no activation) can't
# resolve those by bare name -- notably, some tests' ProcessPoolExecutor
# workers re-import probabilities.py fresh and fall back to its default
# EBI_EXECUTABLE='ebi', which only resolves via PATH.
_SUBPROCESS_ENV = {**os.environ, "PATH": str(Path(sys.executable).parent) + os.pathsep + os.environ.get("PATH", "")}


def _run(cmd, cwd=PROJECT_ROOT):
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, env=_SUBPROCESS_ENV)


def _pyproject_field(name: str) -> str:
    text = (PROJECT_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(rf'^{name}\s*=\s*"([^"]+)"', text, re.MULTILINE)
    if match is None:
        raise ValueError(f"pyproject.toml: could not find a top-level {name} = \"...\" line")
    return match.group(1)


def regenerate_requirements() -> str:
    """
    Overwrites requirements.txt with the current venv's `pip freeze`,
    stripping this package's own self-reference line (pip freeze emits a
    `-e git+...#egg=<name>` entry for whatever's installed editable in the
    same venv it's run from -- that's the package being released, not a
    dependency of it, so it doesn't belong in a dependency-environment
    record).
    """
    package_name = _pyproject_field("name")
    result = _run([sys.executable, "-m", "pip", "freeze"])
    if result.returncode != 0:
        raise RuntimeError(f"pip freeze failed:\n{result.stderr}")

    self_reference = re.compile(rf"#egg={re.escape(package_name)}\b")
    lines = [line for line in result.stdout.splitlines() if not self_reference.search(line)]

    requirements_path = PROJECT_ROOT / "requirements.txt"
    requirements_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return f"requirements.txt regenerated ({len(lines)} entries, self-reference stripped)"


def check_version_consistency() -> tuple[bool, str]:
    """
    pyproject.toml's version must match CHANGELOG.md's topmost entry, and
    that entry must be a real dated release, not a lingering [Unreleased]
    (or missing a date -- Keep a Changelog's convention this project
    follows, e.g. `## [0.2.2] - 2026-09-10`).
    """
    pyproject_version = _pyproject_field("version")

    changelog_text = (PROJECT_ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    heading = re.search(r"^## \[([^\]]+)\](?:\s*-\s*(\d{4}-\d{2}-\d{2}))?", changelog_text, re.MULTILINE)
    if heading is None:
        return False, "CHANGELOG.md: no '## [version] - date' heading found at all"

    changelog_version, changelog_date = heading.group(1), heading.group(2)

    if changelog_version == "Unreleased":
        return False, (
            f"CHANGELOG.md's top entry is still [Unreleased] -- pyproject.toml is at "
            f"{pyproject_version}; move the release notes under a dated "
            f"'## [{pyproject_version}] - YYYY-MM-DD' heading before releasing"
        )
    if changelog_date is None:
        return False, f"CHANGELOG.md's top entry ([{changelog_version}]) has no date"
    if changelog_version != pyproject_version:
        return False, (
            f"version mismatch: pyproject.toml says {pyproject_version}, "
            f"CHANGELOG.md's top entry says {changelog_version}"
        )
    return True, f"pyproject.toml and CHANGELOG.md agree: {pyproject_version} ({changelog_date})"


def build_check_and_test() -> tuple[bool, str]:
    """Clean rebuild, twine check the result, run the full test suite."""
    for stale in ("dist", "build"):
        stale_path = PROJECT_ROOT / stale
        if stale_path.is_dir():
            import shutil
            shutil.rmtree(stale_path)
    for egg_info in (PROJECT_ROOT / "src").glob("*.egg-info"):
        import shutil
        shutil.rmtree(egg_info)

    build_result = _run([sys.executable, "-m", "build"])
    if build_result.returncode != 0:
        return False, f"python -m build failed:\n{build_result.stdout}\n{build_result.stderr}"

    dist_files = list((PROJECT_ROOT / "dist").glob("*"))
    if not dist_files:
        return False, "python -m build produced no artifacts in dist/"

    twine_result = _run([sys.executable, "-m", "twine", "check", *[str(f) for f in dist_files]])
    if twine_result.returncode != 0:
        return False, f"twine check failed:\n{twine_result.stdout}\n{twine_result.stderr}"

    test_result = _run([sys.executable, "-m", "unittest", "discover", "-s", "tests"])
    if test_result.returncode != 0:
        # unittest's own summary is on stderr
        tail = "\n".join(test_result.stderr.strip().splitlines()[-15:])
        return False, f"test suite failed:\n{tail}"

    return True, f"build + twine check + full test suite all passed ({len(dist_files)} dist artifacts)"


def check_stray_files() -> tuple[bool, str]:
    """
    Untracked files are easy to either forget (they never make it into a
    commit and quietly vanish) or accidentally sweep in (a careless
    `git add -A`). Neither is what you want right before a release --
    flag them so it's a decision, not an accident.
    """
    result = _run(["git", "status", "--porcelain"])
    if result.returncode != 0:
        return False, f"git status failed:\n{result.stderr}"

    untracked = [line[3:] for line in result.stdout.splitlines() if line.startswith("??")]
    if untracked:
        return False, "untracked files present: " + ", ".join(untracked)
    return True, "no untracked files"


def _next_steps() -> str:
    version = _pyproject_field("version")
    python = sys.executable
    return "\n".join([
        "All checks passed. Remaining steps are yours to run:",
        "",
        f"  {python} -m twine upload dist/*",
        "  git add -A",
        "  git commit -m \"Release v" + version + "\"",
        f"  git tag v{version}",
        "  git push && git push --tags",
    ])


def main() -> int:
    ok = True

    # requirements.txt regeneration and the stray-files check both act on
    # git-visible state, so run stray-files last: regenerating
    # requirements.txt only ever *modifies* a tracked file (never adds an
    # untracked one), and build_check_and_test's dist/build/*.egg-info are
    # all gitignored, so neither changes what counts as "stray" -- but
    # checking last still means the report reflects the repo exactly as
    # this run leaves it.
    print(regenerate_requirements())

    for check in (check_version_consistency, build_check_and_test, check_stray_files):
        passed, message = check()
        print(("PASS: " if passed else "FAIL: ") + message)
        ok = ok and passed

    print()
    print(_next_steps() if ok else "One or more checks failed -- not printing release commands.")

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
