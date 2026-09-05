"""
Central on/off switch for this package's tqdm progress bars, independent of
logging configuration -- for non-interactive/batch use (e.g. an automated
experiment harness) where progress bars are noise in captured output, not
signal, and remembering to configure every module's logger individually
(see execution.py's progress_logger) is error-prone: setting a parent
logger to DEBUG re-enables any child logger that has no level of its own,
bringing bars back unintentionally.

Defaults to enabled (bars show), matching interactive use. Call
disable_progress_bars() once, explicitly, at the top of a script/harness to
turn every progress bar in this package off in one step -- this always wins
over any per-module logger configuration. Deliberately not
environment-variable-driven: config should be visible in the code that
calls it, not implicit in whatever happens to be set in the shell.
"""
_disabled = False


def disable_progress_bars() -> None:
    global _disabled
    _disabled = True


def enable_progress_bars() -> None:
    global _disabled
    _disabled = False


def progress_bars_disabled() -> bool:
    return _disabled
