"""Severity thresholds for the SDevice output log."""

FILE_LOG_LEVELS = ("debug", "info", "warning", "error")


def log_enabled(threshold, category, message=""):
    if category in ("HEADER", "SUMMARY"):
        return True
    if category == "FATAL":
        severity = "error"
    elif category == "WARNING":
        severity = "warning"
    elif category == "RUNTIME":
        # Native runtime messages start with a timestamp and the spdlog severity.
        tag = message.partition("] [")[2].partition("]")[0]
        severity = {"warning": "warning", "error": "error", "critical": "error"}.get(tag, "debug")
    elif category == "RUNNER":
        severity = "debug"
    else:
        severity = "info"
    return FILE_LOG_LEVELS.index(severity) >= FILE_LOG_LEVELS.index(threshold)
