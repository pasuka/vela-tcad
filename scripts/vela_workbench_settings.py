"""Read and update per-user SDS preferences."""
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import tempfile

from vela_sdevice.log_levels import FILE_LOG_LEVELS


@dataclass(frozen=True)
class Settings:
    language: str = "en"
    font_size: int = 10
    log_level: str = "debug"


def load_settings(path):
    try:
        return read_settings(path)[0]
    except FileNotFoundError:
        defaults = Settings()
        write_settings(path, asdict(defaults))
        return defaults


def read_settings(path):
    values = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(values, dict):
        raise ValueError("Expected a JSON object")
    settings = Settings(**{key: values[key] for key in asdict(Settings()) if key in values})
    if settings.language not in ("en", "zh-CN"):
        raise ValueError('language must be "en" or "zh-CN"')
    if type(settings.font_size) is not int or not 9 <= settings.font_size <= 16:
        raise ValueError("font_size must be an integer between 9 and 16")
    if settings.log_level not in FILE_LOG_LEVELS:
        raise ValueError("log_level must be debug, info, warning or error")
    return settings, values


def save_language(path, language):
    _, values = read_settings(path)
    values["language"] = language
    write_settings(path, values)


def write_settings(path, values):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name)
        try:
            json.dump(values, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
            stream.close()
            temporary.replace(path)
        finally:
            temporary.unlink(missing_ok=True)
