
from typing import Any, Dict, List
import yaml


def _parse_value(value: Any) -> Any:
    """
    Convert strings like `'1,16,1'` → `[1, 16, 1]`
    Preserves floats, `None`, etc.
    """
    if isinstance(value, str) and "," in value:
        parts = value.split(",")
        try:
            out = []
            for x in parts:
                x = x.strip()
                if x.isdigit():
                    out.append(int(x))
                else:
                    out.append(float(x))
            return out
        except ValueError:
            return parts
    return value


def loadSidecarYaml(fPath: str, keys: List[str]) -> Dict[str, Any]:
    """
    Reads only the specified literal keys (including dots).
    """
    with open(fPath, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)

    result = {}
    for key in keys:
        value = data.get(key, None)
        value = _parse_value(value)
        result[key] = value

    return result
