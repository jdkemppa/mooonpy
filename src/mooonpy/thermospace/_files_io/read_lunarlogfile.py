from typing import Union
from mooonpy import Path

def read_lunarlogfile(file: Union[Path, str], keywords: Union[str, list]) -> dict:
    keywords_set = set([keywords] if isinstance(keywords, str) else keywords)

    log_data = {}
    with open(file) as f:
        for line in f:
            key, _, rest = line.partition(":")
            key = key.strip()
            if key in keywords_set:
                value = float(rest.split()[0])
                log_data[key] = value

    return log_data