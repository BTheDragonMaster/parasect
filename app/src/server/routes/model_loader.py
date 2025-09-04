import threading, logging
from dataclasses import dataclass
from typing import Callable, Dict, Optional, Any


try:
    import psutil
    _proc = psutil.Process()
except Exception:  # psutil optional
    _proc = None 


# Prefer joblig for large sklearn models (supports mmap arrays)
try:
    from joblib import load as joblib_load
except Exception:
    joblib_load = None


@dataclass
class ModelSpec:
    name: str
    path: str
    loader: Optional[Callable[[str], Any]] = None  # custom loader
    mmap: bool = True  # if True and joblib available, mmap arrays (lower peak RAM)


class MultiModelLoader:
    def __init__(self, specs: Dict[str, ModelSpec]):
        self._specs = specs
        self._models: Dict[str, Any] = {}
        self._locks: Dict[str, threading.Lock] = {k: threading.Lock() for k in specs}

    def _rss(self) -> str:
        """Return the resident set size (RSS) memory usage of the process."""
        if not _proc:
            return "n/a"
        return f"{_proc.memory_info().rss/1024/1024:.1f} MB"
    
    def _default_load(self, path: str, mmap: bool) -> Any:
        # Try joblib with mmap first
        if joblib_load:
            mmap_mode = "r" if mmap else None
            return joblib_load(path, mmap_mode=mmap_mode)
        # Fallback to pickle
        import pickle
        with open(path, "rb") as fh:
            return pickle.load(fh)
        
    def get(self, key: str):
        if key in self._models:
            return self._models[key]
        if key not in self._specs:
            raise KeyError(f"Unknown model key: {key}")

        spec = self._specs[key]
        lock = self._locks[key]

        with lock:
            if key in self._models:  # double check
                return self._models[key]
            logging.info(f"Loading model '%s' from %s (rss=%s)", spec.name, spec.path, self._rss())
            loader = spec.loader or (lambda p: self._default_load(p, spec.mmap))
            model = loader(spec.path)
            self._models[key] = model
            logging.info(f"Loaded model '%s' from (rss=%s)", spec.name, self._rss())
            return model
