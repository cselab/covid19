from pathlib import Path
import functools
import json
import pickle

def cache(func):
    """Caches the result of a given function, depending on the function arguments.

    A decorated function accepts only hashable types, e.g. tuples and not lists.

    Example:
        @cache
        def func(x):
            print(x)
            return x * x

        a = func(10)    # Prints 10.
        b = func(20)    # Prints 20.
        c = func(20)    # Prints nothing (value is cached).
        print(a, b, c)  # 100, 400, 400
    """
    _cache = {}
    def inner(*args, **kwargs):
        key = (args, frozenset(kwargs.items()))
        if key in _cache:
            return _cache[key]
        _cache[key] = out = func(*args, **kwargs)
        return out
    return functools.wraps(func)(inner)


def cache_to_file(target, dependencies=[]):
    """Factory for a decorator that caches the result of a no-argument function and stores them in a target file.

    Arguments:
        target: The target cache filename.
        dependencies: The list of files that the result depends on.

    If the target file exists but is older than any of the dependencies, it will be recomputed.
    """
    target = Path(target)
    dependencies = [Path(d) for d in dependencies]

    if target.suffix == '.json':
        def load(path):
            with open(path) as f:
                return json.load(f)

        def save(content, path):
            with open(path, 'w') as f:
                json.dump(content, f)

    elif target.suffix == '.pickle':
        def load(path):
            with open(path, 'rb') as f:
                return pickle.load(f)

        def save(content, path):
            with open(path, 'wb') as f:
                pickle.dump(content, f)

    else:
        raise ValueError(f"Unrecognized extension '{target.suffix}'. "
                         f"Only .json and .pickle supported.")

    def decorator(func):
        def inner():
            try:
                modified_time = target.stat().st_mtime
            except FileNotFoundError:
                pass
            else:
                if all(modified_time < d.stat().st_mtime
                       for d in dependencies):
                    return load(target)

            result = func()
            target.parent.mkdir(parents=True, exist_ok=True)
            save(result, target)
            return result
        return functools.wraps(func)(inner)
    return decorator
