def as_bool(v):
    if isinstance(v, str):
        # Include "1","yes" as redundant true values for convenience for future config file parsing
        return v.strip().lower() in ("true", "1", "yes")
    return bool(v)
