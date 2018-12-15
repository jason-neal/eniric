import yaml
import os


class Config:
    """
    This is a class to allow lazy-loading of config files. This means that the atributes are only read from the file
    when accessed in the code, allowing for changes in the file without having to restart the python instance
    """

    def __init__(self, filename):
        self.filename = filename
        # Format string for saving/reading orders
        self.specfmt = "s{}_o{}"

    def __getitem__(self, item):
        with open(self.filename) as f:
            base = yaml.safe_load(f)
            return base[item]

    def __getattribute__(self, item):
        # Short circuit on filename or else infinite recursion
        if item == 'filename':
            return super().__getattribute__(item)
        with open(self.filename) as f:
            base = yaml.safe_load(f)
            if item == 'paths':
                paths = base['paths']
                for key, value in paths.items():
                    if isinstance(value, list):
                        paths[key] = os.path.join(*value)
                return paths
            elif item == 'bands':
                return base['bands']
            elif item == 'custom_bands':
                return base['custom_bands']
            elif item == 'atmmodel':
                return base['atmmodel']
            elif item == 'cache':
                cache = base['cache']
                if (cache["location"] is None) or (cache["location"] == "None"):
                    cache["location"] = None  # Disables caching
                else:
                    if isinstance(cache["location"], list):
                        cache["location"] = os.path.join(*cache["location"])
                return cache
            else:
                return super().__getattribute__(item)