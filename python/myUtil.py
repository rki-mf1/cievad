import os


def trim_trailing_slash(s:str):
    if s[-1] == '/':
        return s[:-1]
    return s


def mkdir_if_not_present(cfg_path:str):
    dir_exists = os.path.exists(cfg_path)
    if not dir_exists:
        os.makedirs(cfg_path)
        print("New config directory is created at " + cfg_path + "!")
