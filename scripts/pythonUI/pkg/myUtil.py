
def trim_trailing_slash(s:str):
    if s[-1] == '/':
        return s[:-1]
    return s