import re
def extract(s):
    s = s.replace("'","")
    pattern = r"(\d+)\s+(ENERGY|Energy) = ([-.\d]+)\s+([.\w]+)"
    matchObj = re.match(pattern,s) 
    length = matchObj.group(1)
    energy = matchObj.group(3)
    name = matchObj.group(4)
    return name+"|"+length+"|"+energy
