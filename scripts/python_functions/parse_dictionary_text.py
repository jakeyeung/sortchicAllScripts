import json
import gzip
# import yaml

def process_header(line):
  if line.startswith("@"):
    line = line[1:]  # remove @ at beginning
  return(line)

def inf2dic(inf):
  # inf from computeMatrix
  with gzip.open(inf, "rt") as f:
    line = next(f).rstrip()  # remove \n at end
    line = process_header(line)
    dic = text2dic(line)
  return(dic)

def text2dic(line):
  # dic = json.loads(line)
  # print(line)
  dic = json.loads(line)
  return(dic)
