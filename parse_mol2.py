

def process_meat(meat,repls):
  out_l = []
  i = 0
  for line in meat.split("\n"):
    if line.strip() == "":
      out_l.append(line)
      continue
    string_ind = repls[i]
    buffer = (4-len(string_ind))*" "
    line = "".join((line[:50],string_ind+buffer,line[54:]))
    out_l.append(line)
    i+=1
  return "\n".join(out_l)

def join_mol2(start,middle,end):
  """Join the pieces"""
  top = "@<TRIPOS>ATOM".join([start,middle])
  return "@".join([top,end])

def parse_mol2(input_file,repls):
  """Parse the Mol2 file"""
  input = input_file.read()
  spl = input.split("@<TRIPOS>ATOM")
  b4 = spl[0]
  after = "@".join(spl[1].split("@")[1:])
  meat = spl[1].split("@")[0]
  meat = process_meat(meat,repls)
  return b4,meat,after

def test():
  repls = ['N', 'CAb', 'CAb', 'CAb', 'CAb', 'CAb', 'CT', 'CT', 'CAb', 'Cb', 'O', 'OH', 'H', 'HA', 'HA', 'HA', 'HC', 'HC', 'HC', 'HC', 'HO']
  run("input.mol2",repls,"output.mol2")

def run(input_fp,repls,output_fp):
  data = open(input_fp)
  b,m,e = parse_mol2(data,repls)
  output = join_mol2(b,m,e)
  out_w = open(output_fp,"w")
  out_w.write(output)

if __name__ == "__main__":
  test()
