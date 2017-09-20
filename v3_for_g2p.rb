
#Trim sequence to V3 region(non-alignment), translate to amino acid sequence, and seperate in 50s
#"AHC" can be added at the end of sequence
# count identical sequence numbers.
# aa sequence V3 loop only.

class Sequence
  def initialize (name = ">sequence",dna_sequence ="")
    @name = name
    @dna_sequence = dna_sequence.upcase
    @aa_sequence = ""
    @aa_array = []
  end
  attr_accessor :name, :dna_sequence, :aa_sequence, :aa_array
  def rev_complement
    @dna_sequence.reverse.upcase.tr('ATCG','TAGC')
  end
  def rev_complement!
    @dna_sequence = @dna_sequence.reverse.upcase.tr('ATCG','TAGC')
  end
  def get_aa_sequence(initial_position = 0)
    @aa_sequence = ""
    require_sequence = @dna_sequence[initial_position..-1]
    base_array = []
    require_sequence.each_char {|base| base_array<<base}
    while (base_array.length>=3) do
      base_3= ""
      3.times{base_3 += base_array.shift}
      @aa_sequence<< amino_acid(base_3)
    end
  end
  
  #get amino acid calls, return a array.keep ambiguity calls. 
  def get_aa_array(initial_position = 0)
    @aa_array = []
    require_sequence = @dna_sequence[initial_position..-1].tr('-','N')
    base_array = []
    require_sequence.each_char {|base| base_array<<base}
    while (base_array.length>=3) do
      base_3= ""
      3.times{base_3 += base_array.shift}
      @aa_array<< amino_acid_2(base_3)
    end
  end
  def dna_length
    @dna_sequence.length
  end
  def aa_length
    @aa_sequence.length
  end
end

def amino_acid (bases)
  case bases
  when /^TT[TCY]$/
    return "F"
  when /^TT[AGR]$/
    return "L"
  when /^CT.$/
    return "L"
  when /^AT[TCAHYWM]$/
    return "I"
  when "ATG"
    return "M"
  when /^GT.$/
    return "V"
  when /^TC.$/
    return "S"
  when /^CC.$/
    return "P"
  when /^AC.$/
    return "T"
  when /^GC.$/
    return "A"
  when /^TA[TCY]$/
    return "Y"
  when /^TA[AGR]$/ 
    return "*"
  when /^T[GR]A$/
    return "*"
  when /^CA[TCY]$/
    return "H"
  when /^CA[AGR]$/
    return "Q"
  when /^AA[TCY]$/
    return "N"
  when /^AA[AGR]$/
    return "K"
  when /^GA[TCY]$/
    return "D"
  when /^GA[AGR]$/
    return "E"
  when /^TG[TCY]$/
    return "C"
  when "TGG"
    return "W"
  when /^CG.$/
    return "R"
  when /^AG[TCY]$/
    return "S"
  when /^[AM]G[AGR]$/
    return "R"
  when /^GG.$/
    return "G"
  when /^[ATW][CGS][CTY]$/
    return "S"
  when /^[TCY]T[AGR]$/
    return "L"
  else
    return "#"
  end
end

def fasta_to_hash(infile)
  f=File.open(infile,"r")
  return_hash = {}
  name = ""
  while line = f.gets do
    if line =~ /^\>/
      name = line.chomp
      return_hash[name] = ""
    else
      return_hash[name] += line.chomp
    end
  end
  f.close
  return return_hash
end

terminal = "AHC"
indir = ARGV[0]

outdir2 = indir + "_v3"

Dir.mkdir(outdir2) unless File.directory?(outdir2)

samples = []
Dir.chdir(indir) {samples = Dir.glob("*")}

samples.each do |lib|
  path = indir + "/" + lib

  out_lib = outdir2 + "/" + lib
  Dir.mkdir(out_lib) unless File.directory?(out_lib)
  out_file2 = out_lib + "/" + "v3_aa"
  out_file3 = out_lib + "/" + "v3_aa_list"
  out_file5 = out_lib + "/" + lib
  out2 = File.open(out_file2,"w")
  out3 = File.open(out_file3,"w")
  out5 = File.open(out_file5,"w")
  sequences = fasta_to_hash(path)
  v3_aa_seq = {}
  sequences.each do |k,v|
    tags = k.split("_")
    new_name = tags[0..-2].join("_").gsub("index","RSC0")

    out5.puts new_name + "\n" + v
    v3_nt = v[-105..-1]
    v3 = Sequence.new(k,v3_nt)
    v3.get_aa_sequence
    next if v3.aa_sequence.match(/\*/)
    v3_aa_seq[new_name] = v3.aa_sequence
    out2.puts new_name
    out2.puts v3.aa_sequence
  end
  v3_uni_name = {}
  v3_aa_seq.each do |k,v|
    if v3_uni_name[v]
      v3_uni_name[v] << k
    else
      v3_uni_name[v] = []
      v3_uni_name[v] << k
    end
  end
  v3_uni_name2 = {}
  num = 0
  v3_uni_name.each do |k,v|
    num += 1
    aa_seq_name = ">V3_aa_" + num.to_s + "_" + v.size.to_s
    v3_uni_name2[aa_seq_name] = [k,v]
  end
  num2 = 0
  num3 = 1
  out_file4 = out_lib + "/v3_aa_g2p_" + num3.to_s
  out4 = File.open(out_file4,"w")
  v3_uni_name2.each do |k,v|
    out3.puts k + "\t" + v[0] + "\t" + v[1].join("+")
    num2 += 1
    if num2 < 50
      out4.puts k
      out4.puts v[0] + terminal
    elsif num2 == 50
      out4.puts k
      out4.puts v[0] + terminal
      out4.close
      num2 = 0
      num3 += 1
      out_file4 = out_lib + "/v3_aa_g2p_" + num3.to_s
      out4 = File.open(out_file4,"w")
    end
  end
  
  puts lib + "\t" + sequences.size.to_s + "\t" + v3_aa_seq.size.to_s + "\t" + v3_uni_name.size.to_s
  
end

