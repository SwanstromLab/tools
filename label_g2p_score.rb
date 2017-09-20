
#label sequence name with g2p scores
#count lineages above 1% abundance

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

indir = ARGV[0]
outdir = indir + "_g2p"
Dir.mkdir(outdir) unless File.directory?(outdir)
libs = []
Dir.chdir(indir){libs = Dir.glob("*")}
lineage_g2p = {}
libs.each do |lib|
  path = indir + "/" + lib
  g2p_file = path + "/g2p.csv"
  sequence_file = path + "/" + lib
  list_file = path + "/v3_aa_list"
  
  g2p_line = File.readlines(g2p_file)
 
  g2p_line.shift
  g2p_hash = {}
  g2p_line.each do |v|
    line = v.split(",")
    aa_name = ">" + line[0]
    g2p_hash[aa_name] = line[3].to_f
  end

  labeled_name = {}
  list_line = File.readlines(list_file)
  list_line.each do |v|
    line = v.chomp.split("\t")
    aa_name = line[0]
    seq_names = line[2].split("+")
    g2p_hash.each do |k,v1|
      if aa_name == k
        tag = v1.to_s + "_"
        if v1 > 20
          tag += "@"
        elsif v1 <= 2
          tag += "$"
        elsif v1 <= 5
          tag += "%"
        elsif v1 <= 10
          tag += "*"
        elsif v1 <= 20
          tag += "#"
        end
        seq_names.each do |seq_name|
          labeled_name[seq_name] = seq_name + "_" + tag
        end
      end
    end
  end
  new_sequence_file = sequence_file + "_labeled"
  new_sequence_file2 = outdir + "/" + lib
  out = File.open(new_sequence_file,"w")
  out2 = File.open(new_sequence_file2,"w")
  sequences = fasta_to_hash(sequence_file)
  sequences.each do |k,v|
    next unless labeled_name[k]
    out.puts labeled_name[k]
    out2.puts labeled_name[k]
    out.puts v
    out2.puts v
  end
  out.close
  out2.close
  
  lineage_g2p[lib] = []
  total_seq_number = sequences.size
  g2p_hash.each do |k,v|
    seq_number = k.split("_")[-1].to_f
    if seq_number/total_seq_number.to_f > 0.001
      lineage_g2p[lib] << [v,seq_number/total_seq_number.to_f]
    end
  end
end


lineage_g2p.each do |k,v|
  v.each do |v1|
    puts k + "\t" + v1[0].to_s + "\t" + v1[1].round(3).to_s
  end
end

