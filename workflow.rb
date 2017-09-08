#!/usr/bin/env ruby

require 'rbbt-util'
require 'rbbt/workflow'

Workflow.require_workflow "Sequence"

module Annovar
  extend Workflow
  Rbbt.claim Rbbt.share.software[".source"]["annovar.tar.gz"], :proc do |filename|
    raise "Please get Annovar from 'http://www.openbioinformatics.org/annovar/annovar_download_form.php' and place it here #{filename}"
  end

  Rbbt.claim Rbbt.share.software.annovar, :proc do |dir|
    dir = dir.find
    tar = Rbbt.share.software[".source"]["annovar.tar.gz"].produce.find
    Misc.in_dir(dir) do
      CMD.cmd("tar xvfz '#{tar.find}'; mv annovar/* .;rmdir annovar")
      CMD.cmd("./annotate_variation.pl --buildver hg19 --downdb ensGene '#{dir.humandb.find}'")
      CMD.cmd("./annotate_variation.pl --buildver hg19 --downdb seq '#{dir.humandb.hg19_seq.find}'")
      CMD.cmd("./retrieve_seq_from_fasta.pl '#{dir.humandb["hg19_ensGene.txt"]}' -seqdir '#{dir.humandb.hg19_seq.find}' -format ensGene -outfile '#{dir.humandb["hg19_ensGeneMrna.fa"]}'")
      CMD.cmd("./annotate_variation.pl --buildver hg18 --downdb ensGene '#{dir.humandb.find}'")
      CMD.cmd("./annotate_variation.pl --buildver hg18 --downdb seq '#{dir.humandb.hg18_seq.find}'")
      CMD.cmd("./retrieve_seq_from_fasta.pl '#{dir.humandb["hg18_ensGene.txt"]}' -seqdir '#{dir.humandb.hg18_seq.find}' -format ensGene -outfile '#{dir.humandb["hg18_ensGeneMrna.fa"]}'")
    end
  end

  SOFTWARE_DIR=Rbbt.share.software.annovar.produce.find

  dep Sequence, :reference
  task :prepare => :array do |mutations|
    TSV.traverse step(:reference).grace, :type => :array, :into => :stream do |line|
      next if line =~ /^#/
        mutation, ref = line.split "\t"
      chr, pos, mut = mutation.split(":")
      [chr, pos, pos, ref, mut, ""]  * "\t"
    end
  end

  dep :prepare
  task :analysis => :tsv do 
    script = SOFTWARE_DIR["annotate_variation.pl"].find
    db = SOFTWARE_DIR["humandb"].find
    log :annovar, "Running annovar script"
    out = file(:out)
    FileUtils.mkdir_p files_dir
    start = Time.now
    `#{script} -out #{out} -build hg19 #{step(:prepare).join.path.find} #{db} -dbtype ensGene`
    set_info :process_time, Time.now - start
    Open.read(file(:out). + '.exonic_variant_function')
  end

  dep :analysis
  task :mutated_isoforms => :tsv do 
    organism = step(:analysis).step(:prepare).step(:reference).inputs["organism"]
    dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => ["Mutated Isoform"], :type => :flat, :organism => organism
    dumper.init
    enst2enp = Organism.transcripts(organism).index :target => "Ensembl Protein ID", :fields => ["Ensembl Transcript ID"], :persist => true
    TSV.traverse step(:analysis), :type => :array, :into => dumper do |line|
      consequence, type, transcripts, chr, start, eend, ref, alt = line.split "\t"
      mutation = [chr,start,alt] * ":"
      next if transcripts == "UNKNOWN"
      mis = transcripts.split(",").collect do |info|
        gene, t, ex, b, change = info.split ":"
        if change
          change.sub!('p.','')
          change.sub!('X','*')
          change.sub!('fs','FrameShift')
        else
          change = b
        end
        protein = enst2enp[t]
        [protein, change] * ":"
      end
      [mutation, mis]
    end
  end

  dep :prepare
  input :database, :string, "Database to annotate with", "phastConsElements46way"
  input :release, :string, "Hg release", "hg19"
  task :annotate => :tsv do |database,release|
    script = SOFTWARE_DIR["annotate_variation.pl"].find
    db = SOFTWARE_DIR["humandb"].find

    `#{script} -build hg19 -downdb #{database} #{db}` unless db.glob(release + '_' + database+".*").any?
    log :annovar, "Running annovar script"
    out = file(:out)
    FileUtils.mkdir_p files_dir
    `#{script}  -regionanno -build hg19 -out #{out} -dbtype #{database} #{step(:prepare).join.path.find} #{db}`
    Open.read(file(:out). + '.hg19_phastConsElements46way')
  end
end
