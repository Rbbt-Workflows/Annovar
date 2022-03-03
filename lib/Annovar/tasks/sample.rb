module Sample

  dep :vcf_file
  dep :organism
  dep_task :annovar, Annovar, :analysis, :annovar_vcf => :vcf_file, :positions => :placeholder, :vcf => true, :organism => :organism, :full_reference_sequence => false do |jobname,options,dependencies| 
    organism = dependencies.flatten.select{|d| d.task_name.to_s == "organism" }.first.load.strip
    options[:reference_build] = Organism.hg_build(organism)
    {:inputs => options}
  end
end
