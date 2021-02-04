#!/usr/bin/env Rscript
#!/Library/Frameworks/R.framework/Resources/bin/Rscript
	#dada2 requires R >= 3.3

###############################################################################
#
#    exec_dada2.R
#
#	 See usage by calling exec_dada2.R -h 
#    
#    Copyright (C) 2018 Mikihiko Kawai
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

args0<-commandArgs(trailingOnly=T)
	# both args[1] and args[[1]] works
	##################
	# -flag test -set this file -q 2 -p=1
	# argv: test, this, file, 2 (args['q'] is 1 (flag) not 2)
	##################
argv<-c()
keys <-c()
vals <-c()
for (arg in args0){
	print(arg)
	if (length(grep('^-',arg))){
		if (length(grep('^-.+=.+',arg))){
			kv<-strsplit(arg,'=')[[1]]
			key<-kv[1]
			val<-kv[2]
		}else{	# -flag
			key<-arg
			val<-1
		}
		key<-sub('^-','',key)
		keys<-c(keys, key)
		vals <-c(vals,val)
	}else{
		argv<-c(argv,arg)
	}
}
#print(vals)
args_my <- as.list(vals)
#print(args_my)
names(args_my)<-keys

version <- 20190816
#version <- 20190815
#version <- 20190207
version_history <- paste(
"20190816
  - re-set sample.names, fnFs, fnRs, filtFs, and filtRs when some samples have no reads after filtering. 
  - not exit even when some samples have no reads after filtering. 
  - not fail even when part of 'Five seqtab files to check' are executed again, by setting from rds files of each step. 
  - bug fix: save dadaFs correctly (before, both dadaFs and dadaRs files have dadaRs data. each sample rds dadaF and dadaR files were correct). 
  - save errF and errR rds files.
  - output formatted summary.derepFastq.txt.
20190815
  - check point for f_filter is separated from f_construct_sequence_table. 
  - when some samples have no reads after filtering, exit before f_plot_filter to avoid making truncated Rplot.qual_profile.filt.pdf. 
  - fixed: not to fail in f_remove_chimeras, when zero-read samples are removed in the second run after failure of the first run after f_filt (by using join instead of cbind for making summary.seqtab.txt, even when qual_profile.filt.rds is made by the first run with records of zero-read samples). 
20190207
  - phyloseq barplot without facet_wrap
",
		sep = "\n"
)
if (
	!is.null(args_my['v'][[1]])
		||
	!is.null(args_my['version'][[1]])
){
	cat('version ', version, "\n", "\n", "version history: \n", "ver. ", version_history, "\n", sep = "")
	quit()
}

msg <- paste(
"key_list:
	dir, 
		ext, 
	tax, 
	trimLefts, truncLens
	
options:

example:
exec_dada2.R -path -tax=
",
		sep = "\n"
)
if (
	!is.null(args_my['h'][[1]])
		||
	is.null(args_my['dir'][[1]])
){
	cat(msg, "\n")
	quit()
}

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dirname <- dirname(script.name)
#print(script.dirname)

data_d <- paste(script.dirname, 'data', sep = '/')
#data_d <- "/Volumes/HDD_QuadRAIDerEZ/recog_data/etc/src/_RECOG-server-MK/R/RNAseq/bam2gene/data"

obj <- list(
	ext = '.fastq',	#'.fastq.gz$'
	
	trimLefts = '0,0', 
	#truncLens = '290,215', 
	ees = '2,2', 
	#ee = 1,
	n = 0,
	truncQ = 2,
	
	tax = '/home/local/db/dada2/silva_nr_v123_train_set.fa.gz',
		#$HOME/local/db/dada2/silva_nr_v132_train_set.fa.gz
	dummy = 1
		# Remember to eliminate comma for the last entry. 
)

	## key_list for obj
key_list_str_chr <- c()
key_list_str_num <- c()
#key_list_str_num <- c('ylim')
key_list <- c(
	'dir', 
		'ext',	#extension for fastq files (".fastq.gz$")
	'trimLefts', 
	'truncLens', 
	'ees', 
	#'ee', 
	'n', 'truncQ', 
	'tax', 
	#'in',
	#'out',
	'outdir', 'outfiltdir', 
	'dummy'
		# Remember to eliminate comma for the last entry.
) 

for (key in key_list){
	if (exists("args_my") && !is.null(args_my[key][[1]])){
		obj[[key]] <- as.character(args_my[key])
	} 
}
for (key in key_list_str_chr){
	key_str <- paste(key, 'str', sep = '_')
	#obj[[key_str]] <- obj[[key]]
	obj[[key]] <- strsplit(obj[[key_str]], ',')[[1]]
}
for (key in key_list_str_num){
	key_str <- paste(key, 'str', sep = '_')
	#obj[[key_str]] <- obj[[key]]
	obj[[key]] <- as.numeric(strsplit(obj[[key_str]], ',')[[1]])
}
f_template <- function(obj){
		#obj <- f_derep(obj)
	#outpath <- obj[['outdir']]
	outfiltdir <- obj[['outfiltdir']]
	sample.names <- obj[['sample.names']]
	fnFs <- obj[['fnFs']]
	#fnRs <- obj[['fnRs']]
	filtFs <- obj[['filtFs']]
	filtRs <- obj[['filtRs']]
	
	...
	
	obj[['derepFs']] <- derepFs
	obj[['derepRs']] <- derepRs
	
	return(obj)
}

f_init_in_fq_files <- function(obj){
		#obj <- f_init_in_fq_files(obj)
		#in case some samples have no reads after filtering (f_filter), sample.names, fnFs, fnRs, filtFs, and filtRs are re-set in f_filter. 
	path <- obj[['dir']]
	ext_fq <- obj[['ext']]
	
	
	fns <- list.files(path)
	cat("File names under the given dir, ", path, "\n")
	print(fns)
	cat("\n")

	
	fastqs <- fns[grepl(ext_fq, fns)]
	#fastqs <- fns[grepl(".fastq.gz$", fns)]
	fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
	fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
	fnRs <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files
	# Get sample names from the first part of the forward read filenames
	#sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
	
	cat("File names for forward reads (assuming given ext_fq and '_R1')\n")
	print(fnFs)
	cat("\n")
	
	cat("File names for reverse reads (assuming given ext_fq and '_R2')\n")
	print(fnRs)
	cat("\n")
	

	# Fully specify the path for the fnFs and fnRs
	fnFs_path <- file.path(path, fnFs)
	fnRs_path <- file.path(path, fnRs)
	
	#obj[['sample.names']] <- sample.names
	obj[['baseFs']] <- fnFs
	obj[['baseRs']] <- fnRs
	obj[['fnFs']] <- fnFs_path
	obj[['fnRs']] <- fnRs_path
	
	return(obj)
}

f_init_sample_names <- function(obj){
		#obj <- f_init_sample_names(obj)
	fnFs <- obj[['baseFs']]
	#fnFs <- obj[['fnFs']]
	
	sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
	
	cat("Sample names taken from file names for forward reads)\n")
	print(sample.names)
	cat("Note. filterAndTrim does not use above sample names, so samples names for filter count are independently set in f_init_from_filt_rds by the same way as done here.\n")
	#cat("Note. filterAndTrim does not use above sample names, so samples names for filter count are independently set in f_summarize by the same way as done here.\n")
	cat("Note. in case some samples have no reads after filtering (f_filter), sample.names, fnFs, fnRs, filtFs, and filtRs will be re-set later via f_init_from_filt_rds, (1) from two if clauses of main to skip f_filter, or (2) in f_filter after filtering of f_filter. \n")
		#now should be set via f_init_from_filt_rds, (1) from two if clauses of main to skip f_filter, or (2) in f_filter after filtering of f_filter. 
	cat("\n")
	
	obj[['sample.names']] <- sample.names
	
	return(obj)
}

f_init_qual_profile <- function(obj){
		#obj <- f_init_qual_profile(obj)
	outpath <- obj[['outdir']]
	
	f_qual_profile_raw <- file.path(outpath, "Rplot.qual_profile.raw.pdf")
	
	obj[['file_qual_profile_raw']] <- f_qual_profile_raw
	
	return(obj)
}
f_qual_profile <- function(obj){
		#obj <- f_qual_profile(obj)
	outpath <- obj[['outdir']]
	sample.names <- obj[['sample.names']]
	fnFs <- obj[['fnFs']]
	fnRs <- obj[['fnRs']]
	f_qual_profile_raw <- obj[['file_qual_profile_raw']]
	
	pdf(f_qual_profile_raw, onefile = T)
	#pdf(file.path(outpath, "Rplot.qual_profile.raw.pdf"), onefile = T)
	#pdf(file.path(path, "Rplot.qual_profile.raw.pdf"), onefile = T)
	#pdf("Rplot.qual_profile.pdf", onefile = T)
	for(i in seq_along(fnFs)) {
	#for (i in c(1:4)){
		sample.name <- sample.names[[i]]
		cat(paste0(sample.name,'.. '))
		#cat(paste0(sample.name,'..',"\n"))
		
		gF <- plotQualityProfile(fnFs[[i]], n = 10000)
		gF <- gF + labs(title = paste0(sample.name, ' ', '(Fwd)'))
		gR <- plotQualityProfile(fnRs[[i]], n = 10000)
		gR <- gR + labs(title = paste0(sample.name, ' ', '(Rev)'))
		print(gF)
		print(gR)
		#print(plotQualityProfile(fnFs[[i]], n = 10000))
		#print(plotQualityProfile(fnRs[[i]], n = 10000))
		#print(plotQualityProfile(fnFs[[i]]))
		#print(plotQualityProfile(fnRs[[i]]))
	}
	dev.off()
	# plotQualityProfile(fnFs[[1]])
	# plotQualityProfile(fnFs[[2]])
	
	# plotQualityProfile(fnRs[[1]])
	# plotQualityProfile(fnRs[[2]])
	
	return(obj)
}
f_init_filter_parameters <- function(obj){
		#obj <- f_init_filter_parameters(obj)
	#outpath <- obj[['outdir']]
	
	keys_trim <- c('trimLefts', 'truncLens', 'ees')
	for (key in keys_trim){
		if (
			is.null(obj[[key]])
		){
			cat(msg, "\n")
			cat(paste0("!!!\n", "key ", key, " is not set. \n", "!!!\n"))
			quit()
		}
	}
	
	trimLeft <- as.numeric(strsplit(obj[['trimLefts']], ',')[[1]])
	truncLen <- as.numeric(strsplit(obj[['truncLens']], ',')[[1]])
	ee <- as.numeric(strsplit(obj[['ees']], ',')[[1]])
	
	trimLeft
	truncLen
	ee
	
	obj[['trimLeft']] <- trimLeft
	obj[['truncLen']] <- truncLen
	obj[['ee']] <- ee
	
	return(obj)
}
f_init_outfiltdir <- function(obj){
		#obj <- f_init_filtpath(obj)
	#ee <- obj[['ee']]
	n <- obj[['n']]
	truncQ <- obj[['truncQ']]
	
	obj <- f_init_filter_parameters(obj)	#comma-separated string 'trimLefts' to array of numerics 'trimLeft'. 
	trimLeft <- obj[['trimLeft']]	#array of numerics 'trimLeft'
	truncLen <- obj[['truncLen']]
	ee <- obj[['ee']]	#array of numerics 'ee'
	
	str_trimLeft <- paste(trimLeft, collapse = '_')
	str_truncLen <- paste(truncLen, collapse = '_')
	str_ee <- paste(ee, collapse = '_')
	#str_trimLeft <- paste(trimLeft, collapse = ',')
	#str_truncLen <- paste(truncLen, collapse = ',')
	base <- paste(
		"filtered", 
		".tleft", str_trimLeft, '.tlen', str_truncLen, 
		".ee", str_ee, 
		#".ee", ee, 
		".tq", truncQ, ".n", n, 
		sep = ""
	)
	#base <- paste("filtered", ".tleft", str_trimLeft, '.tlen', str_truncLen, ".ee", ee, ".tq2.n", n, sep = "")
		#tq for truncQ, fixed as 2 now. 
	#base <- paste("filtered", ".tleft16,20.tlen290,215.ee", ee, ".tq2.n", n, sep = "")
	outfiltdir <- file.path(outpath, base)
	#outfiltdir <- file.path(path, base)
	if(!file_test("-d", outfiltdir)) dir.create(outfiltdir)
	
	obj[['outfiltdir']] <- outfiltdir
	
	return(obj)
}

f_init_filter <- function(obj){
		#obj <- f_init_filter(obj)
	#outpath <- obj[['outdir']]
	outfiltdir <- obj[['outfiltdir']]
	sample.names <- obj[['sample.names']]
	# fnFs <- obj[['fnFs']]
	# fnRs <- obj[['fnRs']]
	# trimLeft <- obj[['trimLeft']]
	# truncLen <- obj[['truncLen']]
	# ee <- obj[['ee']]
	# n <- obj[['n']]
	# truncQ <- obj[['truncQ']]
	
	#ee <- 1
	#ee <- 2
	#n <- 0
	#truncQ <- 2
	
	# str_trimLeft <- paste(trimLeft, collapse = '_')
	# str_truncLen <- paste(truncLen, collapse = '_')
	# #str_trimLeft <- paste(trimLeft, collapse = ',')
	# #str_truncLen <- paste(truncLen, collapse = ',')
	# base <- paste(
		# "filtered", 
		# ".tleft", str_trimLeft, '.tlen', str_truncLen, 
		# ".ee", ee, ".tq", truncQ, ".n", n, 
		# sep = ""
	# )
	# #base <- paste("filtered", ".tleft", str_trimLeft, '.tlen', str_truncLen, ".ee", ee, ".tq2.n", n, sep = "")
		# #tq for truncQ, fixed as 2 now. 
	# #base <- paste("filtered", ".tleft16,20.tlen290,215.ee", ee, ".tq2.n", n, sep = "")
	# outfiltdir <- file.path(outpath, base)
	# #outfiltdir <- file.path(path, base)
	# if(!file_test("-d", outfiltdir)) dir.create(outfiltdir)
	
	path_filt_fq <- file.path(outfiltdir, 'filt.fq')
	if(!file_test("-d", path_filt_fq)) dir.create(path_filt_fq)
	#cat("filtFs..\n")
	#print(sample.names)
	#paste0(sample.names, "_F_filt.fastq.gz")
	filtFs <- file.path(path_filt_fq, paste0(sample.names, "_F_filt.fastq.gz"))
	filtRs <- file.path(path_filt_fq, paste0(sample.names, "_R_filt.fastq.gz"))
	names(filtFs) <- sample.names
	names(filtRs) <- sample.names
	
	#cat("f_qual_profile_filt..\n")
	#f_qual_profile_filt <- file.path(outfiltdir, "Rplot.qual_profile.filt.pdf")
	
	# obj[['ee']] <- ee
	# obj[['n']] <- n
	#obj[['outfiltdir']] <- outfiltdir
	obj[['filtFs']] <- filtFs
	obj[['filtRs']] <- filtRs
	#obj[['file_qual_profile_filt']] <- f_qual_profile_filt
	
	return(obj)
}


f_init_construct_sequence_table <- function(obj){
		#obj <- f_init_construct_sequence_table(obj)
	#outpath <- obj[['outdir']]
	outfiltdir <- obj[['outfiltdir']]
	
	f_qual_profile_filt <- file.path(outfiltdir, "Rplot.qual_profile.filt.pdf")
	file_qual_profile_filt_rds <- file.path(outfiltdir, "qual_profile.filt.rds")
	file_derepFs_rds <- file.path(outfiltdir, "derepFs.rds")
	file_derepRs_rds <- file.path(outfiltdir, "derepRs.rds")
	file_dadaFs_rds <- file.path(outfiltdir, "dadaFs.rds")
	file_dadaRs_rds <- file.path(outfiltdir, "dadaRs.rds")
	file_errF_rds <- file.path(outfiltdir, "errF.rds")
	file_errR_rds <- file.path(outfiltdir, "errR.rds")
	file_mergers_rds <- file.path(outfiltdir, "mergers.rds")
	
	f_seqtab <- file.path(outfiltdir, "seqtab.txt")
	f_seqtab_fa <- file.path(outfiltdir, "seqtab.fa")
	f_seqtab_asis <- file.path(outfiltdir, "seqtab.asis.txt")
	f_seqtab_otu <- file.path(outfiltdir, "seqtab.otu_table.txt")
	f_seqtab_rds <- file.path(outfiltdir, "seqtab.rds")
	
	
	obj[['file_qual_profile_filt']] <- f_qual_profile_filt
	obj[['file_qual_profile_filt_rds']] <- file_qual_profile_filt_rds
	obj[['file_derepFs_rds']] <- file_derepFs_rds
	obj[['file_derepRs_rds']] <- file_derepRs_rds
	
	obj[['file_dadaFs_rds']] <- file_dadaFs_rds
	obj[['file_dadaRs_rds']] <- file_dadaRs_rds
	obj[['file_errF_rds']] <- file_errF_rds
	obj[['file_errR_rds']] <- file_errR_rds
	
	obj[['file_mergers_rds']] <- file_mergers_rds
	
	obj[['file_seqtab']] <- f_seqtab
	obj[['file_seqtab_fa']] <- f_seqtab_fa
	obj[['file_seqtab_asis']] <- f_seqtab_asis
	obj[['file_seqtab_otu']] <- f_seqtab_otu
	obj[['file_seqtab_rds']] <- f_seqtab_rds
	
	return(obj)
}

f_init_from_filt_rds <- function(obj){
		#obj <- f_init_from_filt_rds(obj)
	path <- obj[['dir']]
	file_qual_profile_filt_rds <- obj[['file_qual_profile_filt_rds']]
	sample.names <- obj[['sample.names']]
	baseFs0 <- obj[['baseFs']]
	baseRs0 <- obj[['baseRs']]
	#out <- obj[['filterAndTrim']]
	
	out <- readRDS(file_qual_profile_filt_rds)

	out2 <- as.data.frame(out)
	samples_out <- row.names(out2)
	# Get sample names from the first part of the forward read filenames
	samples_out <- sapply(strsplit(samples_out, "_"), `[`, 1)
	#print(samples_out)
	out3 <- data.frame(samples = samples_out, out2)
	colnames(out3) <- c('sample', 'input', 'filtered')
	
	out_zero <- out3[out3[,'filtered'] == 0,]
	#out_zero <- out[out[,'reads.out'] == 0,]
	if (0){
	if (length(out_zero[,1]) >= 1){
		cat("some samples are zero-reads after filter. check results, remove corresponding fastq files from dir of raw reads, then execute again.\n")
		print(out_zero)
		cat("be cautious that if you add new samples under dir of raw reads for second run, rather than removing samples from there, you have to manually remove 'qual_profile.filt.rds' file.\n")
		quit()
	}
	}
	if (length(out_zero[,1]) >= 1){
		msg1 <- "Some samples do not have reads after filter. \n  they are omitted from sample.names parameter here, to avoid failure at later steps. \n  see 'samples_with_zero_read_after_filter.txt' files for samples with no reads and omitted here.\n"
		cat(msg1)
		#cat("some samples do not have reads after filter. \n  they are omitted from sample.names parameter here, to avoid failure at later steps. \n  see 'samples_with_zero_read_after_filter.txt' files for samples with no reads and omitted here.\n")
		#cat("some samples are zero-reads after filter. check results, remove corresponding fastq files from dir of raw reads, then execute again.\n")
		print(out_zero)
		msg2 <- "Be cautious that if you add new samples under dir of raw reads for second run, rather than removing samples from there, you have to manually remove 'qual_profile.filt.rds' file.\n"
		cat(msg2, "\n", sep = "")
		
		# f_warning_message_zero_aft_filt <- file.path(outfiltdir, "warninng_msg_samples_with_zero_read_after_filter.txt")
		# es_msg <- data.frame(msg = c(msg1, msg))
		# # f_seqtab_nochim_otu <- file.path(outfiltdir, "seqtab_nochim.fa.otu_table.txt")
		# write.table(es_msg, file = f_warning_message_zero_aft_filt, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
		
		# f_list_zero_aft_filt <- file.path(outfiltdir, "samples_with_zero_read_after_filter.txt")
		# #out_zero1 <- as.data.frame(sample = row.names(out_zero), out_zero)
		# write.table(out_zero, file = f_list_zero_aft_filt, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
		
		
		sample.names.wo_zero <- as.vector(out3[out3[,'filtered'] != 0, 'sample'])
		#sample.names.wo_zero <- row.names(out[out[,'reads.out'] != 0, ])
		
		cat("Re-set sample.names.\n")
		cat("Update sample names, taken from samples which have reads after filtering\n")
		print(sample.names.wo_zero)
		#obj[['sample.names']] <- sample.names.wo_zero
		cat("Note. filterAndTrim does not use above sample names, so samples names for filter count are independently set in f_summarize by the same way as done here.\n")
		cat("\n")
		
		cat("Update fnFs and fnRs.\n")
		# baseFs0 <- obj[['baseFs']]
		# baseRs0 <- obj[['baseRs']]
		es_fn0 <- data.frame(sample = sample.names, baseF = baseFs0, baseR = baseRs0)
			#original sample.names are taken from baseFs, and the order of them has not been changed. 
		es_fn <- es_fn0[es_fn0[, 'sample'] %in% sample.names.wo_zero, ]
		fnFs <- as.vector(es_fn[,'baseF'])
		fnRs <- as.vector(es_fn[,'baseR'])
		#fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
		#fnRs <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files
		# Get sample names from the first part of the forward read filenames
		#sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
		
		cat("File names for forward reads (assuming given ext_fq and '_R1', with new set of sample names)\n", "i.e. without zero-read samples, ", row.names(out_zero), "\n", sep = "")
		print(fnFs)
		cat("\n")
		
		cat("File names for reverse reads (assuming given ext_fq and '_R2', with new set of sample names)\n", "i.e. without zero-read samples, ", row.names(out_zero), "\n", sep = "")
		print(fnRs)
		cat("\n")
		
		# Fully specify the path for the fnFs and fnRs
		fnFs_path <- file.path(path, fnFs)
		fnRs_path <- file.path(path, fnRs)
		
		obj[['sample.names']] <- sample.names.wo_zero
		#obj[['sample.names']] <- sample.names
		obj[['baseFs']] <- fnFs
		obj[['baseRs']] <- fnRs
		obj[['fnFs']] <- fnFs_path
		obj[['fnRs']] <- fnRs_path
		
		# print("fnFs_path:\n")
		# print(fnFs_path)
		# print(obj[['fnFs']])
	
		cat("Update filtFs and filtRs by calling f_init_filter again.\n", "\n", sep = "")
		obj <- f_init_filter(obj)
			# obj[['filtFs']] <- filtFs
			# obj[['filtRs']] <- filtRs
	
		#cat(msg2)
		#quit()
		# print("fnFs_path from obj:\n")
		# print(obj[['fnFs']])
		
	}else{
		cat(
			"All samples have reads after filter.\n", 
			"Use values set before for parameters listed below.\n", 
			"  sample.names (set in f_init_sample_names)\n",
			"  baseFs, baseRs, fnFs, fnRs (set in f_init_in_fq_files)\n", 
			"  filtFs, filtRs (set in f_init_filter)\n", 
			"Parameters set here from rds are: filterAndTrim_my, filterAndTrim_my_zero\n", 
			"\n", 
			sep = ""
		)
		cat("File names for forward reads (as set before)\n", "\n", sep = "")
		print(baseFs0)
		cat("\n")
		
		cat("File names for reverse reads (as set before)\n", "\n", sep = "")
		print(baseRs0)
		cat("\n")
		
		#print(obj[['fnFs']])
	}
	
	obj[['filterAndTrim_my']] <- out3
	obj[['filterAndTrim_my_zero']] <- out_zero
	#obj[['filterAndTrim']] <- out
	
	return(obj)
}
	
f_filter <- function(obj){
		#obj <- f_filter(obj)
	#outpath <- obj[['outdir']]
	outfiltdir <- obj[['outfiltdir']]
	sample.names <- obj[['sample.names']]
	fnFs <- obj[['fnFs']]
	fnRs <- obj[['fnRs']]
	filtFs <- obj[['filtFs']]
	filtRs <- obj[['filtRs']]
	trimLeft <- obj[['trimLeft']]
	truncLen <- obj[['truncLen']]
	ee <- obj[['ee']]
	n <- obj[['n']]
	truncQ <- obj[['truncQ']]
	file_qual_profile_filt <- obj[['file_qual_profile_filt']]
	file_qual_profile_filt_rds <- obj[['file_qual_profile_filt_rds']]
	
	cat("Filter..\n")
	# Filter
	out <- filterAndTrim(
		fnFs, filtFs, 
		fnRs, filtRs, 
		truncLen = truncLen,
		maxN = n, maxEE = ee, truncQ = truncQ, 
		rm.phix = TRUE,
		compress = TRUE, multithread = FALSE
		#compress = TRUE, multithread = TRUE
	)
		#filterAndTrim is a multithreaded convenience interface for the fastqFilter and fastqPairedFilter filtering functions. 
		#Note that error messages and tracking are not handled gracefully when using the multithreading functionality. 
		#If errors arise, it is recommended to re-run without multithreading to troubleshoot the issue.
	#out
	print(head(out))
	#head(out)
	saveRDS(out, file = file_qual_profile_filt_rds)
		# On Windows set multithread=FALSE

	obj <- f_init_from_filt_rds(obj)
		#read file_qual_profile_filt_rds, set filterAndTrim_my, filterAndTrim_my_zero
	
	out_zero <- obj[['filterAndTrim_my_zero']]
	if (length(out_zero[,1]) >= 1){
		# msg1 <- "some samples do not have reads after filter. \n  they are omitted from sample.names parameter here, to avoid failure at later steps. \n  see 'samples_with_zero_read_after_filter.txt' files for samples with no reads and omitted here.\n"
		# cat(msg1)
		# #cat("some samples do not have reads after filter. \n  they are omitted from sample.names parameter here, to avoid failure at later steps. \n  see 'samples_with_zero_read_after_filter.txt' files for samples with no reads and omitted here.\n")
		# #cat("some samples are zero-reads after filter. check results, remove corresponding fastq files from dir of raw reads, then execute again.\n")
		# print(out_zero)
		# msg2 <- "be cautious that if you add new samples under dir of raw reads for second run, rather than removing samples from there, you have to manually remove 'qual_profile.filt.rds' file.\n"
		# cat(msg2)
		
		f_warning_message_zero_aft_filt <- file.path(outfiltdir, "warning_msg_samples_with_zero_read_after_filter.txt")
		es_msg <- data.frame(msg = c(msg1, msg))
		# f_seqtab_nochim_otu <- file.path(outfiltdir, "seqtab_nochim.fa.otu_table.txt")
		write.table(es_msg, file = f_warning_message_zero_aft_filt, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
		
		f_list_zero_aft_filt <- file.path(outfiltdir, "samples_with_zero_read_after_filter.txt")
		#out_zero1 <- as.data.frame(sample = row.names(out_zero), out_zero)
		write.table(out_zero, file = f_list_zero_aft_filt, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
	}

	
	# for(i in seq_along(fnFs)) {
		# sample.name <- sample.names[[i]]
		# cat(paste0(sample.name,'.. '))
		
		# fastqPairedFilter(
			# c(fnFs[i], fnRs[i]), 
			# c(filtFs[i], filtRs[i]),
			# trimLeft = trimLeft,
			# truncLen = truncLen, 
			# #trimLeft = c(16, 20), 
			# #truncLen = c(290,215), 
			# maxN = n, 
			# maxEE = ee, 
			# truncQ = truncQ, 
			# #truncQ = 2, 
			# compress = TRUE, verbose = TRUE
		# )
	# }
	# cat("\n")
	
		# # Learn the Error Rates
	# errF <- learnErrors(filtFs, multithread = FALSE)
	# errR <- learnErrors(filtRs, multithread = FALSE)
	# #errF <- learnErrors(filtFs, nbases = 5e8, multithread = FALSE, randomize=TRUE)
	# #errR <- learnErrors(filtRs, nbases = 5e8, multithread = FALSE, randomize=TRUE)
	# #errF <- learnErrors(filtFs, multithread = TRUE)
	# #errR <- learnErrors(filtRs, multithread = TRUE)
	# #err <- learnErrors(filts, nbases = 1e8, multithread=TRUE, randomize=TRUE)
	
	# cat("pdf plotErrors..\n")
	# pdf(file.path(outfiltdir, "Rplot.estimated_error_rates.pdf"), onefile = T)
	# #pdf(file_plotErrors, onefile = T)
	# g1 <- plotErrors(errF, nominalQ = TRUE)
	# g2 <- plotErrors(errR, nominalQ = TRUE)
	# print(g1)
	# print(g2)
	# dev.off()
	# cat("\n")
	
	# cat("pdf file_qual_profile_filt..\n")
	# pdf(file_qual_profile_filt, onefile = T)
	# #pdf(file.path(outfiltdir, "Rplot.qual_profile.filt.pdf"), onefile = T)
	# #pdf("Rplot.qual_profile.pdf", onefile = T)
	# for(i in seq_along(fnFs)) {
	# #for (i in c(1:4)){
		# sample.name <- sample.names[[i]]
		# cat(paste0(sample.name,'.. '))
		
		# gF <- plotQualityProfile(filtFs[[i]], n = 10000)
		# gF <- gF + labs(title = paste0(sample.name, ' ', '(Fwd)'))
		# gR <- plotQualityProfile(filtFs[[i]], n = 10000)
		# gR <- gR + labs(title = paste0(sample.name, ' ', '(Rev)'))
		# print(gF)
		# print(gR)
		# #print(plotQualityProfile(filtFs[[i]], n = 10000))
		# #print(plotQualityProfile(filtRs[[i]], n = 10000))
		# #print(plotQualityProfile(filtFs[[i]]))
		# #print(plotQualityProfile(filtRs[[i]]))
	# }
	# dev.off()
	# cat("\n")
	
	#obj[['filterAndTrim']] <- out
	# obj[['errF']] <- errF
	# obj[['errR']] <- errR
	# obj[['outfiltdir']] <- outfiltdir
	# obj[['filtFs']] <- filtFs
	# obj[['filtRs']] <- filtRs
	
	return(obj)
}
f_plot_filter <- function(obj){
		#obj <- f_filter(obj)
	#outpath <- obj[['outdir']]
	outfiltdir <- obj[['outfiltdir']]
	sample.names <- obj[['sample.names']]
	fnFs <- obj[['fnFs']]
	fnRs <- obj[['fnRs']]
	filtFs <- obj[['filtFs']]
	filtRs <- obj[['filtRs']]
	trimLeft <- obj[['trimLeft']]
	truncLen <- obj[['truncLen']]
	ee <- obj[['ee']]
	n <- obj[['n']]
	truncQ <- obj[['truncQ']]
	file_qual_profile_filt <- obj[['file_qual_profile_filt']]
	file_qual_profile_filt_rds <- obj[['file_qual_profile_filt_rds']]
	
	# cat("Filter..\n")
	# # Filter
	# out <- filterAndTrim(
		# fnFs, filtFs, 
		# fnRs, filtRs, 
		# truncLen = truncLen,
		# maxN = n, maxEE = ee, truncQ = truncQ, 
		# rm.phix = TRUE,
		# compress = TRUE, multithread = FALSE
		# #compress = TRUE, multithread = TRUE
	# )
		# #filterAndTrim is a multithreaded convenience interface for the fastqFilter and fastqPairedFilter filtering functions. 
		# #Note that error messages and tracking are not handled gracefully when using the multithreading functionality. 
		# #If errors arise, it is recommended to re-run without multithreading to troubleshoot the issue.
	# #out
	# head(out)
	# saveRDS(out, file = file_qual_profile_filt_rds)
		# # On Windows set multithread=FALSE
	
	# for(i in seq_along(fnFs)) {
		# sample.name <- sample.names[[i]]
		# cat(paste0(sample.name,'.. '))
		
		# fastqPairedFilter(
			# c(fnFs[i], fnRs[i]), 
			# c(filtFs[i], filtRs[i]),
			# trimLeft = trimLeft,
			# truncLen = truncLen, 
			# #trimLeft = c(16, 20), 
			# #truncLen = c(290,215), 
			# maxN = n, 
			# maxEE = ee, 
			# truncQ = truncQ, 
			# #truncQ = 2, 
			# compress = TRUE, verbose = TRUE
		# )
	# }
	# cat("\n")
	
		# # Learn the Error Rates
	# errF <- learnErrors(filtFs, multithread = FALSE)
	# errR <- learnErrors(filtRs, multithread = FALSE)
	# #errF <- learnErrors(filtFs, nbases = 5e8, multithread = FALSE, randomize=TRUE)
	# #errR <- learnErrors(filtRs, nbases = 5e8, multithread = FALSE, randomize=TRUE)
	# #errF <- learnErrors(filtFs, multithread = TRUE)
	# #errR <- learnErrors(filtRs, multithread = TRUE)
	# #err <- learnErrors(filts, nbases = 1e8, multithread=TRUE, randomize=TRUE)
	
	# cat("pdf plotErrors..\n")
	# pdf(file.path(outfiltdir, "Rplot.estimated_error_rates.pdf"), onefile = T)
	# #pdf(file_plotErrors, onefile = T)
	# g1 <- plotErrors(errF, nominalQ = TRUE)
	# g2 <- plotErrors(errR, nominalQ = TRUE)
	# print(g1)
	# print(g2)
	# dev.off()
	# cat("\n")
	
	cat("pdf file_qual_profile_filt..\n")
	pdf(file_qual_profile_filt, onefile = T)
	#pdf(file.path(outfiltdir, "Rplot.qual_profile.filt.pdf"), onefile = T)
	#pdf("Rplot.qual_profile.pdf", onefile = T)
	for(i in seq_along(fnFs)) {
	#for (i in c(1:4)){
		sample.name <- sample.names[[i]]
		cat(paste0(sample.name,'.. '))
		
		gF <- plotQualityProfile(filtFs[[i]], n = 10000)
		gF <- gF + labs(title = paste0(sample.name, ' ', '(Fwd)'))
		gR <- plotQualityProfile(filtFs[[i]], n = 10000)
		gR <- gR + labs(title = paste0(sample.name, ' ', '(Rev)'))
		print(gF)
		print(gR)
		#print(plotQualityProfile(filtFs[[i]], n = 10000))
		#print(plotQualityProfile(filtRs[[i]], n = 10000))
		#print(plotQualityProfile(filtFs[[i]]))
		#print(plotQualityProfile(filtRs[[i]]))
	}
	dev.off()
	cat("\n")
	
	# obj[['filterAndTrim']] <- out
	# obj[['errF']] <- errF
	# obj[['errR']] <- errR
	# obj[['outfiltdir']] <- outfiltdir
	# obj[['filtFs']] <- filtFs
	# obj[['filtRs']] <- filtRs
	
	return(obj)
}

f_derep <- function(obj){
		#obj <- f_derep(obj)
	#outpath <- obj[['outdir']]
	outfiltdir <- obj[['outfiltdir']]
	sample.names <- obj[['sample.names']]
	fnFs <- obj[['fnFs']]
	#fnRs <- obj[['fnRs']]
	filtFs <- obj[['filtFs']]
	filtRs <- obj[['filtRs']]
	file_derepFs_rds <- obj[['file_derepFs_rds']]
	file_derepRs_rds <- obj[['file_derepRs_rds']]
	
	#print(filtFs)
	
	derepFs <- vector("list", length(sample.names))
	derepRs <- vector("list", length(sample.names))
	# Name the derep-class objects by the sample names
	names(derepFs) <- sample.names
	names(derepRs) <- sample.names
	for (name in sample.names) {
		cat("Processing:", name, "\n")
		derepFs[[name]] <- derepFastq(filtFs[[name]], verbose = TRUE)
		derepRs[[name]] <- derepFastq(filtRs[[name]], verbose = TRUE)
		#print(filtFs[[name]])
		#derepFs[[name]] <- derepFastq(filtFs[[name]])
		#derepRs[[name]] <- derepFastq(filtRs[[name]])
	}
	#derepFs <- derepFastq(filtFs, verbose = TRUE)
	#derepRs <- derepFastq(filtRs, verbose = TRUE)
	
	f_summary_derep <- file.path(outfiltdir, "summary.derepFastq.txt")
	#f_summary_derep <- "summary.derepFastq.txt"
	for(i in seq_along(fnFs)) {
		strF <- summary(derepFs[[i]])
		strR <- summary(derepRs[[i]])
		write.table(strF, file = f_summary_derep, append = TRUE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
		write.table(strR, file = f_summary_derep, append = TRUE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
		#write(strF, file = f_summary_derep, append = T)
		#write(strR, file = f_summary_derep, append = T)
			#The data (usually a matrix) x are written to file file. If x is a two-dimensional matrix you need to transpose it to get the columns in file the same as those in the internal representation.
	}
	#derepFs[[1]]
	# derep-class: R object describing dereplicated sequencing reads
	# $uniques: 7124 reads in 1860 unique sequences
	  # Sequence lengths: min=230, median=230, max=230
	# $quals: Quality matrix dimension:  1860 230
	  # Consensus quality scores: min=12, median=37.75, max=39
	# $map: Map from reads to unique sequences:  2 49 1 62 1 ...
	
	path_rds_derep <- file.path(outfiltdir, 'rds.derep')
	if(!file_test("-d", path_rds_derep)) dir.create(path_rds_derep)
	for(i in seq_along(fnFs)) {
		f_obj_F <- file.path(path_rds_derep, paste0(sample.names[i], ".derepF.rds"))
		saveRDS(derepFs[[i]], file = f_obj_F)
		f_obj_R <- file.path(path_rds_derep, paste0(sample.names[i], ".derepR.rds"))
		saveRDS(derepRs[[i]], file = f_obj_R)
	}
	
	#f_derepFs_rds <- file.path(outfiltdir, derepFs.rds)
	#f_derepRs_rds <- file.path(outfiltdir, derepRs.rds)
	saveRDS(derepFs, file = file_derepFs_rds)
	saveRDS(derepRs, file = file_derepRs_rds)
	
	obj[['derepFs']] <- derepFs
	obj[['derepRs']] <- derepRs
	
	return(obj)
}

f_infer_sample <- function(obj){
		#obj <- f_infer_sample(obj)
	# Sample Inference
	# very slow
	
	#outpath <- obj[['outdir']]
	outfiltdir <- obj[['outfiltdir']]
	sample.names <- obj[['sample.names']]
	file_dadaFs_rds <- obj[['file_dadaFs_rds']]
	file_dadaRs_rds <- obj[['file_dadaRs_rds']]
	file_errF_rds <- obj[['file_errF_rds']]
	file_errR_rds <- obj[['file_errR_rds']]
	fnFs <- obj[['fnFs']]
		#use seq_along(fnFs) during saving each dada rds files. version <= 20190815 has global fnFs set in main, so did not fail. 
	#fnRs <- obj[['fnRs']]
	filtFs <- obj[['filtFs']]
	filtRs <- obj[['filtRs']]
	# errF <- obj[['errF']]
	# errR <- obj[['errR']]
	derepFs <- obj[['derepFs']]
	derepRs <- obj[['derepRs']]
	
		# Learn the Error Rates
	cat("learnErrors..\n")
	errF <- learnErrors(filtFs, multithread = FALSE)
	errR <- learnErrors(filtRs, multithread = FALSE)
	#errF <- learnErrors(filtFs, nbases = 5e8, multithread = FALSE, randomize=TRUE)
	#errR <- learnErrors(filtRs, nbases = 5e8, multithread = FALSE, randomize=TRUE)
	#errF <- learnErrors(filtFs, multithread = TRUE)
	#errR <- learnErrors(filtRs, multithread = TRUE)
	#err <- learnErrors(filts, nbases = 1e8, multithread=TRUE, randomize=TRUE)
	cat("\n")
	
	cat("pdf plotErrors..\n")
	pdf(file.path(outfiltdir, "Rplot.estimated_error_rates.pdf"), onefile = T)
	#pdf(file_plotErrors, onefile = T)
	g1 <- plotErrors(errF, nominalQ = TRUE)
	g2 <- plotErrors(errR, nominalQ = TRUE)
	print(g1)
	print(g2)
	dev.off()
	cat("\n")
	
		# Infer sequence variants
	# derepFs <- vector("list", length(sample.names))
	# derepRs <- vector("list", length(sample.names))
	dadaFs <- vector("list", length(sample.names))
	dadaRs <- vector("list", length(sample.names))
	# names(derepFs) <- sample.names
	# names(derepRs) <- sample.names
	names(dadaFs) <- sample.names
	names(dadaRs) <- sample.names
	for (name in sample.names) {
		cat("Processing:", name, "\n")

		dadaFs[[name]] <- dada(derepFs[[name]], err = errF, multithread = FALSE)
		dadaRs[[name]] <- dada(derepRs[[name]], err = errR, multithread = FALSE)
		#dadaFs[[name]] <- dada(derepFs[[name]], err = errF, multithread = TRUE)
	}
	cat("done processing.\n", "\n", sep = "")
	#dadaFs <- dada(derepFs, err = errF, multithread = FALSE)
	#dadaRs <- dada(derepRs, err = errR, multithread = FALSE)
	#dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
	#dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
	# dadaFs <- dada(derepFs, err = NULL, selfConsist = FALSE)
	# dadaRs <- dada(derepRs, err = NULL, selfConsist = FALSE)
	# dadaFs <- dada(derepFs, err = NULL, selfConsist = TRUE)
	# dadaRs <- dada(derepRs, err = NULL, selfConsist = TRUE)
	
	print(dadaFs[[1]])
	#dadaFs[[1]]
	cat("\n")
	
	cat("save rds files..\n")
	path_rds_dada2 <- file.path(outfiltdir, 'rds.dada2')
	if(!file_test("-d", path_rds_dada2)) dir.create(path_rds_dada2)
	for(i in seq_along(fnFs)) {
		f_obj_F <- file.path(path_rds_dada2, paste0(sample.names[i], ".dadaF.rda"))
		saveRDS(dadaFs[[i]], file = f_obj_F)
		f_obj_R <- file.path(path_rds_dada2, paste0(sample.names[i], ".dadaR.rda"))
		saveRDS(dadaRs[[i]], file = f_obj_R)
	}
	
	# file_dadaFs_rds <- file.path(outfiltdir, "dadaFs.rds")
	# file_dadaRs_rds <- file.path(outfiltdir, "dadaRs.rds")
	saveRDS(dadaFs, file = file_dadaFs_rds)
	saveRDS(dadaRs, file = file_dadaRs_rds)
	
	saveRDS(errF, file = file_errF_rds)
	saveRDS(errR, file = file_errR_rds)
	
	#save(mod, file = "mymodel.rda")
	#saveRDS
	#foo <- readRDS("/path/to/run1/output/seqtab.rds").

	# Visualize estimated error rates:
	#plotErrors(dadaFs[[1]], nominalQ=TRUE)
	#plotErrors(dadaFs[[2]], nominalQ=TRUE)
	
	# pdf(file.path(outfiltdir, "Rplot.estimated_error_rates.pdf"), onefile = T)
	# #for(i in seq_along(fnFs)) {
	# plotErrors(dadaFs[[1]], nominalQ=TRUE)
	# plotErrors(dadaFs[[2]], nominalQ=TRUE)
	# #}
	# dev.off()
	
	# obj[['errF']] <- errF
	# obj[['errR']] <- errR
	obj[['dadaFs']] <- dadaFs
	obj[['dadaRs']] <- dadaRs
	# obj[['file_dadaFs_rds']] <- file_dadaFs_rds
	# obj[['file_dadaRs_rds']] <- file_dadaRs_rds
	
	return(obj)
}

f_merge_paired_reads <- function(obj){
		#obj <- f_merge_paired_reads(obj)
	#outpath <- obj[['outdir']]
	outfiltdir <- obj[['outfiltdir']]
	sample.names <- obj[['sample.names']]
	file_mergers_rds <- obj[['file_mergers_rds']]
	fnFs <- obj[['fnFs']]
	#fnRs <- obj[['fnRs']]
	dadaFs <- obj[['dadaFs']]
	dadaRs <- obj[['dadaRs']]
	derepFs <- obj[['derepFs']]
	derepRs <- obj[['derepRs']]
	
	mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
	# Inspect the merger data.frame from the first sample
	head(mergers[[1]])
	
	path_mergers <- file.path(outfiltdir, 'mergers')
	if(!file_test("-d", path_mergers)) dir.create(path_mergers)
	
	#sample.names 
	#samples <- names(mergers)
	for(i in seq_along(fnFs)) {
		sample <- sample.names[i]
		names <- rownames(mergers[[i]])
		mergers1 <- cbind(name = c(names), mergers[[i]][c(2:length(mergers[[i]][1,]),1)])
		f_mergers <- file.path(path_mergers, paste("mergers", ".", sample, ".txt", sep = ""))
		write.table(mergers1, file = f_mergers, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
	}
	
	#file_mergers_rds <- file.path(outfiltdir, "mergers.rds")
	saveRDS(mergers, file = file_mergers_rds)
	#file_mergers_rds <- obj[['file_mergers_rds']]
	
	#obj[['file_mergers_rds']] <- file_mergers_rds
	obj[['mergers']] <- mergers
	
	return(obj)
}


f_construct_sequence_table <- function(obj){
		#obj <- f_construct_sequence_table(obj)
	#outpath <- obj[['outdir']]
	#outfiltdir <- obj[['outfiltdir']]
	#sample.names <- obj[['sample.names']]
	#fnFs <- obj[['fnFs']]
	#fnRs <- obj[['fnRs']]
	
	file_mergers_rds <- obj[['file_mergers_rds']]
	
	f_seqtab <- obj[['file_seqtab']]
	f_seqtab_fa <- obj[['file_seqtab_fa']]
	f_seqtab_asis <- obj[['file_seqtab_asis']]
	f_seqtab_otu <- obj[['file_seqtab_otu']]
	f_seqtab_rds <- obj[['file_seqtab_rds']]
	
	mergers <- readRDS(file_mergers_rds)
	
	seqtab <- makeSequenceTable(mergers)
	#seqtab <- makeSequenceTable(mergers[names(mergers) != "Mock"])
	dim(seqtab)
	# [1]    4 2521
	# Inspect distribution of sequence lengths
	table(nchar(getSequences(seqtab)))
	#table(nchar(colnames(seqtab)))
	# 290 294 298 299 312 318 335 351 363 372 375 387 406 407 409 411 412 413 414 
	  # 1   1  27   2   1   1   1   1   1   1   1   1   1  65  20 502  34 232 224 
	# 415 438 440 441 442 443 444 450 451 452 454 455 459 463 465 468 
	 # 42 148 483 237 109   5   2 180  21  81   1  86   1   1   4   3
	
	t.seqtab <- t(seqtab)
	t.seqtab1 <- as.data.frame(cbind(t.seqtab, seq = rownames(t.seqtab)))
	#f_seqtab <- file.path(outfiltdir, "seqtab.txt")
	write.table(t.seqtab1, file = f_seqtab, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
	
	#labels <- paste('dada2', seq(0,length(t.seqtab[,1])-1), sep = "_")
	labels <- paste0("dada2_", seq(1, length(getUniques(seqtab))), ";size=", unname(getUniques(seqtab)), ";")
	#labels <- paste0("dada2_", seq(0, length(getUniques(seqtab))-1), ";size=", unname(getUniques(seqtab)), ";")
	#f_seqtab_fa <- file.path(outfiltdir, "seqtab.fa")
	uniquesToFasta(seqtab, f_seqtab_fa, ids = labels, mode = "w")
	
	#f_seqtab_asis <- file.path(outfiltdir, "seqtab.asis.txt")
	write.table(seqtab, file = f_seqtab_asis, append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
	#seqtab.recv <- read.table(f_seqtab_asis, header = T, sep = "\t", quote = "", fill = TRUE, comment.char = "")
	
	#t.seqtab.nochim <- t(seqtab.nochim)
	t.seqtab.otu <- as.data.frame(cbind('#OTU ID' = labels, t.seqtab))
	#f_seqtab_otu <- file.path(outfiltdir, "seqtab.fa.otu_table.txt")
	write.table(t.seqtab.otu, file = f_seqtab_otu, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
	
	saveRDS(seqtab, f_seqtab_rds) # CHANGE ME to where you want sequence table saved
	#saveRDS(seqtab, "path/to/study/seqtab_final.rds") # CHANGE ME to where you want sequence table saved
	
	obj[['seqtab']] <- seqtab
	
	return(obj)
}

f_init_path_remove_chimeras <- function(obj){
		#obj <- f_init_path_remove_chimeras(obj)
	#outpath <- obj[['outdir']]
	outfiltdir <- obj[['outfiltdir']]
	#sample.names <- obj[['sample.names']]
	#fnFs <- obj[['fnFs']]
	#fnRs <- obj[['fnRs']]
	#seqtab <- obj[['seqtab']]
	
	f_seqtab_nochim <- file.path(outfiltdir, "seqtab_nochim.txt")
	f_seqtab_nochim_fa <- file.path(outfiltdir, "seqtab_nochim.fa")
	f_seqtab_nochim_asis <- file.path(outfiltdir, "seqtab_nochim.asis.txt")
	f_seqtab_nochim_otu <- file.path(outfiltdir, "seqtab_nochim.otu_table.txt")
	f_seqtab_nochim_rds <- file.path(outfiltdir, "seqtab_nochim.rds")
	f_summary_seqtab <- file.path(outfiltdir, "summary.seqtab.txt")
	f_summary_seqlen <- file.path(outfiltdir, "summary.seqlength.txt")
	
	obj[['file_seqtab_nochim']] <- f_seqtab_nochim
	obj[['file_seqtab_nochim_fa']] <- f_seqtab_nochim_fa
	obj[['file_seqtab_nochim_asis']] <- f_seqtab_nochim_asis
	obj[['file_seqtab_nochim_otu']] <- f_seqtab_nochim_otu
	obj[['file_seqtab_nochim_rds']] <- f_seqtab_nochim_rds
	obj[['file_summary_seqtab']] <- f_summary_seqtab
	obj[['file_summary_seqlen']] <- f_summary_seqlen
	
	return(obj)
}
f_init_remove_chimeras <- function(obj){
		#obj <- f_init_remove_chimeras(obj)
	#outpath <- obj[['outdir']]
	#outfiltdir <- obj[['outfiltdir']]
	sample.names <- obj[['sample.names']]
	#fnFs <- obj[['fnFs']]
	#fnRs <- obj[['fnRs']]
	#seqtab <- obj[['seqtab']]
	# file_qual_profile_filt_rds <- obj[['file_qual_profile_filt_rds']]
	# file_dadaFs_rds <- obj[['file_dadaFs_rds']]
	# file_dadaRs_rds <- obj[['file_dadaRs_rds']]
	# file_mergers_rds <- obj[['file_mergers_rds']]
	#f_seqtab_asis <- obj[['file_seqtab_asis']]
	file_seqtab_rds <- obj[['file_seqtab_rds']]
	
	if (is.null(obj[['seqtab']])){
		cat("seqtab object is not set (this case occurs when f_construct_sequence_table step is skipped), try to set it from file..\n")
		#seqtab.recv <- read.table(f_seqtab_asis, header = T, sep = "\t", quote = "", fill = TRUE, comment.char = "")
		print(file_seqtab_rds)
		seqtab <- readRDS(file_seqtab_rds)
		
		obj[['seqtab']] <- seqtab
	}
	
	return(obj)
}
f_remove_chimeras <- function(obj){
		#obj <- f_remove_chimeras(obj)
	#outpath <- obj[['outdir']]
	outfiltdir <- obj[['outfiltdir']]
	#sample.names <- obj[['sample.names']]
	#fnFs <- obj[['fnFs']]
	#fnRs <- obj[['fnRs']]
	seqtab <- obj[['seqtab']]
	# out <- obj[['filterAndTrim']]
	# derepFs <- obj[['derepFs']]
	# derepRs <- obj[['derepRs']]
	# dadaFs <- obj[['dadaFs']]
	# dadaRs <- obj[['dadaRs']]
	# mergers <- obj[['mergers']]
	seqtab.nochim <- obj[['seqtab.nochim']]
	
	f_seqtab_nochim <- obj[['file_seqtab_nochim']]
	f_seqtab_nochim_fa <- obj[['file_seqtab_nochim_fa']]
	f_seqtab_nochim_asis <- obj[['file_seqtab_nochim_asis']]
	f_seqtab_nochim_otu <- obj[['file_seqtab_nochim_otu']]
	f_seqtab_nochim_rds <- obj[['file_seqtab_nochim_rds']]
	# f_summary_seqtab <- obj[['file_summary_seqtab']]
	# f_summary_seqlen <- obj[['file_summary_seqlen']]
	# f_seqtab_nochim <- file.path(outfiltdir, "seqtab_nochim.txt")
	# f_seqtab_nochim_fa <- file.path(outfiltdir, "seqtab_nochim.fa")
	# f_seqtab_nochim_asis <- file.path(outfiltdir, "seqtab_nochim.asis.txt")
	# f_seqtab_nochim_otu <- file.path(outfiltdir, "seqtab_nochim.fa.otu_table.txt")
	# f_summary_seqtab <- file.path(outfiltdir, "summary.seqtab.txt")
	# f_summary_seqlen <- file.path(outfiltdir, "summary.seqlength.txt")
	
	seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE, verbose = TRUE)
	#seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
	#seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)
	#seqtab.nochim
	
	saveRDS(seqtab.nochim, f_seqtab_nochim_rds) # CHANGE ME to where you want sequence table saved
	#saveRDS(seqtab, "path/to/study/seqtab_final.rds") # CHANGE ME to where you want sequence table saved
	
	#summary(seqtab.nochim)
	cat("dim(seqtab.nochim): \n")
	print(dim(seqtab.nochim))
	# [1]   4 162
	
	cat("sum(seqtab.nochim) / sum(seqtab): \n")
	print(sum(seqtab.nochim) / sum(seqtab))
	# [1] 0.1822704
	
	cat("table(nchar(colnames(seqtab.nochim))): \n")
	print(table(nchar(colnames(seqtab.nochim))))
	# 290 294 298 299 312 318 335 351 363 372 375 387 406 407 409 411 412 413 414 
	  # 1   1   1   2   1   1   1   1   1   1   1   1   1   3   1  26  10  24  10 
	# 415 438 440 441 442 443 444 450 451 452 454 455 459 463 465 468 
	  # 2   7  19  11   8   2   2   7   1   2   1   4   1   1   3   3 
	
	t.seqtab.nochim <- t(seqtab.nochim)
	t.seqtab.nochim1 <- as.data.frame(cbind(t.seqtab.nochim, seq = rownames(t.seqtab.nochim)))
	# f_seqtab_nochim <- file.path(outfiltdir, "seqtab_nochim.txt")
	write.table(t.seqtab.nochim1, file = f_seqtab_nochim, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
		#t.seqtab.nochim1.recv <- read.table(f_seqtab_nochim, header = T, sep = "\t", quote = "", fill = TRUE, comment.char = "")
	
	#f_seqtab_nochim_fa <- file.path(outfiltdir, "seqtab_nochim.uchime.fa")
	#uniquesToFasta(seqtab.nochim, f_seqtab_nochim_fa, ids = NULL, mode = "w")
	
	#labels <- paste('dada2', seq(0,length(t.seqtab.nochim[,1])-1), sep = "_")
	labels <- paste0("dada2_", seq(0, length(getUniques(seqtab.nochim))-1), ";size=", unname(getUniques(seqtab.nochim)), ";")
	# f_seqtab_nochim_fa <- file.path(outfiltdir, "seqtab_nochim.fa")
	uniquesToFasta(seqtab.nochim, f_seqtab_nochim_fa, ids = labels, mode = "w")
	
	# f_seqtab_nochim_asis <- file.path(outfiltdir, "seqtab_nochim.asis.txt")
	write.table(seqtab.nochim, file = f_seqtab_nochim_asis, append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
	#seqtab.nochim.recv <- read.table(f_seqtab_nochim_asis, header = T, sep = "\t", quote = "", fill = TRUE, comment.char = "")
	
	#t.seqtab.nochim <- t(seqtab.nochim)
	t.seqtab.nochim.otu <- as.data.frame(cbind('#OTU ID' = labels, t.seqtab.nochim))
	# f_seqtab_nochim_otu <- file.path(outfiltdir, "seqtab_nochim.fa.otu_table.txt")
	write.table(t.seqtab.nochim.otu, file = f_seqtab_nochim_otu, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
	
	obj[['seqtab.nochim']] <- seqtab.nochim
	obj[['t.seqtab.nochim1']] <- t.seqtab.nochim1
	
	return(obj)
}
f_init_summarize <- function(obj){
		#obj <- f_init_summarize(obj)
	#outpath <- obj[['outdir']]
	#outfiltdir <- obj[['outfiltdir']]
	sample.names <- obj[['sample.names']]
	#fnFs <- obj[['fnFs']]
	#fnRs <- obj[['fnRs']]
	#seqtab <- obj[['seqtab']]
	file_qual_profile_filt_rds <- obj[['file_qual_profile_filt_rds']]
	file_dadaFs_rds <- obj[['file_dadaFs_rds']]
	file_dadaRs_rds <- obj[['file_dadaRs_rds']]
	file_mergers_rds <- obj[['file_mergers_rds']]
	#f_seqtab_asis <- obj[['file_seqtab_asis']]
	file_seqtab_rds <- obj[['file_seqtab_rds']]
	file_seqtab_nochim_rds <- obj[['file_seqtab_nochim_rds']]
	
		#now should be set via f_init_from_filt_rds, (1) from two if clauses of main to skip f_filter, or (2) in f_filter after filtering of f_filter. 
	# if (is.null(obj[['filterAndTrim']])){
		# cat("filterAndTrim object is not set (this case occurs when f_qual_profile_filt_rds step is skipped), try to set it from file..\n")
		# print(file_qual_profile_filt_rds)
		# out <- readRDS(file_qual_profile_filt_rds)
		# #out0 <- readRDS(file_qual_profile_filt_rds)
		# #out <- out0[row.names(out0) %in% sample.names,]
		
		# obj[['filterAndTrim']] <- out
	# }
	
	if (is.null(obj[['dadaFs']])){
		cat("dadaFs object is not set (this case occurs when f_qual_profile_filt_rds step is skipped), try to set it from file..\n")
		print(file_dadaFs_rds)
		obj[['dadaFs']] <- readRDS(file_dadaFs_rds)
	}
	if (is.null(obj[['dadaRs']])){
		cat("dadaRs object is not set (this case occurs when f_qual_profile_filt_rds step is skipped), try to set it from file..\n")
		print(file_dadaRs_rds)
		obj[['dadaRs']] <- readRDS(file_dadaRs_rds)
	}
	
	if (is.null(obj[['mergers']])){
		cat("mergers object is not set (this case occurs when f_qual_profile_filt_rds step is skipped), try to set it from file..\n")
		print(file_mergers_rds)
		obj[['mergers']] <- readRDS(file_mergers_rds)
	}
	
	if (is.null(obj[['seqtab']])){
		cat("seqtab object is not set (this case occurs when f_construct_sequence_table step is skipped), try to set it from file..\n")
		#seqtab.recv <- read.table(f_seqtab_asis, header = T, sep = "\t", quote = "", fill = TRUE, comment.char = "")
		print(file_seqtab_rds)
		seqtab <- readRDS(file_seqtab_rds)
		
		obj[['seqtab']] <- seqtab
	}
	
	if (is.null(obj[['seqtab.nochim']])){
		cat("seqtab.nochim object is not set (this case occurs when f_remove_chimeras step is skipped), try to set it from file..\n")
		#seqtab.recv <- read.table(f_seqtab_asis, header = T, sep = "\t", quote = "", fill = TRUE, comment.char = "")
		print(file_seqtab_nochim_rds)
		seqtab.nochim <- readRDS(file_seqtab_nochim_rds)
		
		obj[['seqtab.nochim']] <- seqtab.nochim
	}
	
	return(obj)
}
f_summarize <- function(obj){
		#obj <- f_remove_chimeras(obj)
	#outpath <- obj[['outdir']]
	outfiltdir <- obj[['outfiltdir']]
	sample.names <- obj[['sample.names']]
	#fnFs <- obj[['fnFs']]
	#fnRs <- obj[['fnRs']]
	seqtab <- obj[['seqtab']]
	out3 <- obj[['filterAndTrim_my']]
	#out <- obj[['filterAndTrim']]
		#now should be set via f_init_from_filt_rds, (1) from two if clauses of main to skip f_filter, or (2) in f_filter after filtering of f_filter. 
	# derepFs <- obj[['derepFs']]
	# derepRs <- obj[['derepRs']]
	dadaFs <- obj[['dadaFs']]
	dadaRs <- obj[['dadaRs']]
	mergers <- obj[['mergers']]
	seqtab.nochim <- obj[['seqtab.nochim']]
	
	# f_seqtab_nochim <- obj[['file_seqtab_nochim']]
	# f_seqtab_nochim_fa <- obj[['file_seqtab_nochim_fa']]
	# f_seqtab_nochim_asis <- obj[['file_seqtab_nochim_asis']]
	# f_seqtab_nochim_otu <- obj[['file_seqtab_nochim_otu']]
	# f_seqtab_nochim_rds <- obj[['file_seqtab_nochim_rds']]
	f_summary_seqtab <- obj[['file_summary_seqtab']]
	f_summary_seqlen <- obj[['file_summary_seqlen']]
	# f_seqtab_nochim <- file.path(outfiltdir, "seqtab_nochim.txt")
	# f_seqtab_nochim_fa <- file.path(outfiltdir, "seqtab_nochim.fa")
	# f_seqtab_nochim_asis <- file.path(outfiltdir, "seqtab_nochim.asis.txt")
	# f_seqtab_nochim_otu <- file.path(outfiltdir, "seqtab_nochim.fa.otu_table.txt")
	# f_summary_seqtab <- file.path(outfiltdir, "summary.seqtab.txt")
	# f_summary_seqlen <- file.path(outfiltdir, "summary.seqlength.txt")
	
	# # f_summary_seqtab <- file.path(outfiltdir, "summary.seqtab.txt")
	if (0){
	es <- c()
	e1 <- c(name = 'seqtab', nseq = sum(seqtab), nrep = dim(seqtab)[2])
	e2 <- c(name = 'seqtab.nochim', nseq = sum(seqtab.nochim), nrep = dim(seqtab.nochim)[2])
	es <- rbind(es, e1,e2)
	write.table(es, file = f_summary_seqtab, append = F, quote = F, sep = "\t", row.names = F, col.names = T)
	}
	if (0){
	getN <- function(x) sum(getUniques(x))
	
	#cat("track\n")
	samples <- data.frame(sample = sample.names)
	#print(samples)
	#print(out)
	track <- cbind(samples, out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
	#print(track)
	# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
	colnames(track) <- c("sample", "input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
	#track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
	# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
	#colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
	#rownames(track) <- sample.names
	head(track)
	write.table(track, file = f_summary_seqtab, append = F, quote = F, sep = "\t", row.names = F, col.names = T)
	}
	if (1){	#require package plyr. 
	getN <- function(x) sum(getUniques(x))
	
	#cat("track\n")
	samples <- data.frame(sample = sample.names)
	#print(samples)
	
	#print(out)
	# out2 <- as.data.frame(out)
	# samples_out <- row.names(out2)
	# # Get sample names from the first part of the forward read filenames
	# samples_out <- sapply(strsplit(samples_out, "_"), `[`, 1)
	# #print(samples_out)
	# out3 <- data.frame(samples = samples_out, out2)
	# colnames(out3) <- c('sample', 'input', 'filtered')
	# #print(out3)
	
	df1 <- as.data.frame(sapply(dadaFs, getN))
	df2 <- data.frame(sample = row.names(df1), df1)
	colnames(df2) <- c('sample', 'denoisedF')
	# > head(df2)
	    # sample denoisedF
	# 10M    10M    149983
	# 11M    11M       261
	# 12M    12M     48486
	# 13M    13M     81718
	# 14M    14M     58155
	# 15M    15M     75688
	#print(df2)
	
	dr1 <- as.data.frame(sapply(dadaRs, getN))
	dr2 <- data.frame(sample = row.names(dr1), dr1)
	colnames(dr2) <- c('sample', 'denoisedR')
	#print(dr2)
	
	m1 <- as.data.frame(sapply(mergers, getN))
	m2 <- data.frame(sample = row.names(m1), m1)
	colnames(m2) <- c('sample', 'merged')
	#print(m2)

	nc1 <- rowSums(seqtab.nochim)
	nc2 <- data.frame(sample = names(nc1), nonchim = nc1)
	#print(nc2)

	t1 <- plyr::join(samples, out3, by = 'sample', type = "right");
	t2 <- plyr::join(t1, df2, by = 'sample');
	t3 <- plyr::join(t2, dr2, by = 'sample');
	t4 <- plyr::join(t3, m2, by = 'sample');
	track <- plyr::join(t4, nc2, by = 'sample');
	
	#track <- cbind(samples, out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
	#print(track)
	# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
	#colnames(track) <- c("sample", "input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
	#track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
	# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
	#colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
	#rownames(track) <- sample.names
	head(track)
	write.table(track, file = f_summary_seqtab, append = F, quote = F, sep = "\t", row.names = F, col.names = T)
	}
	
	# f_summary_seqlen <- file.path(outfiltdir, "summary.seqlength.txt")
	con21 <- table(nchar(colnames(seqtab)))
	es21 <- as.data.frame(con21)
	
	con22 <- table(nchar(colnames(seqtab.nochim)))
	es22 <- as.data.frame(con22)
	es2 <- plyr::join(es21, es22, by = c('Var1'))
	colnames(es2) <- c('length', 'seqtab', 'seqtab.nochim')
	write.table(es2, file = f_summary_seqlen, append = F, quote = F, sep = "\t", row.names = F, col.names = T)
	
	
	# obj[['seqtab.nochim']] <- seqtab.nochim
	# obj[['t.seqtab.nochim1']] <- t.seqtab.nochim1
	
	return(obj)
}


f_init_path_assign_taxonomy <- function(obj){
		#obj <- f_init_path_assign_taxonomy(obj)
	#outpath <- obj[['outdir']]
	outfiltdir <- obj[['outfiltdir']]
	#sample.names <- obj[['sample.names']]
	#fnFs <- obj[['fnFs']]
	#fnRs <- obj[['fnRs']]
	
	keys_tax <- c('tax')	#required in f_assign_taxonomy
		#e.g.
		#~/local/db/dada2/silva_nr_v132_train_set.fa.gz
	for (key in keys_tax){
		if (
			is.null(obj[[key]])
		){
			cat(msg, "\n")
			cat(paste0("!!!\n", "key ", key, " is not set. \n", "!!!\n"))
			quit()
		}
	}
	
	base_tax <- sub('\\.gz$', '', basename(obj[['tax']]))
		#silva_nr_v132_train_set.fa.gz -> silva_nr_v132_train_set.fa
	
	f_taxa <- file.path(outfiltdir, paste("seqtab_nochim", base_tax, "txt", sep = "."))
	f_taxa_w_count <- file.path(outfiltdir, paste("seqtab_nochim", base_tax, "w_count.txt", sep = "."))
	f_taxa_asis <- file.path(outfiltdir, paste("seqtab_nochim", base_tax, "asis.txt", sep = "."))
	f_taxa_rds <- file.path(outfiltdir, paste("seqtab_nochim", base_tax, "rds", sep = "."))
	# f_taxa2 <- file.path(outfiltdir, "seqtab_nochim.silva_nr_v123_train_set.fa.txt")
	# f_taxa_w_count <- file.path(outfiltdir, "seqtab_nochim.silva_nr_v123_train_set.fa.w_count.txt")
	# f_taxa2_asis <- file.path(outfiltdir, "seqtab_nochim.silva_nr_v123_train_set.fa.asis.txt")
	
	obj[['file_taxa']] <- f_taxa
	obj[['file_taxa_w_count']] <- f_taxa_w_count
	obj[['file_taxa_asis']] <- f_taxa_asis
	obj[['file_taxa_rds']] <- f_taxa_rds
	
	return(obj)
}

f_init_assign_taxonomy <- function(obj){
		#obj <- f_init_assign_taxonomy(obj)
	#outpath <- obj[['outdir']]
	#outfiltdir <- obj[['outfiltdir']]
	#sample.names <- obj[['sample.names']]
	#fnFs <- obj[['fnFs']]
	#fnRs <- obj[['fnRs']]
	#t.seqtab.nochim1 <- obj[['t.seqtab.nochim1']]
	f_seqtab_nochim_asis <- obj[['file_seqtab_nochim_asis']]
	f_seqtab_nochim <- obj[['file_seqtab_nochim']]
	f_seqtab_nochim_rds <- obj[['file_seqtab_nochim_rds']]
	
	if (is.null(obj[['seqtab.nochim']])){
		cat("seqtab.nochim object is not set (this case occurs when f_remove_chimeras step is skipped), try to set it from file..\n")
		#seqtab.nochim.recv <- read.table(f_seqtab_nochim_asis, header = T, sep = "\t", quote = "", fill = TRUE, comment.char = "")
		seqtab.nochim <- readRDS(f_seqtab_nochim_rds)
		
		obj[['seqtab.nochim']] <- seqtab.nochim
		#obj[['seqtab.nochim']] <- names(seqtab.nochim.recv)
	}
	
	if (is.null(obj[['t.seqtab.nochim1']])){
		cat("t.seqtab.nochim1 object is not set (this case occurs when f_remove_chimeras step is skipped), try to set it from file..\n")
		t.seqtab.nochim1.recv <- read.table(f_seqtab_nochim, header = T, sep = "\t", quote = "", fill = TRUE, comment.char = "")
		
		obj[['t.seqtab.nochim1']] <- t.seqtab.nochim1.recv
	}
	cat("\n")
	
	return(obj)
}

f_assign_taxonomy <- function(obj){
		#obj <- f_assign_taxonomy(obj)
	#outpath <- obj[['outdir']]
	outfiltdir <- obj[['outfiltdir']]
	sample.names <- obj[['sample.names']]
	#fnFs <- obj[['fnFs']]
	#fnRs <- obj[['fnRs']]
	seqtab.nochim <- obj[['seqtab.nochim']]
	t.seqtab.nochim1 <- obj[['t.seqtab.nochim1']]
	tax_file <- obj[['tax']]
	#tax_file <- obj[['tax_file']]
	
	f_taxa <- obj[['file_taxa']]
	f_taxa_w_count <- obj[['file_taxa_w_count']]
	f_taxa_asis <- obj[['file_taxa_asis']]
	f_taxa_rds <- obj[['file_taxa_rds']]
	
	taxa2 <- assignTaxonomy(seqtab.nochim, tax_file, minBoot = 90)
	#http://benjjneb.github.io/dada2/assign.html
	colnames(taxa2) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
	unname(head(taxa2))
	
	taxa22 <- as.data.frame(cbind(taxa2, seq = rownames(taxa2)))
	# f_taxa2 <- file.path(outfiltdir, "seqtab_nochim.silva_nr_v123_train_set.fa.txt")
	write.table(taxa22, file = f_taxa, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
	
	
	labels <- paste0("dada2_", seq(1, length(getUniques(seqtab.nochim))))
	#labels <- paste0("dada2_", seq(0, length(getUniques(seqtab.nochim))-1), ";size=", unname(getUniques(seqtab.nochim)), ";")
	t.seqtab.nochim2 <- as.data.frame(cbind('#OTU ID' = labels, t.seqtab.nochim1))
	taxa_w_count <- plyr::join(t.seqtab.nochim2, taxa22, by = c('seq'))
	taxa_w_count2 <- taxa_w_count[,c('#OTU ID', sample.names, colnames(taxa2), 'seq')]
	#taxa_w_count2 <- taxa_w_count[,c(1:length(sample.names),c(4:length(taxa_w_count[1,])),3)]
	# f_taxa_w_count <- file.path(outfiltdir, "seqtab_nochim.silva_nr_v123_train_set.fa.w_count.txt")
	write.table(taxa_w_count2, file = f_taxa_w_count, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
	
	
	# f_taxa2_asis <- file.path(outfiltdir, "seqtab_nochim.silva_nr_v123_train_set.fa.asis.txt")
	write.table(taxa2, file = f_taxa_asis, append = FALSE, quote = FALSE, sep = "\t", row.names = T, col.names = TRUE)
	
	saveRDS(taxa2, f_taxa_rds)
	#saveRDS(tax, "pathto/study/tax_final.rds") # CHANGE ME ...

	#obj[['derepFs']] <- derepFs
	
	return(obj)
}

f_phyloseq_plot <- function(obj){
		#obj <- f_phyloseq_plot(obj)
	#outpath <- obj[['outdir']]
	outfiltdir <- obj[['outfiltdir']]
	#sample.names <- obj[['sample.names']]
	#fnFs <- obj[['fnFs']]
	#fnRs <- obj[['fnRs']]
	
	# f_seqtab_nochim <- obj[['file_seqtab_nochim']]
	# f_seqtab_nochim_fa <- obj[['file_seqtab_nochim_fa']]
	f_seqtab_nochim_asis <- obj[['file_seqtab_nochim_asis']]
	# f_seqtab_nochim_otu <- obj[['file_seqtab_nochim_otu']]
	# f_summary_seqtab <- obj[['file_summary_seqtab']]
	# f_summary_seqlen <- obj[['file_summary_seqlen']]
	
	# f_taxa <- obj[['file_taxa']]
	# f_taxa_w_count <- obj[['file_taxa_w_count']]
	f_taxa_asis <- obj[['file_taxa_asis']]
	
		
	# Evaluate accuracy
	# DADA2 accuracy on mock community:
	# unqs.mock <- getUniques(removeBimeraDenovo(mergers[["Mock"]], verbose=TRUE))
	
	# cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
	
	# mockRef <- readFasta(file.path(path, "HMP_MOCK.v35.fasta"))
	# match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, as.character(sread(mockRef))))))
	# cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
	
	seqtab.nochim <- read.table(f_seqtab_nochim_asis, header = T, sep = "\t", quote = "", fill = TRUE, comment.char = "", stringsAsFactors = F)
	seqtab.nochim.recv <- as.matrix(seqtab.nochim)
	
	taxa2.recv <- read.table(f_taxa_asis, header = T, sep = "\t", quote = "", fill = TRUE, comment.char = "", stringsAsFactors = F)
	taxa2.recv <- as.matrix(taxa2.recv)
	
	# Bonus: Handoff to phyloseq
	# Import into phyloseq:
	library(phyloseq); packageVersion("phyloseq")
	#[1] ‘1.16.2’
	library(ggplot2); packageVersion("ggplot2")
	#[1] ‘2.1.0’
	
	# Make a data.frame holding the sample data
	samples.out <- rownames(seqtab.nochim)
	
		# #F3D0_S188_L001_R1_001.fastq
		# #F3D5_S193_L001_R1_001.fastq
		# #F3D150_S216_L001_R1_001.fastq
	# subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
	# gender <- substr(subject,1,1)
	# subject <- substr(subject,2,999)
	# day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
	# samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
	# samdf$When <- "Early"
	# samdf$When[samdf$Day>100] <- "Late"
	# rownames(samdf) <- samples.out
	
	samples <- samples.out
		#Hak16S1 -> Hak16S1
		#Hak16S2 -> 	#Hak16S2 -> 2
	#samples <- sapply(strsplit(samples.out, "OD"), `[`, 2)
	#samples <- sapply(strsplit(samples.out, "16S"), `[`, 2)
		#Hak16S1 -> 1
		#Hak16S2 -> 2
	#duplicate <- c(1,2)
	samdf <- data.frame(Sample_name = samples)
	rownames(samdf) <- samples.out
	
	# Construct phyloseq object (straightforward from dada2 outputs)
	cat("Construct phyloseq object..\n")
	ps <- phyloseq(
		otu_table(seqtab.nochim.recv, taxa_are_rows = FALSE), 
		#otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
		sample_data(samdf), 
		tax_table(taxa2.recv)
		#tax_table(taxa2)
	)
	ps
	
	if (0){
	ps <- phyloseq(
		#otu_table(seqtab.nochim.recv, taxa_are_rows = FALSE), 
		otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
		sample_data(samdf), 
		#tax_table(taxa2.recv)
		tax_table(taxa2)
	)
	ps
	}
	
	# Visualize alpha-diversity:
	cat("Visualize alpha-diversity..\n")
	g0 <- plot_richness(ps, x = "Sample_name", measures = c("Shannon", "Simpson"), color = "Sample_name") + theme_bw()
	pdf(file.path(outfiltdir, "Rplot.richness.pdf"), onefile = T)
	print(g0)
	dev.off()
	
	# estimate_richness(physeq, split = TRUE, measures = measures) で: 
	# plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When") + theme_bw()
	 # estimate_richness(physeq, split = TRUE, measures = measures) で: 
	  # The data you have provided does not have
	# any singletons. This is highly suspicious. Results of richness
	# estimates (for example) are probably unreliable, or wrong, if you have already
	# trimmed low-abundance taxa from the data.
	
	# We recommended that you find the un-trimmed data and retry.
	
	# Ordinate:
	#ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")
	#plot_ordination(ps, ord.nmds.bray, color="Sample", title="Bray NMDS")
	#plot_ordination(ps, ord.nmds.bray, color="When", title="Bray NMDS")
	
	# Bar plot:
	cat("Bar plot..\n")
	top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
	ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
	ps.top20 <- prune_taxa(top20, ps.top20)
	plot_bar(ps.top20, x = "Sample_name", fill = "Family") #+ facet_wrap(~ Sample_name, scales = "free_x")
	#plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")
	
	g1 <- plot_bar(ps.top20, x = "Sample_name", fill = "Family") #+ facet_wrap(~ Sample_name, scales = "free_x")
	pdf(file.path(outfiltdir, "Rplot.barplot.family.top20.pdf"), onefile = T)
	print(g1)
	dev.off()
	
	g2 <- plot_bar(ps.top20, x = "Sample_name", fill = "Phylum") #+ facet_wrap(~ Sample_name, scales = "free_x")
	pdf(file.path(outfiltdir, "Rplot.barplot.phylum.top20.pdf"), onefile = T)
	print(g2)
	dev.off()
	
	all <- names(sort(taxa_sums(ps), decreasing=TRUE))
	ps.all <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
	ps.all <- prune_taxa(all, ps.all)
	g3 <- plot_bar(ps.all, x = "Sample_name", fill = "Phylum") #+ facet_wrap(~ Sample_name, scales = "free_x")
	pdf(file.path(outfiltdir, "Rplot.barplot.phylum.all.pdf"), onefile = T)
	print(g3)
	dev.off()
	
	return(obj)
}

x_____main_______________ <- function(){
}
#http://benjjneb.github.io/dada2/tutorial.html

# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install('dada2')
# BiocManager::install('phyloseq')
# 
# deprecated:
# source("http://bioconductor.org/biocLite.R")
# #Bioconductor version 3.3 (BiocInstaller 1.22.3), ?biocLite for help
# biocLite("dada2")
# biocLite("phyloseq")

# Getting ready

library(dada2); packageVersion("dada2")
#[1] '1.10.1'
library(ShortRead); packageVersion("ShortRead")
#[1] '1.40.0'
library(ggplot2); packageVersion("ggplot2")
#[1] ‘3.1.0’
library(plyr); packageVersion("plyr")
#[1] '1.8.4'

print(obj)


# if (!is.null(obj[['methods']])){
	# methods <- c('exec_qual_profile', 'exec_dada2', 'exec_plot_phyloseq')
# } else {
	# method_str <- obj[['methods']]
	# methods <- strsplit(method_str)[[0]]
# }

# for (method in methods){
	# list_method[[method]] <- 1
# }
# if ("x" %in% names(list_method))
# }

#tot_len = as.numeric(obj[['length']])
#dir <- obj[['dir']]
path <- obj[['dir']]
#path <- "~/Dropbox/work/160901_try_dada2_nancy_meta16S3/16S" # CHANGE ME to the directory containing the fastq files after unzipping.
if (
	is.null(obj[['outdir']])
){
	obj[['outdir']] <- path
}

outpath <- obj[['outdir']]
cat("Output under outdir", outpath, "...\n")
if(!file_test("-d", outpath)){
	cat(paste0("dir.create(outpath)", "\n"))
	dir.create(outpath)
}
#if(!file_test("-d", outpath)) dir.create(outpath)
#dir.create(outpath)
ext_fq <- obj[['ext']]

# Set below after f_qual_profile, to avoid error when these parameters are empty. 
#trimLeft <- as.numeric(strsplit(obj[['trimLefts']], ',')[[1]])
#truncLen <- as.numeric(strsplit(obj[['truncLens']], ',')[[1]])
#trimLeft <- c(17, 21)
#truncLen <- c(290, 215)

# Set below
#tax_file <- obj[['tax']]
#tax_file <- "/home/local/db/dada2/silva_nr_v123_train_set.fa.gz"
#tax_file <- "~/Dropbox/work/160901_try_dada2_tutorial/assignTaxonomy/silva_nr_v123_train_set.fa.gz"

cat('setwd', path, '...', "\n")
setwd(path)
# ln -s /Users/mickey/Dropbox/work/160901_try_dada2_tutorial/assinTaxonomy /Users/mickey/Dropbox/work/160901_try_dada2_nancy_meta16S 
cat("\n")

obj <- f_init_in_fq_files(obj)
obj <- f_init_sample_names(obj)

# # fns <- list.files(path)
# cat("File names under the given dir, ", path, "\n")
# fns
# cat("\n")


# fastqs <- fns[grepl(ext_fq, fns)]
# #fastqs <- fns[grepl(".fastq.gz$", fns)]
# fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
# fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
# fnRs <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files
# # Get sample names from the first part of the forward read filenames
# sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)

# cat("File names for forward reads (assuming given ext_fq and '_R1')\n")
# fnFs
# cat("\n")

# cat("File names for reverse reads (assuming given ext_fq and '_R2')\n")
# fnRs
# cat("\n")

# cat("Sample names taken from file names for forward reads)\n")
# sample.names
# cat("Note. filterAndTrim does not use above sample names, so samples names for filter count are independently set in f_summarize by the same way as done here.\n")
# cat("\n")

# # Fully specify the path for the fnFs and fnRs
# fnFs <- file.path(path, fnFs)
# fnRs <- file.path(path, fnRs)

# obj[['sample.names']] <- sample.names
# obj[['fnFs']] <- fnFs
# obj[['fnRs']] <- fnRs


# Examine quality profiles of forward and reverse reads
start_time <- Sys.time()
obj <- f_init_qual_profile(obj)
f_qual_profile_raw <- obj[['file_qual_profile_raw']]
if (file_test("-f", f_qual_profile_raw)){
	cat("file_qual_profile_raw found, nothing done", " ", "(", f_qual_profile_raw, ")", ".\n", sep = "")
}else{
	# Requirement: 
	# outpath <- obj[['outdir']]
	# sample.names <- obj[['sample.names']]
	# fnFs <- obj[['fnFs']]
	# fnRs <- obj[['fnRs']]
	
	cat("Not found: ", f_qual_profile_raw, "\n", "Do f_qual_profile..\n", sep = "")
	#cat("file_qual_profile_raw not found, do f_qual_profile..\n")
	obj <- f_qual_profile(obj)
}
#f_qual_profile(obj)
end_time <- Sys.time()
print(end_time - start_time)
# Time difference by minute
# equivalent to 
#print(difftime(end_time, start_time, units = c("auto")))

if (is.null(obj[['outfiltdir']])){	#can be set as argument
	obj <- f_init_outfiltdir(obj)
}

obj <- f_init_construct_sequence_table(obj)
	# f_seqtab <- obj[['file_seqtab']]
	# f_seqtab_fa <- obj[['file_seqtab_fa']]
	# f_seqtab_asis <- obj[['file_seqtab_asis']]
	# f_seqtab_otu <- obj[['file_seqtab_otu']]

cat(
	"Two quality-filtered files to check..\n", 
		"\t", "dir: ", dirname(obj[['file_qual_profile_filt']]), "\n",
		"\t", basename(obj[['file_qual_profile_filt']]), "\n",
		"\t", basename(obj[['file_qual_profile_filt_rds']]), "\n",
	sep = ""
)
if (
	file_test("-f", obj[['file_qual_profile_filt']])
		&&
	file_test("-f", obj[['file_qual_profile_filt_rds']])
){
	cat("All of the above two quality-filtered files found, nothing done. \n", "Instead, call f_init_from_filt_rds to set parameters (sample.names, baseFs, baseRs, fnFs, fnRs, filtFs, filtRs, filterAndTrim_my, filterAndTrim_my_zero) in case some samples have no reads after filtering.\n", "\n", sep = "")
	#cat("All of the above two quality-filtered files found, nothing done.\n", "\n", sep = "")
	
	obj <- f_init_filter(obj)
	obj <- f_init_from_filt_rds(obj)
	# if (is.null(obj[['filterAndTrim']])){
		# cat("filterAndTrim object is not set (this case occurs when f_qual_profile_filt_rds step is skipped), try to set it from file..\n")
		# print(file_qual_profile_filt_rds)
		# out <- readRDS(file_qual_profile_filt_rds)
		# #out0 <- readRDS(file_qual_profile_filt_rds)
		# #out <- out0[row.names(out0) %in% sample.names,]
		
		# obj[['filterAndTrim']] <- out
	# }
}else{
	cat("Either of the above two quality-filtered files not found, do f_filter..\n", "\n", sep = "")
	
	
	# Filtering and Trimming
	cat("f_init_filter_parameters..\n")
	obj <- f_init_filter_parameters(obj)
	
	# keys_trim <- c('trimLefts', 'truncLens')
	# for (key in keys_trim){
		# if (
			# is.null(obj[[key]])
		# ){
			# cat(msg, "\n")
			# cat(paste0("!!!\n", "key ", key, " is not set. \n", "!!!\n"))
			# quit()
		# }
	# }
	
	# trimLeft <- as.numeric(strsplit(obj[['trimLefts']], ',')[[1]])
	# truncLen <- as.numeric(strsplit(obj[['truncLens']], ',')[[1]])
	# trimLeft
	# truncLen
	# obj[['trimLeft']] <- trimLeft
	# obj[['truncLen']] <- truncLen
	
	# # Examine quality profiles of forward and reverse reads
	# obj <- f_init_qual_profile(obj)
	# f_qual_profile_raw <- obj[['f_qual_profile_raw']]
	# if (file_test("-f", f_qual_profile_raw)){
	# }else{
		# obj <- f_qual_profile(obj)
	# }
	# #f_qual_profile(obj)
	
	# Filter
	start_time <- Sys.time()
	
	cat("f_init_filter..\n", "\n", sep = "")
	obj <- f_init_filter(obj)
	
	file_qual_profile_filt_rds <- obj[['file_qual_profile_filt_rds']]
	if (file_test("-f", file_qual_profile_filt_rds)){
		cat("file_qual_profile_filt_rds found, nothing done.\n", "Instead, call f_init_from_filt_rds to set parameters (sample.names, baseFs, baseRs, fnFs, fnRs, filtFs, filtRs, filterAndTrim_my, filterAndTrim_my_zero) in case some samples have no reads after filtering.\n", sep = "")
		#cat("file_qual_profile_filt_rds found, nothing done.\n")
		
		obj <- f_init_from_filt_rds(obj)
		# if (is.null(obj[['filterAndTrim']])){
			# cat("filterAndTrim object is not set (this case occurs when f_qual_profile_filt_rds step is skipped), try to set it from file..\n")
			# print(file_qual_profile_filt_rds)
			# out <- readRDS(file_qual_profile_filt_rds)
			# #out0 <- readRDS(file_qual_profile_filt_rds)
			# #out <- out0[row.names(out0) %in% sample.names,]
			
			# obj[['filterAndTrim']] <- out
		# }
	} else {
		cat("file_qual_profile_filt_rds not found, do f_filter..\n")
		obj <- f_filter(obj)
	}
	
	file_qual_profile_filt <- obj[['file_qual_profile_filt']]
	if (file_test("-f", file_qual_profile_filt)){
		cat("file_qual_profile_filt found, nothing done.\n")
	} else {
		cat("file_qual_profile_filt not found, do f_plot_filter..\n")
		obj <- f_plot_filter(obj)
	}
	
	#obj <- f_filter(obj)
	end_time <- Sys.time()
	print(end_time - start_time)
	cat("\n")
}

cat(
	"Five seqtab files to check..\n", 
		"\t", "dir: ", dirname(obj[['file_seqtab']]), "\n",
		"\t", basename(obj[['file_seqtab']]), "\n",
		"\t", basename(obj[['file_seqtab_fa']]), "\n",
		"\t", basename(obj[['file_seqtab_asis']]), "\n",
		"\t", basename(obj[['file_seqtab_otu']]), "\n",
		"\t", basename(obj[['file_seqtab_rds']]), "\n",
	sep = ""
)
if (
	file_test("-f", obj[['file_seqtab']])
		&&
	file_test("-f", obj[['file_seqtab_fa']])
		&&
	file_test("-f", obj[['file_seqtab_asis']])
		&&
	file_test("-f", obj[['file_seqtab_otu']])
		&&
	file_test("-f", obj[['file_seqtab_rds']])
){
	cat("All of the above five seqtab files found, nothing done.\n", "\n", sep = "")
}else{
	cat("Either of the above five seqtab files not found, do f_qual_profile..\n", "\n", sep = "")
	
	# Filtering and Trimming
	cat("f_init_filter_parameters..\n")
	obj <- f_init_filter_parameters(obj)
	
	# keys_trim <- c('trimLefts', 'truncLens')
	# for (key in keys_trim){
		# if (
			# is.null(obj[[key]])
		# ){
			# cat(msg, "\n")
			# cat(paste0("!!!\n", "key ", key, " is not set. \n", "!!!\n"))
			# quit()
		# }
	# }
	
	# trimLeft <- as.numeric(strsplit(obj[['trimLefts']], ',')[[1]])
	# truncLen <- as.numeric(strsplit(obj[['truncLens']], ',')[[1]])
	# trimLeft
	# truncLen
	# obj[['trimLeft']] <- trimLeft
	# obj[['truncLen']] <- truncLen
	
	# # Examine quality profiles of forward and reverse reads
	# obj <- f_init_qual_profile(obj)
	# f_qual_profile_raw <- obj[['f_qual_profile_raw']]
	# if (file_test("-f", f_qual_profile_raw)){
	# }else{
		# obj <- f_qual_profile(obj)
	# }
	# #f_qual_profile(obj)
	
	# # Filter
	# start_time <- Sys.time()
	cat("f_init_filter..\n")
	obj <- f_init_filter(obj)
	# file_qual_profile_filt <- obj[['file_qual_profile_filt']]
	# if (file_test("-f", file_qual_profile_filt)){
		# cat("file_qual_profile_filt found, nothing done.\n")
	# } else {
		# cat("file_qual_profile_filt not found, do f_filter..\n")
		# obj <- f_filter(obj)
	# }
	# #obj <- f_filter(obj)
	# end_time <- Sys.time()
	# print(end_time - start_time)
	
	# Dereplication
	start_time <- Sys.time()
	file_derepFs_rds <- obj[['file_derepFs_rds']]
	file_derepRs_rds <- obj[['file_derepRs_rds']]
	if (
		file_test("-f", file_derepFs_rds)
			&&
		file_test("-f", file_derepRs_rds)
	){
		cat("file_derepFs_rds and file_derepRs_rds found, nothing done.\n")
		cat("Set derepFs and derepRs from the rds files.\n", sep = "")
		
		obj[['derepFs']] <- readRDS(file_derepFs_rds)
		obj[['derepRs']] <- readRDS(file_derepRs_rds)
	} else {
		cat("file_derepFs_rds or file_derepRs_rds not found, do f_derep..\n")
		obj <- f_derep(obj)
	}
	#cat("f_derep..\n")
	#obj <- f_derep(obj)
	end_time <- Sys.time()
	print(end_time - start_time)
	cat("\n")
	
	# Sample Inference
	# very slow
	start_time <- Sys.time()
	file_dadaFs_rds <- obj[['file_dadaFs_rds']]
	file_dadaRs_rds <- obj[['file_dadaRs_rds']]
	if (
		file_test("-f", file_dadaFs_rds)
			&&
		file_test("-f", file_dadaRs_rds)
	){
		cat("file_dadaFs_rds and file_dadaRs_rds found, nothing done.\n")
		cat("Set dadaFs and dadaRs from the rds files.\n", sep = "")
		
		obj[['dadaFs']] <- readRDS(file_dadaFs_rds)
		obj[['dadaRs']] <- readRDS(file_dadaRs_rds)
	} else {
		cat("file_dadaFs_rds or file_dadaRs_rds not found, do f_infer_sample..\n")
		obj <- f_infer_sample(obj)
	}
	#cat("f_infer_sample..\n")
	#obj <- f_infer_sample(obj)
	end_time <- Sys.time()
	print(end_time - start_time)
	cat("\n")
	
	# Merge paired reads
	# Merge the denoised forward and reverse reads:
	start_time <- Sys.time()
	file_mergers_rds <- obj[['file_mergers_rds']]
	if (file_test("-f", file_mergers_rds)){
		cat("file_mergers_rds found, nothing done.\n")
		cat("Set mergers from the rds file.\n", sep = "")
		
		obj[['mergers']] <- readRDS(file_mergers_rds)
	} else {
		cat("file_mergers_rds not found, do f_merge_paired_reads..\n")
		obj <- f_merge_paired_reads(obj)
	}
	#cat("f_merge_paired_reads..\n")
	#obj <- f_merge_paired_reads(obj)
	end_time <- Sys.time()
	print(end_time - start_time)
	cat("\n")
	
	# Constructing the sequence table
	start_time <- Sys.time()
	cat("f_construct_sequence_table..\n")
	obj <- f_construct_sequence_table(obj)
	end_time <- Sys.time()
	print(end_time - start_time)
	cat("\n")
}

obj <- f_init_path_remove_chimeras(obj)
	# f_seqtab_nochim <- obj[['file_seqtab_nochim']]
	# f_seqtab_nochim_fa <- obj[['file_seqtab_nochim_fa']]
	# f_seqtab_nochim_asis <- obj[['file_seqtab_nochim_asis']]
	# f_seqtab_nochim_otu <- obj[['file_seqtab_nochim_otu']]
	# f_summary_seqtab <- obj[['file_summary_seqtab']]
	# f_summary_seqlen <- obj[['file_summary_seqlen']]
	
cat(
	"Five seqtab_nochim files to check..\n", 
		"\t", "dir: ", dirname(obj[['file_seqtab']]), "\n",
		"\t", basename(obj[['file_seqtab_nochim']]), "\n",
		"\t", basename(obj[['file_seqtab_nochim_fa']]), "\n",
		"\t", basename(obj[['file_seqtab_nochim_asis']]), "\n",
		"\t", basename(obj[['file_seqtab_nochim_otu']]), "\n",
		"\t", basename(obj[['file_seqtab_nochim_rds']]), "\n",
	sep = ""
)

if (
	file_test("-f", obj[['file_seqtab_nochim']])
		&&
	file_test("-f", obj[['file_seqtab_nochim_fa']])
		&&
	file_test("-f", obj[['file_seqtab_nochim_asis']])
		&&
	file_test("-f", obj[['file_seqtab_nochim_otu']])
		&&
	file_test("-f", obj[['file_seqtab_nochim_rds']])
){
	cat("All of the above five seqtab_nochim files found, nothing done.\n", "\n", sep = "")
}else{
	cat("Either of the above five seqtab_nochim files not found, do f_remove_chimeras..\n", "\n", sep = "")
	# Remove chimeras
	# Set seqtab object if not (when f_construct_sequence_table step is skipped)
	obj <- f_init_remove_chimeras(obj)
	# outfiltdir <- obj[['outfiltdir']]
		# f_seqtab_nochim <- file.path(outfiltdir, "seqtab_nochim.txt")
		# f_seqtab_nochim_fa <- file.path(outfiltdir, "seqtab_nochim.fa")
		# f_seqtab_nochim_asis <- file.path(outfiltdir, "seqtab_nochim.asis.txt")
		# f_seqtab_nochim_otu <- file.path(outfiltdir, "seqtab_nochim.fa.otu_table.txt")
		# obj[['path_seqtab_nochim']] <- f_seqtab_nochim
		# obj[['path_seqtab_nochim_fa']] <- f_seqtab_nochim_fa
		# obj[['path_seqtab_nochim_asis']] <- f_seqtab_nochim_asis
		# obj[['path_seqtab_nochim_otu']] <- f_seqtab_nochim_otu
	# Remove chimeric sequences:
	start_time <- Sys.time()
	cat("f_remove_chimeras..\n")
	# very slow
	obj <- f_remove_chimeras(obj)
	end_time <- Sys.time()
	print(end_time - start_time)
	cat("\n")
}


cat(
	"Two summary files to check..\n", 
		"\t", "dir: ", dirname(obj[['file_summary_seqtab']]), "\n",
		"\t", basename(obj[['file_summary_seqtab']]), "\n",
		"\t", basename(obj[['file_summary_seqlen']]), "\n",
	sep = ""
)

if (
	file_test("-f", obj[['file_summary_seqtab']])
		&&
	file_test("-f", obj[['file_summary_seqlen']])
){
	cat("All of the above two summary files found, nothing done.\n", "\n", sep = "")
}else{
	cat("Either of the above two summary files not found, do f_summarize..\n", "\n", sep = "")
	# Remove chimeras
	# Set seqtab object if not (when f_construct_sequence_table step is skipped)
	obj <- f_init_summarize(obj)
	# outfiltdir <- obj[['outfiltdir']]
		# f_seqtab_nochim <- file.path(outfiltdir, "seqtab_nochim.txt")
		# f_seqtab_nochim_fa <- file.path(outfiltdir, "seqtab_nochim.fa")
		# f_seqtab_nochim_asis <- file.path(outfiltdir, "seqtab_nochim.asis.txt")
		# f_seqtab_nochim_otu <- file.path(outfiltdir, "seqtab_nochim.fa.otu_table.txt")
		# f_summary_seqtab <- file.path(outfiltdir, "summary.seqtab.txt")
		# f_summary_seqlen <- file.path(outfiltdir, "summary.seqlength.txt")
		# obj[['path_seqtab_nochim']] <- f_seqtab_nochim
		# obj[['path_seqtab_nochim_fa']] <- f_seqtab_nochim_fa
		# obj[['path_seqtab_nochim_asis']] <- f_seqtab_nochim_asis
		# obj[['path_seqtab_nochim_otu']] <- f_seqtab_nochim_otu
		# obj[['path_summary_seqtab']] <- f_summary_seqtab
		# obj[['path_summary_seqlen']] <- f_summary_seqlen
	# Remove chimeric sequences:
	start_time <- Sys.time()
	cat("f_summarize..\n")
	obj <- f_summarize(obj)
	end_time <- Sys.time()
	print(end_time - start_time)
	cat("\n")
}


# Assign taxonomy
# keys_tax <- c('tax')	#required in f_assign_taxonomy
# for (key in keys_trim){
	# if (
		# is.null(obj[[key]])
	# ){
		# cat(msg, "\n")
		# cat(paste0("!!!\n", "key ", key, " is not set. \n", "!!!\n"))
		# quit()
	# }
# }
obj <- f_init_path_assign_taxonomy(obj)
	# file_taxa <- obj[['file_taxa']]
	# file_taxa_w_count <- obj[['file_taxa_w_count']]
	# file_taxa_asis <- obj[['file_taxa_asis']]

cat(
	"Four taxa files to check..\n", 
		"\t", "dir: ", dirname(obj[['file_taxa']]), "\n",
		"\t", basename(obj[['file_taxa']]), "\n",
		"\t", basename(obj[['file_taxa_w_count']]), "\n",
		"\t", basename(obj[['file_taxa_asis']]), "\n",
		"\t", basename(obj[['file_taxa_rds']]), "\n",
	sep = ""
)

if (
	file_test("-f", obj[['file_taxa']])
		&&
	file_test("-f", obj[['file_taxa_w_count']])
		&&
	file_test("-f", obj[['file_taxa_asis']])
		&&
	file_test("-f", obj[['file_taxa_rds']])
){
	cat("All of the above four taxa files found, nothing done.\n", "\n", sep = "")
}else{
	cat("Either of the above four taxa files not found, do f_assign_taxonomy..\n", "\n", sep = "")
	
	# Set seqtab_nochim object if not (when f_remove_chimeras step is skipped)
	obj <- f_init_assign_taxonomy(obj)
	
	start_time <- Sys.time()
	cat("f_assign_taxonomy..\n")
	obj <- f_assign_taxonomy(obj)
	end_time <- Sys.time()
	print(end_time - start_time)
	cat("\n")
}

start_time <- Sys.time()
cat("f_phyloseq_plot..\n")
obj <- f_phyloseq_plot(obj)
end_time <- Sys.time()
print(end_time - start_time)


x_____For_copy_paste_______________ <- function(){
}

####
if (0){

setwd(path)
# ln -s /Users/mickey/Dropbox/work/160901_try_dada2_tutorial/assinTaxonomy /Users/mickey/Dropbox/work/160901_try_dada2_nancy_meta16S 
fns <- list.files(path)
fns

# Filtering and Trimming

fastqs <- fns[grepl(ext_fq, fns)]
#fastqs <- fns[grepl(".fastq.gz$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
fnRs <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files
# Get sample names from the first part of the forward read filenames
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Fully specify the path for the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

# Examine quality profiles of forward and reverse reads
pdf(file.path(outpath, "Rplot.qual_profile.raw.pdf"), onefile = T)
#pdf(file.path(path, "Rplot.qual_profile.raw.pdf"), onefile = T)
#pdf("Rplot.qual_profile.pdf", onefile = T)
for(i in seq_along(fnFs)) {
#for (i in c(1:4)){
	sample.name <- sample.names[[i]]
	cat(paste0(sample.name,'.. '))
	#cat(paste0(sample.name,'..',"\n"))
	
	gF <- plotQualityProfile(fnFs[[i]], n = 10000)
	gF <- gF + labs(title = paste0(sample.name, ' ', '(Fwd)'))
	gR <- plotQualityProfile(fnRs[[i]], n = 10000)
	gR <- gR + labs(title = paste0(sample.name, ' ', '(Rev)'))
	print(gF)
	print(gR)
	#print(plotQualityProfile(fnFs[[i]], n = 10000))
	#print(plotQualityProfile(fnRs[[i]], n = 10000))
	#print(plotQualityProfile(fnFs[[i]]))
	#print(plotQualityProfile(fnRs[[i]]))
}
dev.off()
# plotQualityProfile(fnFs[[1]])
# plotQualityProfile(fnFs[[2]])

# plotQualityProfile(fnRs[[1]])
# plotQualityProfile(fnRs[[2]])

ee <- 1
#ee <- 2
n <- 0
str_trimLeft <- paste(trimLeft, collapse = '_')
str_truncLen <- paste(truncLen, collapse = '_')
#str_trimLeft <- paste(trimLeft, collapse = ',')
#str_truncLen <- paste(truncLen, collapse = ',')
base <- paste("filtered", ".tleft", str_trimLeft, '.tlen', str_truncLen, ".ee", ee, ".tq2.n", n, sep = "")
#base <- paste("filtered", ".tleft16,20.tlen290,215.ee", ee, ".tq2.n", n, sep = "")
outfiltdir <- file.path(outpath, base)
#outfiltdir <- file.path(path, base)
if(!file_test("-d", outfiltdir)) dir.create(outfiltdir)

path_filt_fq <- file.path(outfiltdir, 'filt.fq')
if(!file_test("-d", path_filt_fq)) dir.create(path_filt_fq)
filtFs <- file.path(path_filt_fq, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path_filt_fq, paste0(sample.names, "_R_filt.fastq.gz"))

# Filter
for(i in seq_along(fnFs)) {
	sample.name <- sample.names[[i]]
	cat(paste0(sample.name,'.. '))
	
	fastqPairedFilter(
		c(fnFs[i], fnRs[i]), 
		c(filtFs[i], filtRs[i]),
		trimLeft = trimLeft,
		truncLen = truncLen, 
		#trimLeft = c(16, 20), 
		#truncLen = c(290,215), 
		maxN = n, 
		maxEE = ee, 
		truncQ = 2, 
		compress = TRUE, verbose = TRUE
	)
}

pdf(file.path(outfiltdir, "Rplot.qual_profile.filt.pdf"), onefile = T)
#pdf("Rplot.qual_profile.pdf", onefile = T)
for(i in seq_along(fnFs)) {
#for (i in c(1:4)){
	sample.name <- sample.names[[i]]
	cat(paste0(sample.name,'.. '))
	
	gF <- plotQualityProfile(filtFs[[i]], n = 10000)
	gF <- gF + labs(title = paste0(sample.name, ' ', '(Fwd)'))
	gR <- plotQualityProfile(filtFs[[i]], n = 10000)
	gR <- gR + labs(title = paste0(sample.name, ' ', '(Rev)'))
	print(gF)
	print(gR)
	#print(plotQualityProfile(filtFs[[i]], n = 10000))
	#print(plotQualityProfile(filtRs[[i]], n = 10000))
	#print(plotQualityProfile(filtFs[[i]]))
	#print(plotQualityProfile(filtRs[[i]]))
}
dev.off()

# Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

f_summary_derep <- file.path(outfiltdir, "summary.derepFastq.txt")
#f_summary_derep <- "summary.derepFastq.txt"
for(i in seq_along(fnFs)) {
	strF <- summary(derepFs[[i]])
	strR <- summary(derepRs[[i]])
	write(strF, file = f_summary_derep, append = T)
	write(strR, file = f_summary_derep, append = T)
}
#derepFs[[1]]
# derep-class: R object describing dereplicated sequencing reads
# $uniques: 7124 reads in 1860 unique sequences
  # Sequence lengths: min=230, median=230, max=230
# $quals: Quality matrix dimension:  1860 230
  # Consensus quality scores: min=12, median=37.75, max=39
# $map: Map from reads to unique sequences:  2 49 1 62 1 ...

path_rda_derep <- file.path(outfiltdir, 'rda.derep')
if(!file_test("-d", path_rda_derep)) dir.create(path_rda_derep)
for(i in seq_along(fnFs)) {
	f_obj_F <- file.path(path_rda_derep, paste0(sample.names[i], ".derepF.rda"))
	saveRDS(derepFs[[i]], file = f_obj_F)
	f_obj_R <- file.path(path_rda_derep, paste0(sample.names[i], ".derepR.rda"))
	saveRDS(derepRs[[i]], file = f_obj_R)
}

# Sample Inference
# very slow
dadaFs <- dada(derepFs, err=NULL, selfConsist = TRUE)
dadaRs <- dada(derepRs, err=NULL, selfConsist = TRUE)
dadaFs[[1]]

path_rda_dada2 <- file.path(outfiltdir, 'rda.dada2')
if(!file_test("-d", path_rda_dada2)) dir.create(path_rda_dada2)
for(i in seq_along(fnFs)) {
	f_obj_F <- file.path(path_rda_dada2, paste0(sample.names[i], ".dadaF.rda"))
	saveRDS(dadaFs[[i]], file = f_obj_F)
	f_obj_R <- file.path(path_rda_dada2, paste0(sample.names[i], ".dadaR.rda"))
	saveRDS(dadaRs[[i]], file = f_obj_R)
}
#save(mod, file = "mymodel.rda")
#saveRDS

# Visualize estimated error rates:
#plotErrors(dadaFs[[1]], nominalQ=TRUE)
#plotErrors(dadaFs[[2]], nominalQ=TRUE)

pdf(file.path(outfiltdir, "Rplot.estimated_error_rates.pdf"), onefile = T)
#for(i in seq_along(fnFs)) {
plotErrors(dadaFs[[1]], nominalQ=TRUE)
plotErrors(dadaFs[[2]], nominalQ=TRUE)
#}
dev.off()

# Merge paired reads
# Merge the denoised forward and reverse reads:
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

path_mergers <- file.path(outfiltdir, 'mergers')
if(!file_test("-d", path_mergers)) dir.create(path_mergers)

#sample.names 
#samples <- names(mergers)
for(i in seq_along(fnFs)) {
	sample <- sample.names[i]
	names <- rownames(mergers[[i]])
	mergers1 <- cbind(name = c(names), mergers[[i]][c(2:length(mergers[[i]][1,]),1)])
	f_mergers <- file.path(path_mergers, paste("mergers", ".", sample, ".txt", sep = ""))
	write.table(mergers1, file = f_mergers, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

# Constructing the sequence table
seqtab <- makeSequenceTable(mergers)
#seqtab <- makeSequenceTable(mergers[names(mergers) != "Mock"])
dim(seqtab)
# [1]    4 2521
# Inspect distribution of sequence lengths
table(nchar(colnames(seqtab)))
# 290 294 298 299 312 318 335 351 363 372 375 387 406 407 409 411 412 413 414 
  # 1   1  27   2   1   1   1   1   1   1   1   1   1  65  20 502  34 232 224 
# 415 438 440 441 442 443 444 450 451 452 454 455 459 463 465 468 
 # 42 148 483 237 109   5   2 180  21  81   1  86   1   1   4   3

t.seqtab <- t(seqtab)
t.seqtab1 <- as.data.frame(cbind(t.seqtab, seq = rownames(t.seqtab)))
f_seqtab <- file.path(outfiltdir, "seqtab.txt")
write.table(t.seqtab1, file = f_seqtab, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#labels <- paste('dada2', seq(0,length(t.seqtab[,1])-1), sep = "_")
labels <- paste0("dada2_", seq(0, length(getUniques(seqtab))-1), ";size=", unname(getUniques(seqtab)), ";")
f_seqtab_fa <- file.path(outfiltdir, "seqtab.fa")
uniquesToFasta(seqtab, f_seqtab_fa, ids = labels, mode = "w")

f_seqtab_asis <- file.path(outfiltdir, "seqtab.asis.txt")
write.table(seqtab, file = f_seqtab_asis, append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
#seqtab.recv <- read.table(f_seqtab_asis, header = T, sep = "\t", quote = "", fill = TRUE, comment.char = "")

#t.seqtab.nochim <- t(seqtab.nochim)
t.seqtab.otu <- as.data.frame(cbind('#OTU ID' = labels, t.seqtab))
f_seqtab_otu <- file.path(outfiltdir, "seqtab.fa.otu_table.txt")
write.table(t.seqtab.otu, file = f_seqtab_otu, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Remove chimeras
# Remove chimeric sequences:
# very slow
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)
dim(seqtab.nochim)
# [1]   4 162

sum(seqtab.nochim)/sum(seqtab)
# [1] 0.1822704

table(nchar(colnames(seqtab.nochim)))
# 290 294 298 299 312 318 335 351 363 372 375 387 406 407 409 411 412 413 414 
  # 1   1   1   2   1   1   1   1   1   1   1   1   1   3   1  26  10  24  10 
# 415 438 440 441 442 443 444 450 451 452 454 455 459 463 465 468 
  # 2   7  19  11   8   2   2   7   1   2   1   4   1   1   3   3 

t.seqtab.nochim <- t(seqtab.nochim)
t.seqtab.nochim1 <- as.data.frame(cbind(t.seqtab.nochim, seq = rownames(t.seqtab.nochim)))
f_seqtab_nochim <- file.path(outfiltdir, "seqtab_nochim.txt")
write.table(t.seqtab.nochim1, file = f_seqtab_nochim, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#f_seqtab_nochim_fa <- file.path(outfiltdir, "seqtab_nochim.uchime.fa")
#uniquesToFasta(seqtab.nochim, f_seqtab_nochim_fa, ids = NULL, mode = "w")

#labels <- paste('dada2', seq(0,length(t.seqtab.nochim[,1])-1), sep = "_")
labels <- paste0("dada2_", seq(0, length(getUniques(seqtab.nochim))-1), ";size=", unname(getUniques(seqtab.nochim)), ";")
f_seqtab_nochim_fa <- file.path(outfiltdir, "seqtab_nochim.fa")
uniquesToFasta(seqtab.nochim, f_seqtab_nochim_fa, ids = labels, mode = "w")

f_seqtab_nochim_asis <- file.path(outfiltdir, "seqtab_nochim.asis.txt")
write.table(seqtab.nochim, file = f_seqtab_nochim_asis, append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
#seqtab.nochim.recv <- read.table(f_seqtab_nochim_asis, header = T, sep = "\t", quote = "", fill = TRUE, comment.char = "")

#t.seqtab.nochim <- t(seqtab.nochim)
t.seqtab.nochim.otu <- as.data.frame(cbind('#OTU ID' = labels, t.seqtab.nochim))
f_seqtab_nochim_otu <- file.path(outfiltdir, "seqtab_nochim.fa.otu_table.txt")
write.table(t.seqtab.nochim.otu, file = f_seqtab_nochim_otu, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

f_summary_seqtab <- file.path(outfiltdir, "summary.seqtab.txt")
es <- c()
e1 <- c(name = 'seqtab', nseq = sum(seqtab), nrep = dim(seqtab)[2])
e2 <- c(name = 'seqtab.nochim', nseq = sum(seqtab.nochim), nrep = dim(seqtab.nochim)[2])
es <- rbind(es, e1,e2)
write.table(es, file = f_summary_seqtab, append = F, quote = F, sep = "\t", row.names = F, col.names = T)

f_summary_seqlen <- file.path(outfiltdir, "summary.seqlength.txt")
con21 <- table(nchar(colnames(seqtab)))
es21 <- as.data.frame(con21)

con22 <- table(nchar(colnames(seqtab.nochim)))
es22 <- as.data.frame(con22)
es2 <- plyr::join(es21,es22,by=c('Var1'))
colnames(es2) <- c('length','seqtab','seqtab.nochim')
write.table(es2, file = f_summary_seqlen, append = F, quote = F, sep = "\t", row.names = F, col.names = T)

# Assign taxonomy

taxa2 <- assignTaxonomy(seqtab.nochim, tax_file, minBoot = 90)
#http://benjjneb.github.io/dada2/assign.html
colnames(taxa2) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
unname(head(taxa2))

taxa22 <- as.data.frame(cbind(taxa2, seq = rownames(taxa2)))
f_taxa2 <- file.path(outfiltdir, "seqtab_nochim.silva_nr_v123_train_set.fa.txt")
write.table(taxa22, file = f_taxa2, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

taxa_w_count <- plyr::join(t.seqtab.nochim1, taxa22, by = c('seq'))
taxa_w_count2 <- taxa_w_count[,c(1:length(sample.names),c(4:length(taxa_w_count[1,])),3)]
f_taxa_w_count <- file.path(outfiltdir, "seqtab_nochim.silva_nr_v123_train_set.fa.w_count.txt")
write.table(taxa_w_count2, file = f_taxa_w_count, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

f_taxa2_asis <- file.path(outfiltdir, "seqtab_nochim.silva_nr_v123_train_set.fa.asis.txt")
write.table(taxa2, file = f_taxa2_asis, append = FALSE, quote = FALSE, sep = "\t", row.names = T, col.names = TRUE)

# Evaluate accuracy
# DADA2 accuracy on mock community:
# unqs.mock <- getUniques(removeBimeraDenovo(mergers[["Mock"]], verbose=TRUE))

# cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

# mockRef <- readFasta(file.path(path, "HMP_MOCK.v35.fasta"))
# match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, as.character(sread(mockRef))))))
# cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")

seqtab.nochim.recv <- read.table(f_seqtab_nochim_asis, header = T, sep = "\t", quote = "", fill = TRUE, comment.char = "", stringsAsFactors = F)
seqtab.nochim.recv <- as.matrix(seqtab.nochim.recv)

taxa2.recv <- read.table(f_taxa2_asis, header = T, sep = "\t", quote = "", fill = TRUE, comment.char = "", stringsAsFactors = F)
taxa2.recv <- as.matrix(taxa2.recv)

# Bonus: Handoff to phyloseq
# Import into phyloseq:
library(phyloseq); packageVersion("phyloseq")
#[1] ‘1.16.2’
library(ggplot2); packageVersion("ggplot2")
#[1] ‘2.1.0’

# Make a data.frame holding the sample data
samples.out <- rownames(seqtab.nochim)

	# #F3D0_S188_L001_R1_001.fastq
	# #F3D5_S193_L001_R1_001.fastq
	# #F3D150_S216_L001_R1_001.fastq
# subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
# gender <- substr(subject,1,1)
# subject <- substr(subject,2,999)
# day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
# samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
# samdf$When <- "Early"
# samdf$When[samdf$Day>100] <- "Late"
# rownames(samdf) <- samples.out

samples <- samples.out
	#Hak16S1 -> Hak16S1
	#Hak16S2 -> 	#Hak16S2 -> 2
#samples <- sapply(strsplit(samples.out, "OD"), `[`, 2)
#samples <- sapply(strsplit(samples.out, "16S"), `[`, 2)
	#Hak16S1 -> 1
	#Hak16S2 -> 2
#duplicate <- c(1,2)
samdf <- data.frame(Sample = samples)
rownames(samdf) <- samples.out

# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(
	#otu_table(seqtab.nochim.recv, taxa_are_rows = FALSE), 
	otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
	sample_data(samdf), 
	#tax_table(taxa2.recv)
	tax_table(taxa2)
)
ps

if (0){
ps <- phyloseq(
	#otu_table(seqtab.nochim.recv, taxa_are_rows = FALSE), 
	otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
	sample_data(samdf), 
	#tax_table(taxa2.recv)
	tax_table(taxa2)
)
ps
}

# Visualize alpha-diversity:
pdf(file.path(outfiltdir, "Rplot.richness.pdf"), onefile = T)
plot_richness(ps, x = "Sample", measures = c("Shannon", "Simpson"), color = "Sample") + theme_bw()
dev.off()

# estimate_richness(physeq, split = TRUE, measures = measures) で: 
# plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When") + theme_bw()
 # estimate_richness(physeq, split = TRUE, measures = measures) で: 
  # The data you have provided does not have
# any singletons. This is highly suspicious. Results of richness
# estimates (for example) are probably unreliable, or wrong, if you have already
# trimmed low-abundance taxa from the data.

# We recommended that you find the un-trimmed data and retry.

# Ordinate:
#ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")
#plot_ordination(ps, ord.nmds.bray, color="Sample", title="Bray NMDS")
#plot_ordination(ps, ord.nmds.bray, color="When", title="Bray NMDS")

# Bar plot:
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Sample", fill="Family") + facet_wrap(~ Sample, scales = "free_x")
#plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")

pdf(file.path(outfiltdir, "Rplot.barplot.family.top20.pdf"), onefile = T)
plot_bar(ps.top20, x = "Sample", fill = "Family") + facet_wrap(~ Sample, scales = "free_x")
dev.off()

pdf(file.path(outfiltdir, "Rplot.barplot.phylum.top20.pdf"), onefile = T)
plot_bar(ps.top20, x = "Sample", fill = "Phylum") + facet_wrap(~ Sample, scales = "free_x")
dev.off()

all <- names(sort(taxa_sums(ps), decreasing=TRUE))
ps.all <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.all <- prune_taxa(all, ps.all)
pdf(file.path(outfiltdir, "Rplot.barplot.phylum.all.pdf"), onefile = T)
plot_bar(ps.all, x = "Sample", fill = "Phylum") + facet_wrap(~ Sample, scales = "free_x")
dev.off()

## w.o. facet_wrap
pdf(file.path(outfiltdir, "Rplot.barplot.family.top20.all.pdf"), onefile = T)
plot_bar(ps.top20, x = "Sample", fill = "Family") #+ facet_wrap(~ Sample, scales = "free_x")
dev.off()

pdf(file.path(outfiltdir, "Rplot.barplot.phylum.top20.all.pdf"), onefile = T)
plot_bar(ps.top20, x = "Sample", fill = "Phylum") #+ facet_wrap(~ Sample, scales = "free_x")
dev.off()

pdf(file.path(outfiltdir, "Rplot.barplot.phylum.all.all.pdf"), onefile = T)
plot_bar(ps.all, x = "Sample", fill = "Phylum") #+ facet_wrap(~ Sample, scales = "free_x")
dev.off()

pdf(file.path(outfiltdir, "Rplot.barplot.phylum.top20.all.log.pdf"), onefile = T)
plot_bar(ps.top20, x = "Sample", fill = "Phylum") + scale_y_log10() #+ facet_wrap(~ Sample, scales = "free_x")
dev.off()

}
