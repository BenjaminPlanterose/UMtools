library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(data.table)
library(seqinr)

data(Manifest); head(Manifest)
data(Locations); head(Locations)
annotation <- cbind(Locations, Manifest); rm(Manifest, Locations); gc()
sum(!(rownames(annotation) == annotation$Name)) # 0
dim(annotation) # 865859     11

###################################### bed preparation ######################################

####### CpG sites
chr = as.vector(annotation$chr)
chr = unlist(lapply(strsplit(chr, split = "chr"), function(x) x[2]))
bed_1 = data.frame(chr = chr, start = as.integer(annotation$pos)-1L,
                   end = as.integer(annotation$pos)+1L, ID = rownames(annotation)) # (a, b]
dim(bed_1) # 865859      4
fwrite(bed_1, "CpG_sites.bed", quote = F, sep = "\t", col.names = F)


####### SBE sites
I = annotation[annotation$Type == "I",]; dim(I) # 135476     33
I_p = I[I$strand == "+",]; I_m = I[I$strand == "-",]
chr_p = unlist(lapply(strsplit(I_p$chr, split = "chr"), function(x) x[2]))
chr_m = unlist(lapply(strsplit(I_m$chr, split = "chr"), function(x) x[2]))

bed_p = data.frame(chr = chr_p, start = as.integer(I_p$pos)-2L, # (a,b]
                   end = as.integer(I_p$pos)-1L, ID = rownames(I_p))

bed_m = data.frame(chr = chr_m, start = as.integer(I_m$pos)+1L, # (a,b]
                   end = as.integer(I_m$pos)+2L, ID = rownames(I_m))
bed_2 = rbind(bed_p, bed_m); dim(bed_2) # 142137      4
fwrite(bed_2, "SBE_sites_typeI.bed", quote = F, sep = "\t", col.names = F)

####### Probe sites not CpG/SBE sites

# However, type-I probes include the CpG, type II are one nucleotide short
A_p = annotation[annotation$strand == "+",]; A_m = annotation[annotation$strand == "-",]
l_p = sapply(A_p$ProbeSeqA, function(x) length(s2c(x))); table(l_p) # All probes are 50 bp
l_m = sapply(A_m$ProbeSeqA, function(x) length(s2c(x))); table(l_m) # All probes are 50 bp
L = 50L
A_p_II = A_p[A_p$Type == "II",]; A_p_I = A_p[A_p$Type == "I",]
A_m_II = A_m[A_m$Type == "II",]; A_m_I = A_m[A_m$Type == "I",]
chr_p_I = unlist(lapply(strsplit(A_p_I$chr, split = "chr"), function(x) x[2]))
chr_p_II = unlist(lapply(strsplit(A_p_II$chr, split = "chr"), function(x) x[2]))
chr_m_I = unlist(lapply(strsplit(A_m_I$chr, split = "chr"), function(x) x[2]))
chr_m_II = unlist(lapply(strsplit(A_m_II$chr, split = "chr"), function(x) x[2]))
bed_m_I = data.frame(chr = chr_m_I, start = as.integer(A_m_I$pos)-L+1L,
                   end = as.integer(A_m_I$pos)-1L, ID = rownames(A_m_I)) # # (a, b]
bed_m_II = data.frame(chr = chr_m_II, start = as.integer(A_m_II$pos)-L,
                   end = as.integer(A_m_II$pos)-1L, ID = rownames(A_m_II)) # # (a, b]
bed_p_I = data.frame(chr = chr_p_I, start = as.integer(A_p_I$pos) + 1L,
                     end = as.integer(A_p_I$pos) + L - 1L, ID = rownames(A_p_I)) # # (a, b]
bed_p_II = data.frame(chr = chr_p_II, start = as.integer(A_p_II$pos) + 1L,
                      end = as.integer(A_p_II$pos) + L, ID = rownames(A_p_II)) # # (a, b]

bed_3 = Reduce(rbind, list(bed_p_I, bed_p_II, bed_m_I, bed_m_II)); dim(bed_3) # 485512      4
fwrite(bed_3, "probe_sites_notCpG.bed", quote = F, sep = "\t", col.names = F)



