library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(data.table)

data(Manifest); head(Manifest)
data(Locations); head(Locations)
annotation <- cbind(Locations, Manifest); rm(Manifest, Locations); gc()
sum(!(rownames(annotation) == annotation$Name)) # 0

###################################### CpG ####################################

setwd("/home/ben/Documents/Git/UMtools_dat/dat/vcf/proc/")
CpG_SNP = fread("CpG_SNP.vcf")
colnames(CpG_SNP) = c("chr", "start", "rs", "M", "m", "v6", "v7", "metadata", "chr2", "start2", "end2", "cg")
CpG_SNP$real_pos = annotation[CpG_SNP$cg, "pos"]
deltaL = CpG_SNP$start - CpG_SNP$real_pos
table(deltaL, CpG_SNP$M)
# deltaL    C    G
# 0      7716    4
# 1         0 9004

###################################### SBE ####################################

setwd("/home/ben/Documents/Git/UMtools_dat/dat/vcf/proc/")
SBE_SNP = fread("SBE_SNP.vcf")
colnames(SBE_SNP) = c("chr", "start", "rs", "M", "m", "v6", "v7", "metadata", "chr2", "start2", "end2", "cg")
SBE_SNP$strand = annotation[SBE_SNP$cg,"strand"]
SBE_SNP$real_pos = annotation[SBE_SNP$cg, "pos"]
SBE_SNP_p = SBE_SNP[SBE_SNP$strand == "+",]
SBE_SNP_m = SBE_SNP[SBE_SNP$strand == "-",]
deltaL_p = SBE_SNP_p$start - SBE_SNP_p$real_pos
deltaL_m = SBE_SNP_m$start - SBE_SNP_m$real_pos

table(deltaL_p)
#  -1
# 281

table(deltaL_m)
#    2
#  281

###################################### SBE ####################################

setwd("/home/ben/Documents/Git/UMtools_dat/dat/vcf/proc/")
probe_SNP = fread("probe_SNP.vcf")
colnames(probe_SNP) = c("chr", "start", "rs", "M", "m", "v6", "v7", "metadata", "chr2", "start2", "end2", "cg")
probe_SNP$strand = annotation[probe_SNP$cg, "strand"]
probe_SNP$real_pos = annotation[probe_SNP$cg, "pos"]
probe_SNP$Type = annotation[probe_SNP$cg, "Type"]
probe_SNP_p = probe_SNP[probe_SNP$strand == "+",]
probe_SNP_m = probe_SNP[probe_SNP$strand == "-",]
probe_SNP_p_I = probe_SNP_p[probe_SNP_p$Type == "I",]
probe_SNP_p_II = probe_SNP_p[probe_SNP_p$Type == "II",]
probe_SNP_m_I = probe_SNP_m[probe_SNP_m$Type == "I",]
probe_SNP_m_II = probe_SNP_m[probe_SNP_m$Type == "II",]


deltaL_p_I = probe_SNP_p_I$start - probe_SNP_p_I$real_pos
deltaL_p_II = probe_SNP_p_II$start - probe_SNP_p_II$real_pos
deltaL_m_I = probe_SNP_m_I$real_pos - probe_SNP_m_I$start
deltaL_m_II = probe_SNP_m_II$real_pos - probe_SNP_m_II$start

table(deltaL_p_I)
# 2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28
# 255 236 290 259 331 282 319 278 356 411 357 399 369 379 365 400 382 407 393 400 360 357 396 394 353 356 380
# 29  30  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49
# 387 316 366 354 354 396 390 360 419 366 398 361 359 387 366 386 368 361 374 381 358

table(deltaL_p_II)
# 2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28
# 435 509 536 515 540 544 517 531 707 731 749 730 717 748 764 790 797 777 742 689 807 733 723 734 758 743 778
# 29  30  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50
# 748 729 722 744 754 720 756 761 780 786 751 747 697 737 746 724 763 771 704 742 828 809

table(deltaL_m_I)
#   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27
# 237 289 320 321 292 288 313 309 348 301 408 390 390 386 392 400 406 372 392 404 386 379 398 375 383 380 376
# 28  29  30  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48
# 391 405 406 376 381 343 361 359 424 383 378 360 381 358 375 415 374 366 388 373 395

table(deltaL_m_II)
#   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27
# 492 516 499 509 544 517 528 564 521 483 706 708 668 756 722 740 738 735 704 742 726 740 735 716 712 704 717
# 28  29  30  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49
# 748 732 722 784 759 800 767 757 734 759 733 745 766 755 734 716 731 750 723 728 776 806


