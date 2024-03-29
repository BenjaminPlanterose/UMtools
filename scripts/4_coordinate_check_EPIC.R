library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(data.table)

data(Manifest); head(Manifest)
data(Locations); head(Locations)
annotation <- cbind(Locations, Manifest); rm(Manifest, Locations); gc()
sum(!(rownames(annotation) == annotation$Name)) # 0
dim(annotation) # 865859     11

###################################### CpG ####################################

setwd("/home/ben/Documents/Git/UMtools_dat/dat_EPIC/vcf/proc/")
CpG_SNP = fread("CpG_SNP.vcf")
colnames(CpG_SNP) = c("chr", "start", "rs", "M", "m", "v6", "v7", "metadata", "chr2", "start2", "end2", "cg")
CpG_SNP$real_pos = annotation[CpG_SNP$cg, "pos"]
deltaL = CpG_SNP$start - CpG_SNP$real_pos
table(deltaL, CpG_SNP$M)
# deltaL     C     G
# 0       14005     4
# 1           0 15182

###################################### SBE ####################################

setwd("/home/ben/Documents/Git/UMtools_dat/dat_EPIC/vcf/proc/")
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
# 277

table(deltaL_m)
#    2
#  279

###################################### SBE ####################################

setwd("/home/ben/Documents/Git/UMtools_dat/dat_EPIC/vcf/proc/")
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
#   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28
# 257 245 289 271 343 290 347 316 385 423 383 405 386 393 386 394 380 417 401 409 379 363 405 403 368 362 390
# 29  30  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49
# 387 321 377 383 357 411 397 383 425 383 397 390 389 413 363 393 395 376 386 374 362

table(deltaL_p_II)
# 2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22
# 949 1055 1101 1071 1151 1302 1286 1329 1485 1577 1481 1513 1483 1497 1536 1574 1554 1469 1533 1457 1548
# 23   24   25   26   27   28   29   30   31   32   33   34   35   36   37   38   39   40   41   42   43
# 1474 1411 1464 1489 1484 1537 1492 1500 1475 1541 1513 1514 1526 1532 1552 1538 1562 1543 1445 1440 1520
# 44   45   46   47   48   49   50
# 1435 1493 1554 1473 1485 1530 1571

table(deltaL_m_I)
#   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27
# 236 279 311 330 295 304 344 335 351 317 402 414 378 404 401 422 417 401 406 425 383 400 417 404 399 398 396
# 28  29  30  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48
# 386 418 419 392 391 355 363 373 434 392 388 373 385 378 382 418 382 375 401 390 391

table(deltaL_m_II)
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21
# 977 1049 1035 1115 1094 1256 1319 1371 1317 1308 1474 1490 1428 1579 1517 1487 1497 1467 1482 1457 1472
# 22   23   24   25   26   27   28   29   30   31   32   33   34   35   36   37   38   39   40   41   42
# 1506 1501 1487 1479 1447 1462 1467 1500 1480 1556 1507 1562 1574 1476 1474 1538 1489 1531 1538 1514 1474
# 43   44   45   46   47   48   49
# 1494 1518 1466 1495 1476 1504 1566


