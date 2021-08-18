#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

non_control_region = args[1]
control_region_shifted = args[2]
out = args[3]

shift_back = function(x) {
  if (x < 8570) {
    return(x + 8000)
  } else {
    return (x - 8569)
  }
}

control_region_shifted = read.table(control_region_shifted, header=T)
shifted_back = sapply(control_region_shifted[,"pos"], shift_back)
control_region_shifted[,"pos"] = shifted_back

beginning = subset(control_region_shifted, control_region_shifted[,'pos']<8000)
end = subset(control_region_shifted, control_region_shifted[,'pos']>8000)

non_control_region = read.table(non_control_region, header=T)
combined_table = rbind(beginning, non_control_region, end)
write.table(combined_table, out, row.names=F, col.names=T, quote=F, sep="\t")
