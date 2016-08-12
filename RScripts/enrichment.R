#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
#print(args)
input_mean_file <- ''
ip_mean_file <- ''
input_id <- ''
ip_id <- ''
outdir <- '.'

ip_count <- 0 # as.integer(args$ipsz)
input_count <- 0 # as.integer(args$inputsz)

i <-1
while(i <= length(args)){
    if(args[i] == '--ip-mean-file'){
        i <- i+1
        ip_mean_file <- args[i]
    } else if(args[i] == '--ip-id'){
        i <- i+1
        ip_id <- args[i]
    } else if(args[i] == '--ip-count'){
        i <- i+1
        ip_count <- as.integer(args[i])
    } else if(args[i] == '--input-mean-file'){
        i <- i+1
        input_mean_file <- args[i]
    } else if(args[i] == '--input-id'){
        i <- i+1
        input_id <- args[i]
    } else if(args[i] == '--input-count'){
        i <- i+1
        input_count <- as.integer(args[i])
    } else if(args[i] == '--out'){
        i <- i+1
        outdir  <- args[i]
    }
    i  <- i+1
}

outdir <- paste(outdir,"/",sep="") ## make sure there is a / at the end

dir.create(outdir, showWarnings = FALSE)

txt_file <- paste(outdir,"results.txt",sep="")

x  <- cbind(
	read.delim(ip_mean_file, head=F, colClasses=c(rep("NULL",3), "numeric"), na.strings="nan"), 
	read.delim(input_mean_file, head=F, colClasses=c(rep("NULL",3), "numeric"), na.strings="nan"))

errors = rep(0,3)
x = x[order(x[,1], na.last=FALSE),]
x[is.na(x)] <- 0

######################
## Calculations
######################
n = nrow(x)

cs1 = cumsum(x[,1])
cs2 = cumsum(x[,2])
 
s1 = cs1[n]
s2 = cs2[n]

csdiff = cs2/s2 - cs1/s1;
#if(min(csdiff) < -1e-8) { cat(paste("Input enrichment stronger than IP at bin", which.min(csdiff), "\n")); abline(v = which.min(csdiff)/n, col= "orange", lty = 3, lwd = 2) } 
if(min(csdiff) < -1e-8) { errors[1] = which.min(csdiff) } 

k = which.max(csdiff)
p = cs1[k]/s1
q = cs2[k]/s2
#if(p == 0) {cat(paste("Zero-enriched IP, maximum difference at bin", k, "\n"))}
if(p == 0) {errors[2] = k}

k_pcr = n - floor(0.01 * n)

amplification_value <- cs2[k_pcr]/s2

######################
## Futher calculations
######################
div = -(p + q)/2 * log2((p + q)/2) - (1 - (p + q)/2) * log2(1 - (p + q)/2) + p/2 * log2(p) + (1 - p)/2 * log2(1 - p) + q/2 * log2(q) + (1 - q)/2 * log2(1 - q)
f = sqrt(div)
zscore = NA 

if(ip_count > 0 & input_count > 0) {
   d11 = log2(2 * p / (p + q)) / (4 * f)
   d21 = log2(2 * (1 - p) / (2 - p - q)) / (4 * f)
   d12 = log2(2 * q / (p + q)) / (4 * f)
   d22 = log2(2 * (1 - q) / (2 - p - q)) / (4 * f)
   varf = p * (1 - p) / ip_count * ( d11 * (d11 - d21) + d21 * (d21 - d11) ) + q * (1 - q) / input_count * (d21 * (d21 - d22) + d22 * (d22 - d21) )
   zscore = f/sqrt(varf)
}
pc_enrichment <- 100 * (n - k)/n
input_scaling_factor <- p / q * input_count / input_count
diff_pc_enrichment <- 100 * (q - p)

######################
## PLOT
######################
png_file <- paste(outdir,"enrichment.png",sep="")
png(file = png_file)
plot(y = cs1/s1, x = (1:n)/n, main = "Cumulative percentage enrichment in each channel", ylab = "Percentage of tags", xlab = "Percentage of bins", col = "darkblue", lwd = 3, ty="l")
lines(y = cs2/s2, x = (1:n)/n, col = "red", lwd = 3)
abline(v = k/n, col = "darkgreen", lty = 3, lwd = 2)

if(min(csdiff) < -1e-8) { 
	abline(v = which.min(csdiff)/n, col= "orange", lty = 3, lwd = 2) 
} 

if(amplification_value < 0.75){
#     cat(paste("PCR amplification bias in Input, coverage of 1% of genome", 1 - amplification_value, "\n")); 
     errors[3] <-  1 - amplification_value
     abline(h = amplification_value, col = "turquoise", lty = 3, lwd = 2) 
}
dev.off()

######################
## Collect results
######################
results = round(c(p, q, div, zscore, pc_enrichment, input_scaling_factor, diff_pc_enrichment), 4)
names(results) = c("p", "q", "divergence", "z_score", "percent_genome_enriched", "input_scaling_factor", "differential_percentage_enrichment")
#print(results)

results.df <- data.frame(type=names(results), result=sprintf("%.5f",results))

######################
## Write output file
######################
results.df <- rbind(results.df, data.frame(type=c("ip_id","input_id"), result=sprintf("%s",c(ip_id, input_id))))
results.df <- rbind(results.df, data.frame(type=c("error_1","error_2"), result=sprintf("%d", errors[1:2])))
results.df <- rbind(results.df, data.frame(type=c("error_3"), result=sprintf("%.4f", errors[3])))

print(results.df)

write.table(results.df,file=txt_file ,col.names=F, row.names=F, sep="\t", quote=F)

#FH <- file.create(txt_file)
#writeLines(txt,txt_file)
#close(FH)

#sink(txt_file, append=TRUE)
#cat(sprintf("%s	%s	%.5f	%.5f	%.5f	%.2f	%.3f	%.5f	%.4f	%d	%d	%.4f\n", c(ip_id, input_id, results, errors[1], errors[2], errors[3]))

#library("RSQLite")

#print (c(dbfile,outdir))

#tmp <- dbConnect(SQLite(), dbname=dbfile)
#summary(tmp)
#dbDisconnect(tmp)
