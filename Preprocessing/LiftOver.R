suppressMessages(library(rtracklayer))
suppressMessages(library(stringr))
suppressMessages(library(data.table))

args <- commandArgs(T)

InputFile <- args[1]
ChainFile <- args[2]
out_prefix <- args[3]

chainObject <- import.chain(ChainFile)
input <- fread(file = InputFile,header = T,stringsAsFactors = F)
names(input) <- tolower(names(input))
input$id <- paste(input$chromosome,input$position,input$effect_allele,input$other_allele,sep = ':')

index <- duplicated(input$id)
input <- input[!index,]

grObject <- GRanges(seqnames=paste0("chr",input$chromosome), ranges=IRanges(start=input$position, end=input$position),id=input$id)

results <- as.data.frame(liftOver(grObject, chainObject))

index <- match(results$id , input$id)

output <- cbind.data.frame(input[index,],results$start)

l <- dim(output)[2]

names(output)[l] <- "Position.Liftover"

write.table(output,file = paste0(out_prefix , ".liftover.txt"),row.names = F,col.names = T,quote = F,sep = '\t')
