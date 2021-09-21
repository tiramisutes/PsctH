p<-ggplot(species_tissues, aes(x=Species, y=n, fill=Species)) +
  geom_bar(stat="identity")+theme_minimal() + scale_fill_brewer(palette="Dark2")

p<-ggplot(cellTypes, aes(x=Species, y=n, fill=Species)) +
  geom_bar(stat="identity")+theme_minimal() + scale_fill_brewer(palette="Dark2")

library(stringr)
tissues_count$xlab2 <- str_wrap(tissues_count$xlab, width = 15)

ggplot(tissues_count, aes(x=xlab, y=n, fill=Tissue)) +
  geom_bar(stat="identity")+theme_minimal() + scale_fill_brewer(palette="Paired") +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))