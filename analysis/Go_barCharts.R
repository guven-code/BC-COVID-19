library(ggplot2)




term2=gsub("GO:[0-9]{7}~","",bp_fat10$Term)
bp_fat10$Term<-term2






ggsave("Goplots/downRegulated_david_kegg.pdf")#width=15,height = 10,units="cm")

options(repr.plot.width=6, repr.plot.height=3)
ggplot(GOtop_cc[,1:10], aes(x = Term, y = -log10(PValue), main="GO:CC of DEGs")) +
  geom_bar(stat = "identity") +
  coord_flip() + scale_y_continuous(name="-Log10(PValue)") +
  scale_x_discrete(name="Cellular Component") +
  theme(axis.text.x = element_text(face="bold", color="#008000",
                                   size=8, angle=0),
        axis.text.y = element_text(face="bold", color="#008000",
                                   size=8, angle=0))
dev.off()
#ggplot(bp_fat15, aes(x = Term,
#                     y = -log10(PValue),
 #                    fill = gene_Count)) + geom_bar(stat = "identity",position = "dodge")+ 
 # scale_color_gradient(low="blue",
      #                  high="red")+coord_flip()