library(ggplot2)


bp_fat=david_allDEGs_bpFat

term2=gsub("GO:[0-9]{7}~","",bp_fat$Term)
bp_fat$Term<-term2




data=bp_fat[,1:10]

ggsave("plots/GOplots/bp_david.pdf",width=15,height = 10,units="cm")


options(repr.plot.width=8, repr.plot.height=3)
ggplot(data, aes(x = Term, y = -log10(PValue), main="GO:CC of DEGs")) +
  geom_bar(stat = "identity") +
  coord_flip() + scale_y_continuous(name="-Log10(PValue)") +
  scale_x_discrete(name="Cellular Component") +
  theme(axis.text.x = element_text(face="bold", color="blue",#"#008000",
                                   size=8, angle=0),
        axis.text.y = element_text(face="bold", color="blue",#"#008000",
                                   size=8, angle=0))

options(repr.plot.width=8, repr.plot.height=3)
ggplot(y, aes(x = start_station_name, y = duration, main="Car Distribution")) +
  geom_bar(stat = "identity") +
  coord_flip() + scale_y_continuous(name="Average Trip Duration (in seconds)") +
  scale_x_discrete(name="Start Station") +
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