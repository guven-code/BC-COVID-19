library(ggplot2)

david_bp_fat=kegg_david
david_bp_fat$gene_Count=david_bp_fat$Count
bp_fat10=david_bp_fat

term2=gsub("GO:[0-9]{7}~","",bp_fat10$Term)
bp_fat10$Term<-term2




library(wesanderson)

pal<-wes_palette("Zissou1",100,type='continuous')

ggsave("Goplots/downRegulated_david_kegg.pdf",width=15,height = 10,units="cm")
ggplot(bp_fat10, aes(x = Term,
                 y = -log10(PValue),
                 fill = gene_Count)) + geom_bar(stat = "identity",width=0.25,position = "dodge")+ 
                  scale_fill_gradientn(colours=pal)+scale_x_discrete(expand=c(0,0))+
scale_y_discrete(expand=c(0,0))+coord_flip()+theme(aspect.ratio = 2/1)#+ coord_flip()
dev.off()
#ggplot(bp_fat15, aes(x = Term,
#                     y = -log10(PValue),
 #                    fill = gene_Count)) + geom_bar(stat = "identity",position = "dodge")+ 
 # scale_color_gradient(low="blue",
      #                  high="red")+coord_flip()