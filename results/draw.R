require(ggplot2);require(reshape2);require(scales)

s = read.csv('s100_all_branches.csv')
head(s)
unique(s$Method)

variants = s$Method %in% c("CASTLES (Integrals)", "CASTLES (Taylor 33, Coal)" ,
                    "CASTLES (Taylor 33)" , "CASTLES (Taylor 34, Coal)" ,
                    "CASTLES (Taylor 35, Coal)") 
s[s$Method=="CASTLES (Taylor 33, Handle OG)","Method"]="CASTLES"

ggplot(aes(x=l.true,y=l.est,color=Branch.Type,linetype=Branch.Type),
       data=s[!variants,])+
  facet_grid(sub("_non","bp",sub("fasttree_genetrees_","",Condition))~Method)+
  geom_point(alpha=0.05)+
  scale_x_continuous(trans="log10",name="True length")+
  scale_y_continuous(trans="log10",name="Estimated length")+
  stat_smooth(color="grey30",se=F)+
  scale_color_brewer(palette = "Set2")+
  coord_cartesian(xlim=c(10^-4.4,1),ylim=c(10^-4.4,1))+
  geom_abline(color="blue",linetype=3)+
  theme_bw()+
  theme(legend.position = c(.24,.1))
ggsave("S100-correlation.png",width=12.5,height = 9)


ggplot(aes(x=sub("+","\n",Method,fixed=T),
           y=l.true-l.est,color=Branch.Type),
       data=s[!variants,])+
  facet_wrap(~sub("_non","bp",sub("fasttree_genetrees_","",Condition)),ncol=1)+
  #scale_x_continuous(trans="identity",name="True length")+
  scale_y_continuous(trans="identity",name="True - Estimated length (bias)")+
  stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.5))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  geom_hline(color="blue",linetype=3,yintercept = 0)+
  theme_bw()+
  theme(legend.position = c(.23,.975), legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))
ggsave("S100-bias.pdf",width=5,height = 9.2)
