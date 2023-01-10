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
