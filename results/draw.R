require(ggplot2);require(reshape2);require(scales)

s = read.csv('s100_all_branches.csv')
head(s)
unique(s$Method)
s$se = (s$l.est - s$l.true)^2 
  
variants = s$Method %in% c("CASTLES (Integrals)", "CASTLES (Taylor 33, Coal)" ,
                    "CASTLES (Taylor 33)" , "CASTLES (Taylor 34, Coal)" ,
                    "CASTLES (Taylor 35, Coal)") 
s[s$Method=="CASTLES (Taylor 33, Handle OG)","Method"]="CASTLES"

names(s)=c("X"             ,      "Condition"         ,  "Method"        ,      "replicate"    ,       "Branch.Type"  ,
           "l.true"      ,        "l.est"        ,       "log10.l.true."   ,   
           "log10.l.est."    ,    "abserr"   ,"log10err" ,"se")

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
  facet_wrap(~reorder(sub("_non","bp",sub("fasttree_genetrees_","",Condition)),l.true-l.est),ncol=1)+
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

ggplot(aes(x=reorder(sub("_non","bp",sub("fasttree_genetrees_","",Condition)),l.true-l.est),
           y=abs(l.true-l.est),color=sub("+","\n",Method,fixed=T)),
       data=s[!variants,])+
  #facet_wrap(~reorder(sub("_non","bp",sub("fasttree_genetrees_","",Condition)),l.true-l.est),ncol=1)+
  #scale_x_continuous(trans="identity",name="True length")+
  scale_y_continuous(trans="identity",name="True - Estimated length (bias)")+
  geom_boxplot()+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  geom_hline(color="blue",linetype=3,yintercept = 0)+
  theme_bw()+
  theme(legend.position = c(.23,.975), legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  coord_cartesian(ylim=c(0,0.3))
ggsave("S100-error-all.pdf",width=5,height = 9.2)


ggplot(aes(x=reorder(sub("_non","bp",sub("fasttree_genetrees_","",Condition)),abserr),
           y=abserr,color=sub("+","\n",Method,fixed=T)),
       data=dcast(data=s[!variants,],Condition+Method+replicate~'abserr' ,value.var = "abserr",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_bw()+
  theme(legend.position = c(.38,.8), legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  coord_cartesian(ylim=c(0,0.076))
ggsave("S100-error-perrep.pdf",width=5,height = 4.7)


ggplot(aes(x=reorder(sub("_non","bp",sub("fasttree_genetrees_","",Condition)),se),
           y=sqrt(se),color=sub("+","\n",Method,fixed=T)),
       data=dcast(data=s[!variants,],Condition+Method+replicate~'se' ,value.var = "se",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Root mean square error")+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_bw()+
  theme(legend.position = c(.35,.92), legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  coord_cartesian(ylim=c(0,0.095))
ggsave("S100-error-rmse.pdf",width=5,height = 4.7)

ggplot(aes(x=reorder(sub("_non","bp",sub("fasttree_genetrees_","",Condition)),log10err),
           y=log10err,color=sub("+","\n",Method,fixed=T)),
       data=dcast(data=s[!variants,],Condition+Method+replicate~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))))+
  scale_y_continuous(trans="identity",name="Root mean square error")+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_bw()+
  theme(legend.position = c(.36,.92), legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  coord_cartesian(ylim=c(0,1.5))
ggsave("S100-error-logabs.pdf",width=5,height = 4.7)
