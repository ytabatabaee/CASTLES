require(ggplot2);require(reshape2);require(scales)

s = read.csv('s100_all_branches.csv')
head(s)
unique(s$Method)
s$se = (s$l.est - s$l.true)^2 

im =  c("CASTLES (Integrals)", "CASTLES (Taylor 33, Coal)" ,
        "CASTLES (Taylor 33)" , "CASTLES (Taylor 34, Coal)" ,
        "CASTLES (Taylor 35, Coal)") 
variants = s$Method %in% im
s[s$Method=="CASTLES (Taylor 33, Handle OG)","Method"]="CASTLES"
names(s)=c("X"             ,      "Condition"         ,  "Method"        ,      "replicate"    ,       "Branch.Type"  ,
           "l.true"      ,        "l.est"        ,       "log10.l.true."   ,   
           "log10.l.est."    ,    "abserr"   ,"log10err" ,"se")
s$Method=factor(s$Method,levels=c("CASTLES"        ,           "CASTLES (Integrals)"  ,     "CASTLES (Taylor 33, Coal)" ,"CASTLES (Taylor 33)"  ,    
                         "CASTLES (Taylor 34, Coal)" ,"CASTLES (Taylor 35, Coal)" ,           "ERaBLE"     ,              
                        "Naive"                  ,   "Patristic(AVG)+FastME" ,    "Patristic(MIN)+FastME" ,
                        "Concat+RAxML"   ))
s$Condition =  factor(s$Condition) 
levels(s$Condition) = list("200bp" = "fasttree_genetrees_200_non", 
                           "400bp" = "fasttree_genetrees_400_non", 
                           "800bp" = "fasttree_genetrees_800_non",
                           "1600bp" = "fasttree_genetrees_1600_non",
                           "true gene trees" = "truegenetrees")

q= read.csv('quartets_all_branches.csv')
head(q)
unique(q$Method)
q$se = (q$l.est - q$l.true)^2 
q[q$Method=="CASTLES (Taylor 34, Coal)","Method"]="CASTLES"
qvariants = q$Method %in% im
names(q) =  c(names(s)[1:4],"QuartetType", names(s)[5:12])
q$Condition = factor(q$Condition)
levels(q$Condition) = list("Homogeneous" = "no_variation",
                                 "Sp" = "only_hs",
                                 "Loc" =  "only_hl",
                                 "Sp, Loc" = "hs_hl",
                                 "Sp, Loc, Sp/Loc"  =  "hs_hl_hg" ,
                                 "Sp, Loc, Sp/Loc, highILS"  = "hs_hl_hg_highr")


ggplot(aes(x=l.true,y=l.est,color=Branch.Type,linetype=Branch.Type),
       data=s[!variants,])+
  facet_grid(Condition~Method)+
  geom_point(alpha=0.05)+
  scale_x_continuous(trans="log10",name="True length")+
  scale_y_continuous(trans="log10",name="Estimated length")+
  stat_smooth(color="grey30",se=F,method="glm",formula=y ~ poly(x, 2))+
  scale_color_brewer(palette = "Set2")+
  coord_cartesian(xlim=c(10^-4.4,1),ylim=c(10^-4.4,1))+
  geom_abline(color="blue",linetype=3)+
  theme_bw()+
  theme(legend.position = c(.94,.1)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
ggsave("S100-correlation.png",width=12.5,height = 9)

ggplot(aes(x=l.true,y=l.est,color=Method,linetype),
       data=q[!qvariants,])+
  facet_grid(Branch.Type~Condition)+
  scale_x_continuous(trans="log10",name="True length")+
  scale_y_continuous(trans="log10",name="Estimated length")+
  stat_smooth(se=F,alpha=1,size=0.4,method="glm",formula=y ~ poly(x, 2))+
  scale_color_brewer(palette = "Spectral")+
  coord_cartesian(xlim=c(10^-3.8,0.9),ylim=c(10^-3.8,0.9))+
  geom_abline(color="grey30",linetype=2)+
  geom_point(alpha=0.2,size=0.7)+
  theme_bw()+
  theme(legend.position = "bottom")
ggsave("quartet-correlation-spec.pdf",width=12.5,height = 4.8)


ggplot(aes(x=l.true,y=l.est,color=Method,linetype),
       data=s[!variants,])+
  facet_grid(Branch.Type~Condition)+
  scale_x_continuous(trans="log10",name="True length")+
  scale_y_continuous(trans="log10",name="Estimated length")+
  scale_color_brewer(palette = "Spectral")+
  coord_cartesian(xlim=c(10^-4,0.9),ylim=c(10^-4,0.9))+
  geom_abline(color="grey30",linetype=2)+
  geom_point(alpha=0.1,size=0.5)+
  stat_smooth(se=F,alpha=1,size=0.7,method="glm",formula=y ~ poly(x, 2))+
  theme_bw()+
  theme(legend.position = "bottom")
ggsave("S100-correlation-2.png",width=12,height = 7)


ggplot(aes(x=sub("+","\n",Method,fixed=T),
           y=l.true-l.est,color=Branch.Type),
       data=s[!variants,])+
  facet_wrap(~reorder(sub("_non","bp",sub("fasttree_genetrees_","",Condition)),l.true-l.est),ncol=1)+
  #scale_x_continuous(trans="identity",name="True length")+
  scale_y_continuous(trans="identity",name=expression("True" - "Estimated length (bias)"))+
  stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.5))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  geom_hline(color="blue",linetype=3,yintercept = 0)+
  theme_bw()+
  theme(legend.position = c(.23,.975), legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))
#ggsave("S100-bias.pdf",width=5,height = 9.2)

ggplot(aes(x=sub("+","\n",Method,fixed=T),
           y=l.est/l.true,color=Branch.Type),
       data=q[!qvariants,])+
  facet_wrap(~Condition,ncol=2)+
  #scale_x_continuous(trans="identity",name="True length")+
  scale_y_continuous(trans=scales::trans_new(name="mylog",transform =  function(x) sign(x)*log10(abs(x)),
                                             inverse = function(y) ifelse(y==0,1,sign(y)*10^(abs(y)))),
                     name=expression("Estimated" / "True length (bias)"),
                     breaks = c(1,-10,100,1000,10,-100,-1000))+
  stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.5))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  geom_hline(color="blue",linetype=3,yintercept = 1)+
  theme_bw()+
  theme(legend.position = c(.23,.975), legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))


ggplot(aes(x= Condition,
           y=l.true-l.est,color=Method),
       data=q[!qvariants,])+
  facet_wrap(~reorder(Branch.Type,-l.true),ncol=2)+
  scale_x_discrete(label=function(x) gsub(",","\n",x,fixed=T))+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  #scale_x_continuous(trans="identity",name="True length")+
  scale_y_continuous(name=expression("True" - "Estimated length (bias)"))+
  stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.5))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_bw()+
  theme(legend.position = c(.75,.18), legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  guides(color=guide_legend(nrow=3, byrow=TRUE))
ggsave("quartet-bias.pdf",width=8,height =  4.7)


ggplot(aes(x= Condition,
           y=l.true-l.est,color=Method),
       data=s[!variants,])+
  facet_wrap(~reorder(Branch.Type,-l.true),ncol=2)+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  #scale_x_continuous(trans="identity",name="True length")+
  scale_y_continuous(name=expression("True" - "Estimated length (bias)"))+
  stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.5))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_bw()+
  theme(legend.position = c(.73,.9), legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))
ggsave("S100-bias.pdf",width=11,height =  4.5)

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


ggplot(aes(x=Condition,
           y=abserr,color=Method),
       data=dcast(data=s[!variants,],Condition+Method+replicate~'abserr' ,value.var = "abserr",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  coord_cartesian(ylim=c(0,0.045))
ggsave("S100-error-perrep.pdf",width=6.5,height = 5)

ggplot(aes(x=Condition,
           y=abserr,color=Method),
       data=dcast(data=q[!qvariants,],Condition+Method+replicate~'abserr' ,value.var = "abserr",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  scale_x_discrete(label=function(x) gsub(",","\n",x,fixed=T))+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_bw()+
  theme(legend.position = c(.34,.9), legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  coord_cartesian(ylim=c(0,0.75))+
  guides(color=guide_legend(nrow=2, byrow=TRUE))
ggsave("quartet-error-perrep.pdf",width=7,height = 4.5)

ggplot(aes(x=Condition,
       y=sqrt(se),color=Method),
      data=dcast(data=s[!variants,],Condition+Method+replicate~'se' ,value.var = "se",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Root mean square error")+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  coord_cartesian(ylim=c(0,0.08))
ggsave("S100-error-rmse.pdf",width=6.5,height = 5)

ggplot(aes(x=Condition,
           y=sqrt(se),color=Method),
       data=dcast(data=q[!qvariants,],Condition+Method+replicate~'se' ,value.var = "se",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Root mean square error")+
  scale_x_discrete(label=function(x) gsub(",","\n",x,fixed=T))+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  coord_cartesian(ylim=c(0,0.8))
ggsave("quartet-error-rmse.pdf",width=6.5,height = 5)

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
#ggsave("S100-error-rmse.pdf",width=5,height = 4.7)

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
#ggsave("S100-error-logabs.pdf",width=5,height = 4.7)


ggplot(aes(x=Condition,
           y=log10err,color=Method),
       data=dcast(data=s[!variants,],Condition+Method+replicate~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))))+
  scale_y_continuous(trans="identity",name="Mean log error")+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  coord_cartesian(ylim=c(0,1.5))
ggsave("S100-error-logabs.pdf",width=6.5,height = 5)

ggplot(aes(x=Condition,
           y=log10err,color=Method),
       data=dcast(data=q[!qvariants,],Condition+Method+replicate~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))))+
  scale_y_continuous(trans="identity",name="Mean log error")+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_bw()+
  scale_x_discrete(label=function(x) gsub(",","\n",x,fixed=T))+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  coord_cartesian(ylim=c(0,1.5))
ggsave("S100-error-logabs.pdf",width=6.5,height = 5)
