require(ggplot2);require(reshape2);require(scales)

s = read.csv('s100_all_branches.csv')
head(s)
unique(s$Method)
s$se = (s$l.est - s$l.true)^2 

im =  c("CASTLES (Integrals)", "CASTLES (Taylor 33, Coal)" ,
        "CASTLES (Taylor 33)" , "CASTLES (Taylor 34, Coal)" ,
        "CASTLES (Taylor 35, Coal)", "Naive") 
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
q[q$Method=="CASTLES (Taylor 33)","Method"]="CASTLES"
qvariants = q$Method %in% im
names(q) =  c(names(s)[1:4],"QuartetType", names(s)[5:12])
q$Condition = factor(q$Condition)
levels(q$Condition) = list("Homogeneous" = "no_variation",
                                 "Sp" = "only_hs",
                                 "Loc" =  "only_hl",
                                 "Sp, Loc" = "hs_hl",
                                 "Sp, Loc, Sp/Loc"  =  "hs_hl_hg" ,
                                 "Sp, Loc, Sp/Loc, highILS"  = "hs_hl_hg_highr")

m = read.csv('mvroot_all_branches_estgt.csv')
head(m)
nrow(m)
unique(m$Method)
m$se = (m$l.est - m$l.true)^2 
mvariants = m$Method %in% c("Naive")
names(m) =  c(names(s)[1:4],"AD", "GTEE", names(s)[5:12])
m$outgroup = factor(grepl("outgroup.0", m$Condition))
m$ratevar =  unique(sub(".genes.*","",sub("outgroup.*.species.","",m$Condition)))
m$Method = factor(m$Method, levels=c("CASTLES" ,"ERaBLE", "Naive" , 
                           "Patristic(AVG)+FastME" ,"Patristic(MIN)+FastME", "Concat+RAxML"))

### Comment out to include negative branch lengths.
summary(with(m[m$Method =="CASTLES" ,],l.est < 0))
m$l.est = ifelse(m$l.est <=0, 1e-6, m$l.est)
m$log10err = log10(m$l.est / m$l.true )
m$abserr = abs(m$l.true - m$l.est)
m$se = (m$l.est - m$l.true)^2 
q$l.est = ifelse(q$l.est <=0, 1e-6, q$l.est)
q$log10err = log10(q$l.est / q$l.true )
q$abserr = abs(q$l.true - q$l.est)
q$se = (q$l.est - q$l.true)^2 
s$l.est = ifelse(s$l.est <=0, 1e-6, s$l.est)
s$log10err = log10(s$l.est / s$l.true )
s$abserr = abs(s$l.true - s$l.est)
s$se = (s$l.est - s$l.true)^2 



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
  stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_bw()+
  theme(legend.position = c(.75,.14), legend.direction = "horizontal",
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
  stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))#+
  guides(color=guide_legend(nrow=3, byrow=TRUE))
ggsave("S100-bias.pdf",width=8,height =  4)

ggplot(aes(x= Condition,
           y=l.true-l.est,color=Method),
       data=s[s$Method %in% c("CASTLES","Naive"),])+
  facet_wrap(~reorder(Branch.Type,-l.true),ncol=2)+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  #scale_x_continuous(trans="identity",name="True length")+
  scale_y_continuous(name=expression("True" - "Estimated length (bias)"))+
  stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  stat_summary(aes(group=Method),fun.data = mean_sdl,position = position_dodge(width=0.75),geom="line")+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))#+
guides(color=guide_legend(nrow=3, byrow=TRUE))

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
  coord_cartesian(ylim=c(0,0.045),xlim=c(1,5),clip="off")+
  annotate(geom="text",label="b)", x = 0.17, y = .046, size = 6)
ggsave("S100-error-perrep.pdf",width=7,height = 5)

ggplot(aes(x=Condition,
           y=abserr,color=Method),
       data=dcast(data=q[!qvariants,],Condition+Method+replicate~'abserr' ,value.var = "abserr",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  scale_x_discrete(label=function(x) gsub(",","\n",x,fixed=T))+
  geom_boxplot(outlier.alpha = 0.3,width=0.86,outlier.size=1)+
  stat_summary(position = position_dodge(width=0.86))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_bw()+
  theme(legend.position = c(.3,.9), legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  coord_cartesian(ylim=c(0,0.68),xlim=c(1,6),clip="off")+
  annotate(geom="text",label="b)", x = 0.185, y = .69, size = 5)+
  guides(color=guide_legend(nrow=2, byrow=TRUE))
ggsave("quartet-error-perrep.pdf",width=7,height = 4.5)


ggplot(aes(x=ratevar, y=abserr,color=Method,fill=outgroup),
       data=dcast(data=m[!mvariants,],
                  outgroup+ratevar+Method+replicate~'abserr' ,value.var = "abserr",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  scale_x_discrete(labels=c("High","Med","Low"),name="Clock deviation")+
  geom_boxplot(outlier.alpha = 0.3,width=0.9,outlier.size = 1)+
  stat_summary(position = position_dodge(width=0.9))+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_manual(values=c("white","grey70"),name="",labels=c("With outgroup","No outgroup"))+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_bw()+
  theme(legend.position =  "bottom", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0),
        legend.box.margin = margin(0), legend.margin = margin(0)
        )+
  coord_cartesian(ylim=c(0,0.1) ,xlim=c(1,3), clip = "off")+
  guides(color=guide_legend(nrow=2,  byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))+
  annotate(geom="text",label="a)", x = 0.13, y = 0.104, size = 5) 
ggsave("MV-error-perrep.pdf",width=6.5,height = 4.3)

ggplot(aes(x=Method, y=abserr,fill=ratevar,color=outgroup,shape=outgroup),
       data=dcast(data=m[!mvariants,],
                  outgroup+ratevar+Method+replicate~'abserr' ,value.var = "abserr",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T))+
  geom_boxplot(outlier.alpha = 0.3,width=0.9,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.9))+
  #geom_boxplot(outlier.size = 0)+
  scale_color_manual(values=c("black","grey50"),name="",labels=c("With outgroup","No outgroup"))+
  scale_shape(name="",labels=c("With outgroup","No outgroup"))+
  scale_fill_brewer(palette = 1,labels=c("High","Med","Low"),name="Clock deviation",direction = -1)+
  theme_bw()+
  theme(legend.position =  "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  coord_cartesian(ylim=c(0,0.08))+
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         fill=guide_legend(nrow=1, byrow=TRUE))
ggsave("MV-error-perrep-bymethod.pdf",width=6.2*0.95,height = 4.3*0.95)

ggplot(aes(color=Method, y=log10err,x=cut(AD,4)),
       data=merge(
         dcast(data=m[!mvariants & m$outgroup ==FALSE,],
                  outgroup+ratevar+Method+replicate~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))),
         dcast(data=m[!mvariants  & m$outgroup ==FALSE,], outgroup+replicate+ratevar~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean log error")+
  #facet_wrap(~outgroup,ncol=2,labeller = label_both)+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T),name="True gene tree discordance (ILS)")+
  geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.8))+
  #geom_boxplot(outlier.size = 0)+
  scale_color_manual(values=c("black","grey50"),name="",labels=c("With outgroup","No outgroup"))+
  scale_shape(name="",labels=c("With outgroup","No outgroup"))+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_bw()+
  theme(legend.position =  "bottom", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0,size=11))+
  coord_cartesian()+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-logerr-perrep-ILS-bymethod.pdf",width=6.2*0.95,height = 4.3*0.95)

ggplot(aes(color=Method, y=log10err,x=cut(AD,4)),
       data=merge(
         dcast(data=m[!mvariants & m$outgroup ==FALSE,],
               outgroup+ratevar+Method+replicate~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))),
         dcast(data=m[!mvariants  & m$outgroup ==FALSE,], outgroup+replicate+ratevar~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean log10 error")+
  #facet_wrap(~outgroup,ncol=2,labeller = label_both)+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T),name="True gene tree discordance (ILS)")+
  geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.8))+
  #geom_boxplot(outlier.size = 0)+
  scale_color_manual(values=c("black","grey50"),name="",labels=c("With outgroup","No outgroup"))+
  scale_shape(name="",labels=c("With outgroup","No outgroup"))+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_bw()+
  theme(legend.position =  "none", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0,size=11))+
  coord_cartesian(ylim=c(0.1,1),xlim=c(1,3), clip = "off")+
  annotate(geom="text",label="b)", x = 0.13, y = 1.04, size = 5) + 
  geom_text(aes(color="Patristic(MIN)+FastME",y=0.6,label="Log\nError\n>2"),position = position_nudge(x  = 0.165),size=2.5)
ggsave("MV-logerr-perrep-ILS-bymethod-nopatristic.pdf",width=6.2,height = 4)


ggplot(aes(color=Method, y=log10err,x=AD),
       data=merge(
         dcast(data=m[!mvariants & m$outgroup ==FALSE,],
               outgroup+ratevar+Method+replicate~'log10err' ,value.var = "log10err",fun.aggregate = function(x) mean(abs(x))),
         dcast(data=m[!mvariants  & m$outgroup ==FALSE,], outgroup+replicate+ratevar~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean log10 error",limits = c(0,1))+
  #facet_wrap(~outgroup,ncol=2,labeller = label_both)+
  #scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T),name="True gene tree discordance (ILS)")+
  geom_point(alpha=0.5)+
  stat_smooth(se=F,method="lm")+
  #geom_boxplot(outlier.size = 0)+
  scale_color_manual(values=c("black","grey50"),name="",labels=c("With outgroup","No outgroup"))+
  scale_shape(name="",labels=c("With outgroup","No outgroup"))+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_bw()+
  theme(legend.position =  "right", legend.direction = "vertical",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0,size=11))+
  stat_cor(method = "pearson", 
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  coord_cartesian()
ggsave("MV-logerr-perrep-ILS-correlation.pdf",width=6.5,height = 5)

ggplot(aes(color=Method, y=abserr,x=cut(AD,4)),
       data=merge(
         dcast(data=m[!mvariants & m$outgroup ==FALSE,],
               outgroup+ratevar+Method+replicate~'abserr' ,value.var = "abserr",fun.aggregate = function(x) mean(abs(x))),
         dcast(data=m[!mvariants  & m$outgroup ==FALSE,], outgroup+replicate+ratevar~'AD' ,value.var = "AD",fun.aggregate = mean)))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  #facet_wrap(~outgroup,ncol=2,labeller = label_both)+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T),name="True gene tree discordance (ILS)")+
  geom_boxplot(outlier.alpha = 0.3,width=0.8,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.8))+
  #geom_boxplot(outlier.size = 0)+
  scale_color_manual(values=c("black","grey50"),name="",labels=c("With outgroup","No outgroup"))+
  scale_shape(name="",labels=c("With outgroup","No outgroup"))+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_bw()+
  theme(legend.position =  "bottom", legend.direction = "horizontal",
        legend.box.margin = margin(0), legend.margin = margin(0),
        axis.text.x = element_text(angle=0,size=11))+
  coord_cartesian()+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE))
ggsave("MV-abserror-perrep-ILS-bymethod.pdf",width=6.2*0.95,height = 4.3*0.95)

ggplot(aes(x=Method, y=sqrt(se),fill=ratevar,color=outgroup,shape=outgroup),
       data=dcast(data=m[!mvariants,],
                  outgroup+ratevar+Method+replicate~'se' ,value.var = "se",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Root mean square error")+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T))+
  geom_boxplot(outlier.alpha = 0.3,width=0.9,outlier.size = 0.8)+
  stat_summary(position = position_dodge(width=0.9))+
  #geom_boxplot(outlier.size = 0)+
  scale_color_manual(values=c("black","grey50"),name="",labels=c("With outgroup","No outgroup"))+
  scale_shape(name="",labels=c("With outgroup","No outgroup"))+
  scale_fill_brewer(palette = 1,labels=c("High","Med","Low"),name="Clock deviation",direction = -1)+
  theme_bw()+
  theme(legend.position =  "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  coord_cartesian(ylim=c(0,0.1))+
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         fill=guide_legend(nrow=1, byrow=TRUE))
ggsave("MV-rmse-bymethod.pdf",width=6.4,height = 5)

ggplot(aes(x=Method, y=l.true-l.est,color=ratevar,shape=outgroup),
       data=m[!mvariants,])+
  scale_y_continuous(trans="identity",name=expression("True" - "Estimated length (bias)"))+
  scale_x_discrete(label=function(x) gsub("+","\n",x,fixed=T))+
  stat_summary(position = position_dodge(width=0.9),size=0.8,fun.data = mean_sdl)+
  #geom_boxplot(outlier.size = 0)+
  #scale_color_manual(values=c("black","grey50"),name="",labels=c("With outgroup","No outgroup"))+
  scale_shape(name="",labels=c("With outgroup","No outgroup"))+
  scale_color_brewer(palette = 1,labels=c("High","Med","Low"),name="Clock deviation",direction = -1)+
  theme_bw()+
  theme(legend.position =  "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0))+
  geom_hline(color="grey50",linetype=1,yintercept = 0)+
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         fill=guide_legend(nrow=1, byrow=TRUE))
ggsave("MV-bias-bymethod.pdf",width=6.4,height = 5)



ggplot(aes(x=ratevar, y=abserr,color=Method,linetype=outgroup,shape=outgroup),
       data=dcast(data=m[!mvariants,],
                  outgroup+ratevar+Method+replicate~'abserr' ,value.var = "abserr",fun.aggregate = mean))+
  scale_y_continuous(trans="identity",name="Mean absolute error")+
  scale_x_discrete(labels=c("High","Med","Low"),name="Clock deviation")+
  stat_summary(aes(group=interaction(outgroup,Method)),geom="line")+
  stat_summary()+
  #geom_boxplot(outlier.size = 0)+
  scale_fill_manual(values=c("white","grey70"),name="",labels=c("No outgroup","With outgroup"))+
  scale_color_brewer(palette = "Dark2",name="")+
  theme_bw()+
  theme(legend.position =  "bottom", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0))+
  coord_cartesian()+
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         linetype=guide_legend(nrow=2, byrow=TRUE))

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
ggsave("S100-error-rmse.pdf",width=7,height = 5)

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
  coord_cartesian(ylim=c(0,1))+
  geom_text(aes(color="Patristic(MIN)+FastME",y=0.6,label="Log\nError\n>1.5"),position = position_nudge(x  = 0.165),size=2.5)
ggsave("S100-error-logabs.pdf",width=7,height = 5)

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
  coord_cartesian(ylim=c(0,4.5))
ggsave("quartet-error-logabs.pdf",width=6.5,height = 5)



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




t = read.csv('s100_time_memory.csv')
head(t)
ggplot(aes(x=reorder(Condition,time_s),y=time_s/60,color=Method,group=Method),
       data=t[t$Method %in% c("CASTLES+ASTER","Concat+RAxML","Patristic(ALL)+ERaBLE" ,   "Patristic(AVG)+FastME(AVG)" ,    "Patristic(MIN)+FastME(MIN)" ),])+
  stat_summary(geom="line")+
  stat_summary()+
  scale_fill_brewer(palette = "Dark2",name="")+
  scale_color_brewer(palette = "Dark2",name="")+
  scale_y_continuous(trans="log10",name="time (minutes)" )+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0),
        legend.box.margin = margin(0), legend.margin = margin(0)
        )+
  coord_cartesian(xlim=c(1,5),clip="off")+
  scale_x_discrete(label=function(x) gsub(" ","\n",x,fixed=T))+
  annotate(geom="text",label="a)", x = 0.13, y = 20.04, size = 5) +
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         linetype=guide_legend(nrow=2, byrow=TRUE))
ggsave("S100-time.pdf",width=6.3,height =5)

ggplot(aes(x=reorder(Condition,time_s),y=mem_gb,color=Method,group=Method),
       data=t[t$Method %in% c("CASTLES+ASTER","Concat+RAxML","Patristic(ALL)+ERaBLE" ,   "Patristic(AVG)+FastME(AVG)" ,    "Patristic(MIN)+FastME(MIN)" ),])+
  stat_summary(geom="line")+
  stat_summary()+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  scale_y_continuous(trans="log10",name="Memory (Gb)" )+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0),
        legend.box.margin = margin(0), legend.margin = margin(0))+
  coord_cartesian(xlim=c(1,5),clip="off") +
  scale_x_discrete(label=function(x) gsub(" ","\n",x,fixed=T))+
  annotate(geom="text",label="b)", x = 0.03, y = 10.04, size = 5) +
  guides(color=guide_legend(nrow=2, byrow=TRUE),
         linetype=guide_legend(nrow=2, byrow=TRUE))
ggsave("S100-memory.pdf",width=6.3,height =5)
