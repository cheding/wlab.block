#' A Function to get fitness heatmap
#'
#' @description
#' This function allows you to get fitness heatmap.
#'
#' @param input single mutant fitness data table
#' @param wt_aa wt amino acid sequence
#' @param title heatmap title
#' @param legend_limits range of legend
#'
#' @return fitness heatmap
#' @export
#' @import data.table
fitness_heatmap<-function(
    input,
    wt_aa,
    title="fitness",
    legend_limits=NULL
){
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", "")))
  num<-nchar(wt_aa)+1
  input_single<-input
  input_single[,position:=AA_Pos1]
  input_single[,WT_AA:=wtcodon1]
  #create assistant data table
  # * represent STOP codon
  heatmap_tool_fitness<-data.table(wtcodon1 = rep(unlist(strsplit(wt_aa,"")),each=21),
                                   position = rep(2:num,each=21),
                                   codon1 = c(unlist(aa_list),"*"))
  heatmap_tool_fitness_anno_single<-merge(input_single,heatmap_tool_fitness,by=c("wtcodon1","position","codon1"),all=T)
  heatmap_tool_fitness_anno_single<-within(heatmap_tool_fitness_anno_single,
                                           codon1 <- factor(codon1,
                                                            levels = c("*","D","E","R","H","K","S","T","N","Q","C","G","P","A","V","I","L","M","F","W","Y")))
  heatmap_tool_fitness_anno_single[wtcodon1==codon1,nor_fitness_nooverlap:=0]
  ggplot2::ggplot()+
    ggplot2::theme_classic()+
    ggplot2::geom_tile(data=heatmap_tool_fitness_anno_single[position>1,],ggplot2::aes(x=position,y=codon1,fill=nor_fitness_nooverlap))+
    ggplot2::scale_x_discrete(limits=c(2:num),labels=c(2:num))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, vjust = 0.5,hjust = 0.5,
                                                       color = c(NA,NA,NA,rep(c("black",NA,NA,NA,NA),nchar(wt_aa)%/%5))))+
    ggplot2::scale_fill_gradient2(limits=legend_limits,low="#F4270C",mid="gray",high="#1B38A6",na.value = "white")+
    ggplot2::ylab("Mutant aa")+
    ggplot2::ggtitle(title)+
    ggplot2::labs(fill=NULL)+
    ggplot2::geom_text(data=heatmap_tool_fitness_anno_single[position>1&wtcodon1==codon1,],ggplot2::aes(x=position,y=codon1),label="-",size=3)+
    ggplot2::theme(text = ggplot2::element_text(size=5),
                   axis.ticks.x=ggplot2::element_blank(),
                   axis.ticks.y=ggplot2::element_blank(),
                   legend.position="bottom",
                   title = ggplot2::element_text(size=9,face = "bold"),
                   legend.text = ggplot2::element_text(size=6),
                   axis.title.x = ggplot2::element_text(size =9,face = "plain"),
                   axis.title.y = ggplot2::element_text(size =9,face = "plain"),
                   axis.text.x = ggplot2::element_text(size =8, angle=90, vjust=.5, hjust=1),
                   axis.text.y = ggplot2::element_text(family="Courier",angle=90,size=7.5, vjust = .5,hjust = .5,margin=ggplot2::margin(0,-0.5,0,0,"mm")),
                   legend.key.height= ggplot2::unit(3.1, 'mm'),
                   legend.key.width = ggplot2::unit(3.1, 'mm'),
                   legend.key.size = ggplot2::unit(1,"mm"),
                   plot.margin=ggplot2::margin(0,-0,0,0))+
    ggplot2::coord_fixed()
}
