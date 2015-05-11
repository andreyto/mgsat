library(ComplexHeatmap)
Heatmap(log(as.matrix(div.counts$e)),cluster_columns=F,show_row_names = FALSE, clustering_distance_rows = "pearson") + HeatmapAnnotation(df=df,which="row") + Heatmap(m_a.norm$count)
Heatmap(log(as.matrix(div.counts$e)),cluster_columns=F,show_row_names = FALSE, km=3, clustering_distance_rows = "pearson") + HeatmapAnnotation(df=df,which="row") + Heatmap(m_a.norm$count)
adonis(m_a.norm$count~SampleType,data=m_a.norm$attr,method="euclidean")
xtabs(~Year+SampleType,m_a.norm$attr)
with(m_a.norm$attr,{
Heatmap(log(as.matrix(div.counts$e)),name="Renyi Ind",
        cluster_columns=F,
        show_row_names = FALSE, 
        clustering_distance_rows = "pearson") + 
  Heatmap(WP.collection,name="WP.collection")+
  Heatmap(Hidradenitis,name="Hidradenitis")+
  Heatmap(Year,name="Year")+
  Heatmap(SampleType,name="SampleType")+
  Heatmap(m_a.norm$count,name="Norm. counts")
})