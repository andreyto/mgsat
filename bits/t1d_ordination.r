

taxa.level = "6"
m_a.norm.all = norm.count.m_a(res$res.tasks[[1]]$res.taxa.levels[[taxa.level]]$res.tests$m_a,method="norm.clr")
ranked.feats = get.feature.ranking(res$res.tasks[[1]]$res.taxa.levels[[taxa.level]]$res.tests,only.names = F)
de.res = as.data.frame(res$res.tasks[[1]]$res.taxa.levels[[taxa.level]]$res.tests$deseq2$results[[1]])
de.res = de.res[rownames(de.res) %in% colnames(m_a.norm.all$count),]
de.res$Rank = seq(nrow(de.res))
de.res$Up = ifelse(de.res$log2FoldChange>0,"T1D","Control")
de.res$Increase = ifelse(quant.mask(1/de.res$Rank,0.8),de.res$Up,"Neutral")

de.res = de.res[de.res$Increase != "Neutral",]

m_a.sel = subset.m_a(m_a.norm.all,select.count = rownames(de.res))

if(F) {
ranked.feats = get.feature.ranking(res$res.tasks[[1]]$res.taxa.levels[[taxa.level]]$res.tests)$ranked
m_a.sel = subset.m_a(m_a.norm.all,select.count = ranked.feats[1:10],
                 subset=(! (rownames(m_a.norm.all$count) %in% c("ASOJON","AHYTAH")) ))
GP1 = m_a.to.phyloseq(m_a.sel)
#GP1 = m_a.to.phyloseq(subset.m_a(norm.count.m_a(m_a,method="norm.clr",method.args=list(offset=0.0001)),subset = (! (rownames(m_a.norm$count) %in% c("ASOJON","AHYTAH")) )))
print(plot_ordination(GP1, ordinate(GP1,method="RDA",distance="euclidean"), 
                      type = "biplot", color = "T1D", title = "taxa",label="Feature",axes=c(1,2)) + 
        geom_text(aes(label=Feature),size=rel(6)))
}
#attr.taxa=data.frame(ranked=quant.mask(ranked.feats$ranked,0.8))
attr.taxa=de.res
ph = m_a.to.phyloseq(m_a.sel,attr.taxa=attr.taxa)
#45
#point_label="T1D"
#"kamada.kawai"
ggsave("net.pdf",
       phyloseq:::plot_net(ph,distance = "euclidean",
                           maxdist = 60, 
                           laymeth="kamada.kawai",
                           hjust=-0.2,
                           type="taxa", color="Increase",
                           point_label="Feature")

#+theme(legend.position = "none")
)
