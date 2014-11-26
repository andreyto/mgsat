require(vegan)
## this shows that count-based richness estimates get higher with
## higher depth of sequencing, even though they are designed
## to be independent of that
View(rowMeans(estimateR(m_a$count)-estimateR(rrarefy(m_a$count,400))))
## variability between rarifications to the same count
View(rowMeans(abs(estimateR(rrarefy(m_a$count,400))-estimateR(rrarefy(m_a$count,400)))))
View(rowMeans(abs(estimateR(rrarefy(m_a$count,400)))))


