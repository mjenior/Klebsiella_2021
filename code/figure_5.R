
# Read in data
clinical <- '~/Desktop/repos/Klebsiella_2021/data/clinical_maxfit_reps/flux_samples.tsv'
clinical <- read.delim(clinical, sep='\t', header=TRUE, row.names=1)
clinical <- as.data.frame(apply(clinical, 2, as.numeric))
laboratory <- '~/Desktop/repos/Klebsiella_2021/data/laboratory_maxfit_reps/flux_samples.tsv'
laboratory <- read.delim(laboratory, sep='\t', header=TRUE, row.names=1)
laboratory <- as.data.frame(apply(laboratory, 2, as.numeric))
clinical_VALTA <- clinical$VALTA
laboratory_VALTA <- laboratory$VALTA
rm(clinical, laboratory)

clinical <- read.delim('/home/mjenior/Desktop/repos/Klebsiella_2021/data/clinical_dFBA.tsv', sep='\t', header=TRUE)
laboratory <- read.delim('/home/mjenior/Desktop/repos/Klebsiella_2021/data/laboratory_dFBA.tsv', sep='\t', header=TRUE)

library(vioplot)

pdf(file='~/Desktop/repos/Klebsiella_2021/results/Figure_5AB.pdf', width=7, height=5)
layout(matrix(c(1,2,
                3,4), nrow=2, ncol=2, byrow=TRUE))

par(mar=c(3,3,2,3), xpd=FALSE, las=1, mgp=c(1.5,0.75,0), lwd=2, xaxs='i')

plot(x=clinical$time, y=clinical$valine, type='l', col='white', xlim=c(0,5), ylim=c(0,100), 
     main='Clinical', xaxt='n', yaxt='n', xlab='', ylab='')
lines(x=clinical$time, y=clinical$valine, col='gray30', lwd=3)
lines(x=clinical$time+2.2, y=clinical$biomass, col='#B13AED', lwd=3)
segments(x0=0, y0=0.1, x1=2.2, col='#B13AED', lwd=3)
axis(1, at=seq(0,5,1), cex.axis=0.6, lwd=2) # hours
axis(2, at=seq(0,100,20), cex.axis=0.6, lwd=2) # percent valine
axis(4, at=seq(0,100,20), cex.axis=0.6, lwd=2) # biomass units...
box()
par(xpd=TRUE)
text(x=2.5, y=-30, font=2, labels='Hours', cex=0.8)
text(x=-0.7, y=50, font=2, labels='L-Valine', srt=90, col='gray30', cex=0.8)
text(x=5.7, y=50, font=2, labels='Biomass', srt=270, col='#B13AED', cex=0.8)
text(x=-0.75, y=115, font=2, labels='A')
par(xpd=FALSE)

plot(x=laboratory$time, y=laboratory$valine, type='l', col='white', xlim=c(0,3.5), ylim=c(0,100), 
     main='Laboratory', xaxt='n', yaxt='n', xlab='', ylab='')
lines(x=laboratory$time, y=laboratory$valine, col='gray30', lwd=3)
lines(x=laboratory$time, y=laboratory$biomass, col='#76EEC6', lwd=3)
axis(1, at=seq(0,3.5,0.5), cex.axis=0.6, lwd=2) # hours
axis(2, at=seq(0,100,20), cex.axis=0.6, lwd=2) # percent valine
axis(4, at=seq(0,100,20), cex.axis=0.6, lwd=2) # biomass units...
box()
par(xpd=TRUE)
text(x=1.75, y=-30, font=2, labels='Hours', cex=0.8)
text(x=-0.52, y=50, font=2, labels='L-Valine', srt=90, col='gray30', cex=0.8)
text(x=4, y=50, font=2, labels='Biomass', srt=270, col='#76EEC6', cex=0.8)
text(x=-0.525, y=115, font=2, labels='B')
par(xpd=FALSE)

plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10), axes=FALSE)
par(xpd=TRUE)
text(x=-1.5, y=11.5, font=2, labels='C')
par(xpd=FALSE)

par(mar=c(3,3,2,3), xpd=FALSE, las=1, mgp=c(1.5,0.7,0), lwd=2, xaxs='i')
vioplot(clinical_VALTA, 0, col=c('#B13AED', '#76EEC6'), xaxt='n', yaxt='n',
        ylim=c(0, 1200), ylab='Sampled Reaction Flux', lwd=1.7, drawRect=FALSE, yaxs='i')
axis(2, at=seq(0,1200,400), cex.axis=0.6, lwd=2)
text(x=2, y=100, 'inactive', cex=0.9)
legend('topright', legend='Valine transaminase activity', bty='n')
box()
par(xpd=TRUE)
text(x=c(1,2), y=-100, labels=c('Clinical','Laboratory'), cex=1.1)
text(x=0.2, y=1300, font=2, labels='D')
par(xpd=FALSE)

dev.off()






