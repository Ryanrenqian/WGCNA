# WGCNA
# This is my own Rscript for WGCNA in bioinformatic analysis
#Notice
when you run or refer to my codes, some of them maybe warn you.Hence, I write down some problem you might meet.
## Input
fpkm is recommend for your input. However, when readcounts is input this program. you may get warning, `int` type should be converted into float, you can divide them by float number,
or you can add them 0.0001 which would not change them significantly.
## export data into network plot software
if you choose separate data into several blocks or you plan to load data from block files, there is really a problem. The TOM matrix you obtain from those files is not complete a matrix but 
a condensed and flattened matrix.
you have to converted it to original type.
for example, if the origin is a M X N matrix, the condensed one is M X (N-1) /2. Based on this you can convert it back.
But when you have enough available memory, caculate TOM matrix alone and then you can get an complete matrix so that you needn\`t convert them at all.
Hope this can help you !!
# WGCNA.R
This is my script provided for you to do WGCNA, but there are still some bug when saving pictures. I don\`t know how to cope with it.
this script will provides you modules in you file, and will automatically choose 0.9 when choose softpower.
So, more work remains to be done.
> WGCNA.R -i your expression table -s your sample name
