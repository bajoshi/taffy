pro readifu,inputfile,binput,rinput,nx1,nx2,ny1,ny2

; program to read and plot individual spaxels from a lzifuout file

print,'ploting ',nx1,nx2,ny1,ny2


;usage give variable xy range to plot
; readifu,'HCG95C_1_comp.fits','../data/HCG95C_B.fits','../data/HCG95C_R.fits',nx1,nx2,ny1,ny2


cr=''
bin=mrdfits(binput,0,hb0)
rin=mrdfits(rinput,0,hr0)

B_CONTINUUM=mrdfits(inputfile,1,h1)
R_CONTINUUM=mrdfits(inputfile,2,h2)
B_ADDPOLY=mrdfits(inputfile,3,h3)
R_ADDPOLY=mrdfits(inputfile,4,h4)
B_MPOLY=mrdfits(inputfile,5,h5)
R_MPOLY=mrdfits(inputfile,6,h6)
B_RESIDFIT=mrdfits(inputfile,7,h7)
R_RESIDFIT=mrdfits(inputfile,8,h8)
B_CONT_MASK=mrdfits(inputfile,9,h9)
R_CONT_MASK=mrdfits(inputfile,10,h10)
B_LINE=mrdfits(inputfile,11,h11)
R_LINE=mrdfits(inputfile,12,h12)
;atv,R_LINE
;B_LINE_COMP1=mrdfits(inputfile,13,h13)
;R_LINE_COMP1=mrdfits(inputfile,14,h14)
;V=mrdfits(inputfile,15,h15)
;V_ERR=mrdfits(inputfile,16,h16)
;VDISP=mrdfits(inputfile,17,h17)
;VDISP_ERR=mrdfits(inputfile,18,h18)

;CHI2=mrdfits(inputfile,19,h19)
;DOF=mrdfits(inputfile,20,h20)
;HBETA=mrdfits(inputfile,21,h21)
;HBETA_ERR=mrdfits(inputfile,22,h22)
;OIII5007=mrdfits(inputfile,23,h23)
;OIII5007_ERR=mrdfits(inputfile,24,h24)
;OI6300=mrdfits(inputfile,25,h25)
;OI6300_ERR=mrdfits(inputfile,26,h26)
;OI6364=mrdfits(inputfile,27,h27)
;OI6364_ERR=mrdfits(inputfile,28,h28)
;HALPHA=mrdfits(inputfile,29,h29)
;HALPHA_ERR=mrdfits(inputfile,30,h30)
;NII6583=mrdfits(inputfile,31,h31)
;NII6583_ERR=mrdfits(inputfile,32,h32)
;SET=mrdfits(inputfile,33,h33)
cr=''



print,'Blue spectrum'

for i=nx1,nx2 do begin
for j=ny1,ny2 do begin
;set range of plot
nun=where(finite(B_continuum[i,j,*]))

minp=min(B_continuum[i,j,nun])
minp=minp-minp*2.
maxp=max(B_continuum[i,j,nun])
maxp=maxp+maxp*2.
plot,B_continuum[i,j,*],yrange=[-15,100];,yrange=[minp,2.*maxp]
oplot,b_line[i,j,*],color=2
oplot,B_cont_mask[i,j,*]*maxp,color=4
oplot,bin[i,j,*],color=2
read,cr,prompt='Hit return for model'
plot,bin[i,j,*],yrange=[-15,100];,yrange=[minp,2.*maxp]
oplot,B_cont_mask[i,j,*]*maxp,color=4
oplot,b_line[i,j,*],color=2
oplot,B_continuum[i,j,*]+b_line[i,j,*],color=2
print,i,j
read,cr,prompt='next plot cr '
endfor
endfor

print,'now plot mpoly,resfit '

;for i=nx1,nx2 do begin
;for j=ny1,ny2 do begin
;set range of plot
;nun=where(finite(B_mpoly[i,j,*]))
;minp=min(B_continuum[i,j,nun])
;minp=minp-minp*0.3
;maxp=max(B_continuum[i,j,nun])
;maxp=maxp+maxp*0.3
;plot,B_continuum[i,j,*]
;oplot,B_mpoly[i,j,*],color=4
;plot,B_residfit[i,j,*],color=2
;print,i,j
;read,cr,prompt='next plot cr '
;endfor
;endfor



print,'Blue zoom'


for i=nx1,nx2 do begin
for j=ny1,ny2 do begin
;set range of plot
minp=min(B_continuum[i,j,700:1300])
minp=minp-minp*2.
maxp=max(B_continuum[i,j,700:1300])
maxp=maxp+maxp*2.0
plot,B_continuum[i,j,*],xrange=[500,1300],yrange=[-15,100];,yrange=[minp,2.*maxp]
oplot,B_cont_mask[i,j,*]*maxp,color=4
oplot,bin[i,j,*],color=2
read,cr,prompt='Hit return for model'
plot,bin[i,j,*],xrange=[500,1300],yrange=[-15,100];,yrange=[minp,2.*maxp]
oplot,B_continuum[i,j,*]+b_line[i,j,*],color=2

print,i,j
read,cr,prompt='next plot cr '
endfor
endfor
read,cr,prompt='Now for Red '


for i=nx1,nx2 do begin
for j=ny1,ny2 do begin
;set range of plot
nun=where(finite(R_continuum[i,j,*]))
minp=min(R_continuum[i,j,nun])
minp=minp-minp*0.5
maxp=max(R_continuum[i,j,nun])
maxp=maxp+maxp*0.5
plot,R_continuum[i,j,*],yrange=[-15,100];,yrange=[minp,2.*maxp]
oplot,R_cont_mask[i,j,*]*maxp,color=4
oplot,rin[i,j,*],color=2
read,cr,prompt='Hit return for model'
plot,rin[i,j,*],yrange=[-15,100];,yrange=[minp,2.*maxp]
oplot,R_continuum[i,j,*]+r_line[i,j,*],color=2
print,i,j
read,cr,prompt='next plot cr '
endfor
endfor


;print,'Now plot polynomials and residuals'

;for i=nx1,nx2 do begin
;for j=ny1,ny2 do begin
;set range of plot
;nun=where(finite(B_mploy[i,j,*]))

;minp=min(B_continuum[i,j,nun])
;minp=minp-minp*0.3
;maxp=max(B_continuum[i,j,nun])
;maxp=maxp+maxp*0.3
;plot,R_continuum[i,j,*]
;plot,R_mpoly[i,j,*],color=4
;plot,R_residfit[i,j,*],color=2
;print,i,j
;read,cr,prompt='next plot cr '
;endfor
;endfor




read,cr,prompt='Now for Red H-alpha region'




for i=nx1,nx2 do begin
for j=ny1,ny2 do begin
;set range of plot

nun=where(finite(R_continuum[i,j,*]))


minp=min(R_continuum[i,j,nun])
minp=minp-minp*2.0
maxp=max(R_continuum[i,j,nun])
maxp=maxp+maxp*2.0
plot,R_continuum[i,j,*],xrange=[1000,2250],yrange=[-15,100];,yrange=[minp,2.0*maxp]
oplot,R_cont_mask[i,j,*]*maxp,color=4
oplot,rin[i,j,*],color=2

read,cr,prompt='Hit return for model'
plot,rin[i,j,*],xrange=[1000,2250],yrange=[-15,100]
oplot,R_continuum[i,j,*]+r_line[i,j,*],color=2



print,i,j
read,cr,prompt='next plot cr '
endfor
endfor
read,cr,prompt='Now for Red [OI] region'



for i=nx1,nx2 do begin
for j=ny1,ny2 do begin
;set range of plot

nun=where(finite(R_continuum[i,j,*]))


minp=min(R_continuum[i,j,nun])
minp=minp-minp*2.0
maxp=max(R_continuum[i,j,nun])
maxp=maxp+maxp*2.0
plot,R_continuum[i,j,*],xrange=[1150,1800],yrange=[-15,100];,yrange=[minp,2.0*maxp]
oplot,R_cont_mask[i,j,*]*maxp,color=4
oplot,rin[i,j,*],color=2

read,cr,prompt='Hit return for model'
plot,rin[i,j,*],xrange=[1150,1800],yrange=[-15,100]
oplot,R_continuum[i,j,*]+r_line[i,j,*],color=2



print,i,j
read,cr,prompt='next plot cr '
endfor
endfor


end
