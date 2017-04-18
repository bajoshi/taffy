PRO plot_bpt_hc,inputfile,regionfile



; first divide up the image into several regions which will be
; displayed in different colors

;Three regions will be defined

nnx1=intarr(9)
nnx2=intarr(9)
nny1=intarr(9)
nny2=intarr(9)


readcol,regionfile,nnx1,nnx2,nny1,nny2
for i=0,8 do print,'Region ,i,  is ', nnx1[i],nnx2[i],nny1[i],nny2[i]






;usage    plot_bpt,'HCG95C_1_comp.fits'

	sig_cut = 3.0
	; Load data
	ha     = mrdfits(inputfile,29,h29)
	ha_err = mrdfits(inputfile,30,h30)
	hb     = mrdfits(inputfile,21,h21)
	hb_err = mrdfits(inputfile,22,h22)
	n2     = mrdfits(inputfile,31,h31)
	n2_err = mrdfits(inputfile,32,h32)
	o1     = mrdfits(inputfile,25,h25)
	o1_err = mrdfits(inputfile,26,h26)
	o3     = mrdfits(inputfile,23,h23)
	o3_err = mrdfits(inputfile,24,h24)
 cr=''
 read,cr,prompt='Do you want to create hardcopy ps file? y/n '

 if cr eq 'y' then begin
;Set up hardcopy plot                                                           
      psname=inputfile+'_bpt.eps'
      ifps,psname,7,6,/post
   endif else begin
  set_plot, 'X'
endelse







	; plotting
;	window,xsize=1000,ysize=300
;	!P.MULTI=[0,3,1]
	!P.CHARSIZE=2
	; N2-Halpha BPT
	plot,indgen(10),/nodata,xrange=[-1,1],yrange=[-1,0.8],/xsty,/ysty,ytitle='log([OIII]/Hbeta)',xtitle='log([NII]/Halpha)'
nx1=nnx1[0]
nx2=nnx2[0]
ny1=nny1[0]
ny2=nny2[0]

 	ind = where(n2[nx1:nx2,ny1:ny2,1] gt sig_cut * n2_err[nx1:nx2,ny1:ny2,1] and ha[nx1:nx2,ny1:ny2,1] gt sig_cut * ha_err[nx1:nx2,ny1:ny2,1] and $
        o3[nx1:nx2,ny1:ny2,1] gt sig_cut * o3_err[nx1:nx2,ny1:ny2,1] and hb[nx1:nx2,ny1:ny2,1] gt sig_cut * hb_err[nx1:nx2,ny1:ny2,1] )
	oplot,(alog10(n2[nx1:nx2,ny1:ny2,1]/ha[nx1:nx2,ny1:ny2,1]))[ind],(alog10(o3[nx1:nx2,ny1:ny2,1]/hb[nx1:nx2,ny1:ny2,1]))[ind],psym=4,color=cgcolor('blue'),thick=3
nx1=nnx1[1]
nx2=nnx2[1]
ny1=nny1[1]
ny2=nny2[1]

 	ind = where(n2[nx1:nx2,ny1:ny2,1] gt sig_cut * n2_err[nx1:nx2,ny1:ny2,1] and ha[nx1:nx2,ny1:ny2,1] gt sig_cut * ha_err[nx1:nx2,ny1:ny2,1] and $
        o3[nx1:nx2,ny1:ny2,1] gt sig_cut * o3_err[nx1:nx2,ny1:ny2,1] and hb[nx1:nx2,ny1:ny2,1] gt sig_cut * hb_err[nx1:nx2,ny1:ny2,1] )
	oplot,(alog10(n2[nx1:nx2,ny1:ny2,1]/ha[nx1:nx2,ny1:ny2,1]))[ind],(alog10(o3[nx1:nx2,ny1:ny2,1]/hb[nx1:nx2,ny1:ny2,1]))[ind],psym=4,color=cgcolor('blue'),thick=3
nx1=nnx1[2]
nx2=nnx2[2]
ny1=nny1[2]
ny2=nny2[2]

 	ind = where(n2[nx1:nx2,ny1:ny2,1] gt sig_cut * n2_err[nx1:nx2,ny1:ny2,1] and ha[nx1:nx2,ny1:ny2,1] gt sig_cut * ha_err[nx1:nx2,ny1:ny2,1] and $
        o3[nx1:nx2,ny1:ny2,1] gt sig_cut * o3_err[nx1:nx2,ny1:ny2,1] and hb[nx1:nx2,ny1:ny2,1] gt sig_cut * hb_err[nx1:nx2,ny1:ny2,1] )
	oplot,(alog10(n2[nx1:nx2,ny1:ny2,1]/ha[nx1:nx2,ny1:ny2,1]))[ind],(alog10(o3[nx1:nx2,ny1:ny2,1]/hb[nx1:nx2,ny1:ny2,1]))[ind],psym=4,color=cgcolor('blue'),thick=3
nx1=nnx1[3]
nx2=nnx2[3]
ny1=nny1[3]
ny2=nny2[3]

 	ind = where(n2[nx1:nx2,ny1:ny2,1] gt sig_cut * n2_err[nx1:nx2,ny1:ny2,1] and ha[nx1:nx2,ny1:ny2,1] gt sig_cut * ha_err[nx1:nx2,ny1:ny2,1] and $
        o3[nx1:nx2,ny1:ny2,1] gt sig_cut * o3_err[nx1:nx2,ny1:ny2,1] and hb[nx1:nx2,ny1:ny2,1] gt sig_cut * hb_err[nx1:nx2,ny1:ny2,1] )
	oplot,(alog10(n2[nx1:nx2,ny1:ny2,1]/ha[nx1:nx2,ny1:ny2,1]))[ind],(alog10(o3[nx1:nx2,ny1:ny2,1]/hb[nx1:nx2,ny1:ny2,1]))[ind],psym=5,color=cgcolor('black'),thick=3

nx1=nnx1[4]
nx2=nnx2[4]
ny1=nny1[4]
ny2=nny2[4]

 	ind = where(n2[nx1:nx2,ny1:ny2,1] gt sig_cut * n2_err[nx1:nx2,ny1:ny2,1] and ha[nx1:nx2,ny1:ny2,1] gt sig_cut * ha_err[nx1:nx2,ny1:ny2,1] and $
        o3[nx1:nx2,ny1:ny2,1] gt sig_cut * o3_err[nx1:nx2,ny1:ny2,1] and hb[nx1:nx2,ny1:ny2,1] gt sig_cut * hb_err[nx1:nx2,ny1:ny2,1] )
	oplot,(alog10(n2[nx1:nx2,ny1:ny2,1]/ha[nx1:nx2,ny1:ny2,1]))[ind],(alog10(o3[nx1:nx2,ny1:ny2,1]/hb[nx1:nx2,ny1:ny2,1]))[ind],psym=4,color=cgcolor('red'),thick=3


nx1=nnx1[5]
nx2=nnx2[5]
ny1=nny1[5]
ny2=nny2[5]

 	ind = where(n2[nx1:nx2,ny1:ny2,1] gt sig_cut * n2_err[nx1:nx2,ny1:ny2,1] and ha[nx1:nx2,ny1:ny2,1] gt sig_cut * ha_err[nx1:nx2,ny1:ny2,1] and $
        o3[nx1:nx2,ny1:ny2,1] gt sig_cut * o3_err[nx1:nx2,ny1:ny2,1] and hb[nx1:nx2,ny1:ny2,1] gt sig_cut * hb_err[nx1:nx2,ny1:ny2,1] )
	oplot,(alog10(n2[nx1:nx2,ny1:ny2,1]/ha[nx1:nx2,ny1:ny2,1]))[ind],(alog10(o3[nx1:nx2,ny1:ny2,1]/hb[nx1:nx2,ny1:ny2,1]))[ind],psym=4,color=cgcolor('red'),thick=3





nx1=nnx1[6]
nx2=nnx2[6]
ny1=nny1[6]
ny2=nny2[6]

 	ind = where(n2[nx1:nx2,ny1:ny2,1] gt sig_cut * n2_err[nx1:nx2,ny1:ny2,1] and ha[nx1:nx2,ny1:ny2,1] gt sig_cut * ha_err[nx1:nx2,ny1:ny2,1] and $
        o3[nx1:nx2,ny1:ny2,1] gt sig_cut * o3_err[nx1:nx2,ny1:ny2,1] and hb[nx1:nx2,ny1:ny2,1] gt sig_cut * hb_err[nx1:nx2,ny1:ny2,1] )
	oplot,(alog10(n2[nx1:nx2,ny1:ny2,1]/ha[nx1:nx2,ny1:ny2,1]))[ind],(alog10(o3[nx1:nx2,ny1:ny2,1]/hb[nx1:nx2,ny1:ny2,1]))[ind],psym=4,color=cgcolor('green'),thick=3


nx1=nnx1[7]
nx2=nnx2[7]
ny1=nny1[7]
ny2=nny2[7]

 	ind = where(n2[nx1:nx2,ny1:ny2,1] gt sig_cut * n2_err[nx1:nx2,ny1:ny2,1] and ha[nx1:nx2,ny1:ny2,1] gt sig_cut * ha_err[nx1:nx2,ny1:ny2,1] and $
        o3[nx1:nx2,ny1:ny2,1] gt sig_cut * o3_err[nx1:nx2,ny1:ny2,1] and hb[nx1:nx2,ny1:ny2,1] gt sig_cut * hb_err[nx1:nx2,ny1:ny2,1] )
	oplot,(alog10(n2[nx1:nx2,ny1:ny2,1]/ha[nx1:nx2,ny1:ny2,1]))[ind],(alog10(o3[nx1:nx2,ny1:ny2,1]/hb[nx1:nx2,ny1:ny2,1]))[ind],psym=4,color=cgcolor('green'),thick=3


nx1=nnx1[8]
nx2=nnx2[8]
ny1=nny1[8]
ny2=nny2[8]

 	ind = where(n2[nx1:nx2,ny1:ny2,1] gt sig_cut * n2_err[nx1:nx2,ny1:ny2,1] and ha[nx1:nx2,ny1:ny2,1] gt sig_cut * ha_err[nx1:nx2,ny1:ny2,1] and $
        o3[nx1:nx2,ny1:ny2,1] gt sig_cut * o3_err[nx1:nx2,ny1:ny2,1] and hb[nx1:nx2,ny1:ny2,1] gt sig_cut * hb_err[nx1:nx2,ny1:ny2,1] )
	oplot,(alog10(n2[nx1:nx2,ny1:ny2,1]/ha[nx1:nx2,ny1:ny2,1]))[ind],(alog10(o3[nx1:nx2,ny1:ny2,1]/hb[nx1:nx2,ny1:ny2,1]))[ind],psym=4,color=cgcolor('green'),thick=3




;	ind = where(n2[*,*,2] gt sig_cut * n2_err[*,*,2] and ha[*,*,2] gt sig_cut * ha_err[*,*,2] and $ 
;        o3[*,*,2] gt sig_cut * o3_err[*,*,2] and hb[*,*,2] gt sig_cut * hb_err[*,*,2] )
;	oplot,(alog10(n2[*,*,2]/ha[*,*,2]))[ind],(alog10(o3[*,*,2]/hb[*,*,2]))[ind],psym=4,color=cgcolor('orange')
;	ind = where(n2[*,*,3] gt sig_cut * n2_err[*,*,3] and ha[*,*,3] gt sig_cut * ha_err[*,*,3] and o3[*,*,3] gt sig_cut * o3_err[*,*,3] and hb[*,*,3] gt sig_cut * hb_err[*,*,3] )
;	oplot,(alog10(n2[*,*,3]/ha[*,*,3]))[ind],(alog10(o3[*,*,3]/hb[*,*,3]))[ind],psym=4,color=cgcolor('red'),thick=3
	; BPT curves

	n2bptx1 = FINDGEN(100)/100 - 1.1
	n2bpty1 = 0.61/(n2bptx1-0.05)+1.3
	n2bptx2 = findgen(310)/100 - 2.8
	n2bpty2 = 0.61/(n2bptx2-0.47)+1.19
	oplot,n2bptx1,n2bpty1,linestyle=2 ; Kauffmann et al. (2003) curve
	oplot,n2bptx2,n2bpty2,linestyle=0 ; Kewley et al. (2001) curve

	; S2-Halpha BPT
;	plot,indgen(10),/nodata,xrange=[-1,1],yrange=[-1,0.8],/xsty,/ysty,ytitle='log([OIII]/Hbeta)',xtitle='log([SII]/Halpha)'
;	ind = where(s2[*,*,1] gt sig_cut * s2_err[*,*,1] and ha[*,*,1] gt sig_cut * ha_err[*,*,1] and o3[*,*,1] gt sig_cut * o3_err[*,*,1] and hb[*,*,1] gt sig_cut * hb_err[*,*,1] )
;	oplot,(alog10(s2[*,*,1]/ha[*,*,1]))[ind],(alog10(o3[*,*,1]/hb[*,*,1]))[ind],psym=4,color=cgcolor('blue'),thick=3
;	ind = where(s2[*,*,2] gt sig_cut * s2_err[*,*,2] and ha[*,*,2] gt sig_cut * ha_err[*,*,2] and o3[*,*,2] gt sig_cut * o3_err[*,*,2] and hb[*,*,2] gt sig_cut * hb_err[*,*,2] )
;	oplot,(alog10(s2[*,*,2]/ha[*,*,2]))[ind],(alog10(o3[*,*,2]/hb[*,*,2]))[ind],psym=4,color=cgcolor('orange')
;	ind = where(s2[*,*,3] gt sig_cut * s2_err[*,*,3] and ha[*,*,3] gt sig_cut * ha_err[*,*,3] and o3[*,*,3] gt sig_cut * o3_err[*,*,3] and hb[*,*,3] gt sig_cut * hb_err[*,*,3] )
;	oplot,(alog10(s2[*,*,3]/ha[*,*,3]))[ind],(alog10(o3[*,*,3]/hb[*,*,3]))[ind],psym=4,color=cgcolor('red'),thick=3
;	; BPT curves
;	s2bptx1 = -1.*findgen(200)/100 + 0.32
;	s2bpty1 = 0.72/(s2bptx1-0.32)+1.3
;	s2bptx2 = findgen(40)/100 - 0.3
;	s2bpty2 = 1.89*s2bptx2 + 0.76
;	oplot,s2bptx1,s2bpty1,linestyle=0 ; Kewley et al. (2001) curve
;	oplot,s2bptx2,s2bpty2,linestyle=1 ; Kewley et al. (2006) line

	
	; O1-Halpha BPT
	plot,indgen(10),/nodata,xrange=[-2.5,0],yrange=[-1,0.8],/xsty,/ysty,ytitle='log([OIII]/Hbeta)',xtitle='log([OI]/Halpha)'

nx1=nnx1[0]
nx2=nnx2[0]
ny1=nny1[0]
ny2=nny2[0]

	ind = where(o1[nx1:nx2,ny1:ny2,1] gt sig_cut * o1_err[nx1:nx2,ny1:ny2,1] and ha[nx1:nx2,ny1:ny2,1] gt sig_cut * ha_err[nx1:nx2,ny1:ny2,1] and o3[nx1:nx2,ny1:ny2,1] gt sig_cut * o3_err[nx1:nx2,ny1:ny2,1] and hb[nx1:nx2,ny1:ny2,1] gt sig_cut * hb_err[nx1:nx2,ny1:ny2,1] )
	oplot,(alog10(o1[nx1:nx2,ny1:ny2,1]/ha[nx1:nx2,ny1:ny2,1]))[ind],(alog10(o3[nx1:nx2,ny1:ny2,1]/hb[nx1:nx2,ny1:ny2,1]))[ind],psym=4,color=cgcolor('blue'),thick=3

nx1=nnx1[1]
nx2=nnx2[1]
ny1=nny1[1]
ny2=nny2[1]

	ind = where(o1[nx1:nx2,ny1:ny2,1] gt sig_cut * o1_err[nx1:nx2,ny1:ny2,1] and ha[nx1:nx2,ny1:ny2,1] gt sig_cut * ha_err[nx1:nx2,ny1:ny2,1] and o3[nx1:nx2,ny1:ny2,1] gt sig_cut * o3_err[nx1:nx2,ny1:ny2,1] and hb[nx1:nx2,ny1:ny2,1] gt sig_cut * hb_err[nx1:nx2,ny1:ny2,1] )
	oplot,(alog10(o1[nx1:nx2,ny1:ny2,1]/ha[nx1:nx2,ny1:ny2,1]))[ind],(alog10(o3[nx1:nx2,ny1:ny2,1]/hb[nx1:nx2,ny1:ny2,1]))[ind],psym=4,color=cgcolor('blue'),thick=3

nx1=nnx1[2]
nx2=nnx2[2]
ny1=nny1[2]
ny2=nny2[2]

	ind = where(o1[nx1:nx2,ny1:ny2,1] gt sig_cut * o1_err[nx1:nx2,ny1:ny2,1] and ha[nx1:nx2,ny1:ny2,1] gt sig_cut * ha_err[nx1:nx2,ny1:ny2,1] and o3[nx1:nx2,ny1:ny2,1] gt sig_cut * o3_err[nx1:nx2,ny1:ny2,1] and hb[nx1:nx2,ny1:ny2,1] gt sig_cut * hb_err[nx1:nx2,ny1:ny2,1] )
	oplot,(alog10(o1[nx1:nx2,ny1:ny2,1]/ha[nx1:nx2,ny1:ny2,1]))[ind],(alog10(o3[nx1:nx2,ny1:ny2,1]/hb[nx1:nx2,ny1:ny2,1]))[ind],psym=4,color=cgcolor('blue'),thick=3


nx1=nnx1[3]
nx2=nnx2[3]
ny1=nny1[3]
ny2=nny2[3]


	ind = where(o1[nx1:nx2,ny1:ny2,1] gt sig_cut * o1_err[nx1:nx2,ny1:ny2,1] and ha[nx1:nx2,ny1:ny2,1] gt sig_cut * ha_err[nx1:nx2,ny1:ny2,1] and o3[nx1:nx2,ny1:ny2,1] gt sig_cut * o3_err[nx1:nx2,ny1:ny2,1] and hb[nx1:nx2,ny1:ny2,1] gt sig_cut * hb_err[nx1:nx2,ny1:ny2,1] )
	oplot,(alog10(o1[nx1:nx2,ny1:ny2,1]/ha[nx1:nx2,ny1:ny2,1]))[ind],(alog10(o3[nx1:nx2,ny1:ny2,1]/hb[nx1:nx2,ny1:ny2,1]))[ind],psym=5,color=cgcolor('black'),thick=3


nx1=nnx1[4]
nx2=nnx2[4]
ny1=nny1[4]
ny2=nny2[4]


	ind = where(o1[nx1:nx2,ny1:ny2,1] gt sig_cut * o1_err[nx1:nx2,ny1:ny2,1] and ha[nx1:nx2,ny1:ny2,1] gt sig_cut * ha_err[nx1:nx2,ny1:ny2,1] and o3[nx1:nx2,ny1:ny2,1] gt sig_cut * o3_err[nx1:nx2,ny1:ny2,1] and hb[nx1:nx2,ny1:ny2,1] gt sig_cut * hb_err[nx1:nx2,ny1:ny2,1] )
	oplot,(alog10(o1[nx1:nx2,ny1:ny2,1]/ha[nx1:nx2,ny1:ny2,1]))[ind],(alog10(o3[nx1:nx2,ny1:ny2,1]/hb[nx1:nx2,ny1:ny2,1]))[ind],psym=4,color=cgcolor('red'),thick=3

nx1=nnx1[5]
nx2=nnx2[5]
ny1=nny1[5]
ny2=nny2[5]


	ind = where(o1[nx1:nx2,ny1:ny2,1] gt sig_cut * o1_err[nx1:nx2,ny1:ny2,1] and ha[nx1:nx2,ny1:ny2,1] gt sig_cut * ha_err[nx1:nx2,ny1:ny2,1] and o3[nx1:nx2,ny1:ny2,1] gt sig_cut * o3_err[nx1:nx2,ny1:ny2,1] and hb[nx1:nx2,ny1:ny2,1] gt sig_cut * hb_err[nx1:nx2,ny1:ny2,1] )
	oplot,(alog10(o1[nx1:nx2,ny1:ny2,1]/ha[nx1:nx2,ny1:ny2,1]))[ind],(alog10(o3[nx1:nx2,ny1:ny2,1]/hb[nx1:nx2,ny1:ny2,1]))[ind],psym=4,color=cgcolor('red'),thick=3

nx1=nnx1[6]
nx2=nnx2[6]
ny1=nny1[6]
ny2=nny2[6]


	ind = where(o1[nx1:nx2,ny1:ny2,1] gt sig_cut * o1_err[nx1:nx2,ny1:ny2,1] and ha[nx1:nx2,ny1:ny2,1] gt sig_cut * ha_err[nx1:nx2,ny1:ny2,1] and o3[nx1:nx2,ny1:ny2,1] gt sig_cut * o3_err[nx1:nx2,ny1:ny2,1] and hb[nx1:nx2,ny1:ny2,1] gt sig_cut * hb_err[nx1:nx2,ny1:ny2,1] )
	oplot,(alog10(o1[nx1:nx2,ny1:ny2,1]/ha[nx1:nx2,ny1:ny2,1]))[ind],(alog10(o3[nx1:nx2,ny1:ny2,1]/hb[nx1:nx2,ny1:ny2,1]))[ind],psym=4,color=cgcolor('green'),thick=3

nx1=nnx1[7]
nx2=nnx2[7]
ny1=nny1[7]
ny2=nny2[7]


	ind = where(o1[nx1:nx2,ny1:ny2,1] gt sig_cut * o1_err[nx1:nx2,ny1:ny2,1] and ha[nx1:nx2,ny1:ny2,1] gt sig_cut * ha_err[nx1:nx2,ny1:ny2,1] and o3[nx1:nx2,ny1:ny2,1] gt sig_cut * o3_err[nx1:nx2,ny1:ny2,1] and hb[nx1:nx2,ny1:ny2,1] gt sig_cut * hb_err[nx1:nx2,ny1:ny2,1] )
	oplot,(alog10(o1[nx1:nx2,ny1:ny2,1]/ha[nx1:nx2,ny1:ny2,1]))[ind],(alog10(o3[nx1:nx2,ny1:ny2,1]/hb[nx1:nx2,ny1:ny2,1]))[ind],psym=4,color=cgcolor('green'),thick=3


nx1=nnx1[8]
nx2=nnx2[8]
ny1=nny1[8]
ny2=nny2[8]


	ind = where(o1[nx1:nx2,ny1:ny2,1] gt sig_cut * o1_err[nx1:nx2,ny1:ny2,1] and ha[nx1:nx2,ny1:ny2,1] gt sig_cut * ha_err[nx1:nx2,ny1:ny2,1] and o3[nx1:nx2,ny1:ny2,1] gt sig_cut * o3_err[nx1:nx2,ny1:ny2,1] and hb[nx1:nx2,ny1:ny2,1] gt sig_cut * hb_err[nx1:nx2,ny1:ny2,1] )
	oplot,(alog10(o1[nx1:nx2,ny1:ny2,1]/ha[nx1:nx2,ny1:ny2,1]))[ind],(alog10(o3[nx1:nx2,ny1:ny2,1]/hb[nx1:nx2,ny1:ny2,1]))[ind],psym=4,color=cgcolor('green'),thick=3





;	ind = where(o1[*,*,2] gt sig_cut * o1_err[*,*,2] and ha[*,*,2] gt sig_cut * ha_err[*,*,2] and o3[*,*,2] gt sig_cut * o3_err[*,*,2] and hb[*,*,2] gt sig_cut * hb_err[*,*,2] )
;	oplot,(alog10(o1[*,*,2]/ha[*,*,2]))[ind],(alog10(o3[*,*,2]/hb[*,*,2]))[ind],psym=4,color=cgcolor('orange')
;	ind = where(o1[*,*,3] gt sig_cut * o1_err[*,*,3] and ha[*,*,3] gt sig_cut * ha_err[*,*,3] and o3[*,*,3] gt sig_cut * o3_err[*,*,3] and hb[*,*,3] gt sig_cut * hb_err[*,*,3] )
;	oplot,(alog10(o1[*,*,3]/ha[*,*,3]))[ind],(alog10(o3[*,*,3]/hb[*,*,3]))[ind],psym=4,color=cgcolor('red'),thick=3
	; BPT curves
	o1bptx1 = -1.*findgen(200)/100 - 0.59
	o1bpty1 = 0.73/(o1bptx1+0.59) + 1.33
	o1bptx2 = findgen(60)/100 - 1.1
	o1bpty2 = 1.18 * o1bptx2 + 1.3
	oplot,o1bptx1,o1bpty1,linestyle=0 ; Kewley et al. (2001) curve
	oplot,o1bptx2,o1bpty2,linestyle=1 ; Kewley et al. (2006) line

	; add legend
;	al_legend,['c1','c2','c3'],color=[cgcolor('blue'),cgcolor('orange'),cgcolor('red')],$
;	psym=[4,4,4],charsize=1.5,/top,/left


   if cr eq 'y' then begin
   endps,psname,/post
   print,' Postscript file called ',psname,' has been created '
endif



END

pro ifps,filenm,xsize,ysize,post=post,encapsulated=encapsulated,color=color,$
   r=r,g=g,b=b,landscape=landscape
;
;   FILENM the file you want to create
;   XSIZE  the size in the x-direction ('X': 100pixels 'PS': inches)
;   YSIZE  the size in the y-direction ('X': 100pixels 'PS': inches)
;   POST   set this keyword if you want to create a postscript file
;   ENCAPSULATED  set this keyword if it's an encapsulated ps file
;   COLOR  set this keyword if it's a color file
;   R,G,B  these values will be passed back to the user if you need
;   the
;          color table
;   LANDSCAPE  set this to plot landscape - no guarantees that this
;                                           works
;
;   eg.  ifps,'myfile.ps',6,6,/post,/encapsulated
;
;   This command must be followed by 'endps'

common devs,mydevice

!p.background=255
!p.color=0
!p.multi=0
!p.charsize=1.5
!p.charthick=2.0
!p.thick=1.0

IF keyword_set(ENCAPSULATED) then POST=1
IF keyword_set(POST) THEN BEGIN
   IF keyword_set(LANDSCAPE) then BEGIN
      xoffset = 0.0
      yoffset = 0.0
   ENDIF ELSE BEGIN
      xoffset = 0.0
      yoffset = 0.0
      landscape = 0
   ENDELSE
   mydevice = !D.NAME
   set_plot,'ps'
   device,file=filenm,ysize=ysize,xsize=xsize,xoffset=xoffset,color=color,$
      yoffset=yoffset,/inches,bits_per_pixel=8,encapsulated=encapsulated,$
      landscape=landscape
   device,/symbol,font_index=4
   device,/helvetica,font_index=3
   !p.font=0
   IF keyword_set(COLOR) then TVLCT, r,g,b,/get 
ENDIF ELSE BEGIN
   mydevice = !D.NAME
   set_plot,'x'
   device,decomposed=0
   window,0,xsize=xsize*100,ysize=ysize*100,retain=2
ENDELSE
!p.background=255
!p.color=0
!p.multi=[0,1,1]
!p.charsize=0.8
!p.charthick=2.0
!p.thick=1.0

return
end


pro endps,fname,post=post,pdf=pdf,noview=noview
;;
;;  FNAME  the file that you created
;;  POST   set this keyword if it's a postscript file
;;  PDF    set this keyword if you also want a pdf file
;;  NOVIEW set this keyword if you don't want to view the file
;;
;;  eg.   endps,'myfile.ps',/post,/pdf  
;;        this will create both a 
;;        PS and PDF file, and display the pdf file using Acrobat
;;        Reader
;;  
;;  This command must be preceded by 'ifps'


common devs,mydevice

IF keyword_set(post) THEN BEGIN
   device,/close
   print,'Created PostScript file '+fname
   IF not keyword_set(pdf) THEN BEGIN
      IF not keyword_set(noview) THEN spawn,'ghostview '+fname
   ENDIF ELSE BEGIN
      idot = rstrpos(fname,'.')
      pdffile = strmid(fname,0,idot+1)+'pdf'
      spawn,'epstopdf --outfile='+pdffile+' '+fname
      print,'Created PDF file '+pdffile
      IF not keyword_set(noview) THEN spawn,'acroread '+pdffile
   ENDELSE
print,'leaving ENDPS - DEVICE is '+mydevice
   set_plot,mydevice
   IF MYDEVICE eq 'x' or MYDEVICE eq 'X' THEN device,decomposed=0
ENDIF

return
end

