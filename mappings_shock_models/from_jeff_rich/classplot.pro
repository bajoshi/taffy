PRO classplot,line,xmin,xmax,ymin,ymax,thick=thick,xthick=xthick,ythick=ythick,charthick=charthick,charsize=charsize,ytitle=ytitle,xtitle=xtitle,color=color,background=background,position=position,xstyle=xstlye,ystyle=ystyle,kauff=kauff,oplot=oplot,ccolor=ccolor,lcolor=lcolor,meterr=meterr,font=font,divstyle=divstyle
  
;Plots axes/max starburst line for BPT diagrams for NII (also
;composite line), SII or OI 

if n_params() lt 1 then begin
	print, ' Syntax:    classplot, line, [xmin, xmax, ymin, ymax]'
	print, '                         line = NII, SII or OI'
	print, '                         xmin/max,ymin/max=range'
	print, ''
	return
endif

if ~keyword_set(oplot) then oplot=0
if ~keyword_set(meterr) then meterr=0
if ~keyword_set(font) then font=0


;x/ystyle default to 1
if ~keyword_set(xstyle) then xstyle=1
if ~keyword_set(ystyle) then ystyle=1
if ~keyword_set(kauff)  then kauff=0

;define y range if not defined
if n_elements(ymin) eq 0 then ymin = -1.2
if n_elements(ymax) eq 0 then ymax = 1.6
if n_elements(thick) eq 0 then thick=1
if n_elements(xthick) eq 0 then xthick=1
if n_elements(ythick) eq 0 then ythick=1
if n_elements(charthick) eq 0 then charthick=1
if n_elements(charsize) eq 0 then charsize=1

;div line seyfert/liner
if n_elements(divstyle) eq 0 then divstyle=2


if (oplot eq 0) then begin

if (line eq 'NII') then begin
    
    ;define x range if not defined
    if n_elements(xmin) eq 0 then xmin = -2.0
    if n_elements(xmax) eq 0 then xmax = 1.0

    x = (findgen(1000)*(xmax-xmin)/1000.)+xmin


    ind = where(x lt 0.47)
    ind2= where(x lt 0.05)

    y = (0.61 / (x[ind]-0.47))+1.19
    y2= (0.61 / (x[ind2]-0.05))+1.30
    y3= (0.61 / (x[ind2]-0.15))+1.30

    if ~ keyword_set(xtitle) then xtitle='log[NII/H!7a!3]'
    if ~ keyword_set(ytitle) then ytitle='log[OIII/H!7b!3]'


    plot,x,y,xrange=[xmin,xmax],yrange=[ymin,ymax],xstyle=xstyle,ystyle=ystyle,xtitle=xtitle,ytitle=ytitle,thick=thick,xthick=xthick,ythick=ythick,charthick=charthick,charsize=charsize,xcharsize=charsize,ycharsize=charsize,background=background,color=ccolor,position=position,font=font
    if kauff eq 1 then oplot,x,y2,max_value=0.842,linestyle=2,thick=thick
    if meterr eq 1 then oplot,x,y3,max_value=0.842,linestyle=2,thick=thick

endif

if (line eq 'SII') then begin

    ;define x range if not defined
    if n_elements(xmin) eq 0 then xmin = -1.2
    if n_elements(xmax) eq 0 then xmax = 0.7
    
    x = (findgen(1000)*(xmax-xmin)/1000.)+xmin

    
    ind = where(x lt 0.32)

    y = (0.72 / (x[ind]-0.32))+1.30
    
    y2= (1.89*x)+0.76

    if ~ keyword_set(xtitle) then xtitle='log[SII/H!7a!3]'
    if ~ keyword_set(ytitle) then ytitle='log[OIII/H!7b!3]'


    plot,x,y,xrange=[xmin,xmax],yrange=[ymin,ymax],xstyle=xstyle,ystyle=ystyle,xtitle=xtitle,ytitle=ytitle,thick=thick,xthick=xthick,ythick=ythick,charthick=charthick,charsize=charsize,xcharsize=charsize,ycharsize=charsize,background=background,color=ccolor,position=position,font=font
    oplot,x,y2,min_value=0.17,max_value=0.8,linestyle=divstyle,thick=thick,color=lcolor


endif

if (line eq 'OI') then begin
    
    ;define x range if not defined
    if n_elements(xmin) eq 0 then xmin = -2.2
    if n_elements(xmax) eq 0 then xmax = 0.2

    x = (findgen(1000)*(xmax-xmin)/1000.)+xmin

    ind = where(x lt -0.59)

    y = (0.73 / (x[ind]+0.59))+1.33
    y2= (1.18*x)+1.30

    if ~ keyword_set(xtitle) then xtitle='log[OI/H!7a!3]'
    if ~ keyword_set(ytitle) then ytitle='log[OIII/H!7b!3]'

    plot,x,y,xrange=[xmin,xmax],yrange=[ymin,ymax],xstyle=xstyle,ystyle=ystyle,xtitle=xtitle,ytitle=ytitle,thick=thick,xthick=xthick,ythick=ythick,charthick=charthick,charsize=charsize,xcharsize=charsize,ycharsize=charsize,background=background,color=ccolor,position=position,font=font
    oplot,x,y2,min_value=-0.025,max_value=0.8,linestyle=divstyle,thick=thick,color=lcolor

endif

endif


if (oplot eq 1) then begin

if (line eq 'NII') then begin
    
    ;define x range if not defined
    if n_elements(xmin) eq 0 then xmin = -2.0
    if n_elements(xmax) eq 0 then xmax = 1.0

    x = (findgen(1000)*(xmax-xmin)/1000.)+xmin


    ind = where(x lt 0.47)
    ind2= where(x lt 0.05)

    y = (0.61 / (x[ind]-0.47))+1.19
    y2= (0.61 / (x[ind2]-0.05))+1.30

    if ~ keyword_set(xtitle) then xtitle='log[NII/H!7a!3]'
    if ~ keyword_set(ytitle) then ytitle='log[OIII/H!7b!3]'


    oplot,x,y,color=ccolor
    if kauff eq 1 then oplot,x,y2,max_value=0.842,linestyle=2,thick=thick

endif

if (line eq 'SII') then begin

    ;define x range if not defined
    if n_elements(xmin) eq 0 then xmin = -1.2
    if n_elements(xmax) eq 0 then xmax = 0.7
    
    x = (findgen(1000)*(xmax-xmin)/1000.)+xmin

    
    ind = where(x lt 0.32)

    y = (0.72 / (x[ind]-0.32))+1.30
    
    y2= (1.89*x)+0.76

    if ~ keyword_set(xtitle) then xtitle='log[SII/H!7a!3]'
    if ~ keyword_set(ytitle) then ytitle='log[OIII/H!7b!3]'


    oplot,x,y,color=ccolor
    oplot,x,y2,min_value=0.17,max_value=0.8,linestyle=2,thick=thick,color=lcolor


endif

if (line eq 'OI') then begin
    
    ;define x range if not defined
    if n_elements(xmin) eq 0 then xmin = -2.2
    if n_elements(xmax) eq 0 then xmax = 0.2

    x = (findgen(1000)*(xmax-xmin)/1000.)+xmin

    ind = where(x lt -0.59)

    y = (0.73 / (x[ind]+0.59))+1.33
    y2= (1.18*x)+1.30

    if ~ keyword_set(xtitle) then xtitle='log[OI/H!7a!3]'
    if ~ keyword_set(ytitle) then ytitle='log[OIII/H!7b!3]'

    oplot,x,y,color=ccolor
    oplot,x,y2,min_value=-0.025,max_value=0.8,linestyle=2,thick=thick,color=lcolor

endif

endif



END
