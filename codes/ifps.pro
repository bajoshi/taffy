pro ifps,filenm,xsize,ysize,post=post,encapsulated=encapsulated,color=color,$
   r=r,g=g,b=b,landscape=landscape
;
;   FILENM the file you want to create
;   XSIZE  the size in the x-direction ('X': 100pixels 'PS': inches)
;   YSIZE  the size in the y-direction ('X': 100pixels 'PS': inches)
;   POST   set this keyword if you want to create a postscript file
;   ENCAPSULATED  set this keyword if it's an encapsulated ps file
;   COLOR  set this keyword if it's a color file
;   R,G,B  these values will be passed back to the user if you need the
;          color table
;   LANDSCAPE  set this to plot landscape - no guarantees that this works
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
