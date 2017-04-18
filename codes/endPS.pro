pro endps,fname,post=post,pdf=pdf,noview=noview
;;
;;  FNAME  the file that you created
;;  POST   set this keyword if it's a postscript file
;;  PDF    set this keyword if you also want a pdf file
;;  NOVIEW set this keyword if you don't want to view the file
;;
;;  eg.   endps,'myfile.ps',/post,/pdf  
;;        this will create both a 
;;        PS and PDF file, and display the pdf file using Acrobat Reader
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
