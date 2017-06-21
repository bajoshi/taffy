PRO s2hao3hb_em_m,ps=ps,oplot=oplot,legend=legend,age=age,inst=inst

if ~keyword_set(age) then age='0.0'
if ~keyword_set(inst) then inst='inst'

if ~keyword_set(ps) then ps=0
if ~keyword_set(oplot) then oplot=0
if ~keyword_set(legend) then legend=0

readcol,'/Users/jrich/idl/'+inst+'_'+age+'.txt',mets_ph,qs_ph,oiia_ph,oiib_ph,oiii_ph,oi_ph,niia_ph,ha_ph,niib_ph,siia_ph,siib_ph,arii_ph,ariii_ph,siv_ph,neii_ph,neiii_ph,siii18_ph,oiv_ph,feii_ph,siii33_ph,slii_ph,niii_ph,oi63a_ph,oi63b_ph,oiii88_ph,nii122_ph,oi146a_ph,oi146b_ph,cii_ph,nii205_ph,/silent

sii_ph = siia_ph+siib_ph


q1e7=where(mets_ph ge 8.29 AND qs_ph eq 1e7) & q2e7=where(mets_ph ge 8.29 AND qs_ph eq 2e7) & q4e7=where(mets_ph ge 8.29 AND qs_ph eq 4e7) & q8e7=where(mets_ph ge 8.29 AND qs_ph eq 8e7) & q1e8=where(mets_ph ge 8.29 AND qs_ph eq 1e8) & q2e8=where(mets_ph ge 8.29 AND qs_ph eq 2e8) & q4e8=where(mets_ph ge 8.29 AND qs_ph eq 4e8)

mp739=where(mets_ph eq 7.39) & mp799=where(mets_ph eq 7.99) & mp829=where(mets_ph eq 8.29) & mp869=where(mets_ph eq 8.69) & mp899=where(mets_ph eq 8.99)

zph = ['8.99','8.69','8.29'];,'7.99','7.39']
qph = ['4e8','2e8','1e8','8e7','4e7','2e7','1e7'];,'6.75','6.50']
metcols = [155,140,125,110,95,80,65,50,35]
vcols   = [255,245,235,225,220,215]
metcolp = [125,110,95];,65,35]
qcolp   = [255,250,245,240,235,230,225];,220,215]

OIIINII_ph=alog10(OIII_ph/NIIb_ph)
NIISII_ph = alog10(NIIb_ph/SII_ph)
NIIHA_ph  = alog10(NIIb_ph/Ha_ph)
SIIHA_ph  = alog10(SII_ph/HA_ph)
OIIIHB_ph = alog10(OIII_ph)
OIHA_ph = alog10(OI_ph/Ha_ph)
neiiineii_ph=alog10(neiii_ph/neii_ph)


x=alog10(sii_ph/ha_ph)
y=alog10(oiii_ph)


loadct,13,/silent

alli = [mp829,mp869,mp899]
if oplot then begin
    oplot,x[alli],y[alli],psym=3
endif else begin
    plot,x[alli],y[alli],psym=3
endelse
;plot lines of constant q
;    oplot,x[mp739],y[mp739],thick=5,color=metcolp[4]
;    oplot,x[mp799],y[mp799],thick=5,color=metcolp[3]
    oplot,x[mp829],y[mp829],thick=5,color=metcolp[2]
    oplot,x[mp869],y[mp869],thick=5,color=metcolp[1]
    oplot,x[mp899],y[mp899],thick=5,color=metcolp[0]
;    oplot,x[mp],y[mp],thick=5,color=metcolp[5]
;    oplot,x[mp],y[mp],thick=5,color=metcolp[6]
;    oplot,x[mp],y[mp],thick=5,color=metcolp[7]
;    oplot,x[mp],y[mp],thick=5,color=metcolp[8]
;plot lines of constant z

    loadct,13,/silent
    oplot,x[q4e8],y[q4e8],thick=5,color=qcolp[0]
    oplot,x[q2e8],y[q2e8],thick=5,color=qcolp[1]
    oplot,x[q1e8],y[q1e8],thick=5,color=qcolp[2]
    oplot,x[q8e7],y[q8e7],thick=5,color=qcolp[3]
    oplot,x[q4e7],y[q4e7],thick=5,color=qcolp[4]
    oplot,x[q2e7],y[q2e7],thick=5,color=qcolp[5]
    oplot,x[q1e7],y[q1e7],thick=5,color=qcolp[6]
;    oplot,x[qind675],y[qind675],thick=5,color=qcolp[7]
;    oplot,x[qind650],y[qind650],thick=5,color=qcolp[8]
;


if legend then begin
    legend,zph,/left,/bottom,colors=metcolp,linestyle=(zph*0),thick=(zph*0)+5
    legend,qph,/right,colors=qcolp,linestyle=(qph*0),thick=(qph*0)+5
endif

loadct,0,/silent

END
