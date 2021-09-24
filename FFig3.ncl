load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin

;PDRMIP
models=(/"CanESM2","GISS-E2-R","HadGEM2","HadGEM3","IPSL-CM5A","MIROC-SPRINTARS","MPI-ESM","NCAR-CESM1-CAM5","NorESM1"/)
nmodels=dimsizes(models)


sfilter=30
sfilter2=150
calc_lag=False
Lsfile1="/Users/stergios/lpool/CESM-LEM/SOL/CESM-LEM_SOL_TS_850-2005_fltbwPAD"+sfilter2+"_fldm_ym.nc"
Lsfile2="/Users/stergios/lpool/CESM-LEM/SOL/CESM-LEM_SOL_TS_850-2005_fltbwPAD"+sfilter2+"_ym.nc"
Lsfile3="/Users/stergios/lpool/CESM-LEM/SOL/CESM-LEM_SOL_PREC_850-2005_fltbwPAD"+sfilter2+"_ym.nc"
Lsfile5="/Users/stergios/lpool/CESM-LEM/SOL/CESM-LEM_SOL_PREC_850-2005_ym.nc"

Lcfile1="/Users/stergios/lpool/CESM-LEM/GHG/CESM-LEM_GHG_TS_1850-2005_fltbwPAD30_fldm_ym.nc"
Lcfile2="/Users/stergios/lpool/CESM-LEM/GHG/CESM-LEM_GHG_TS_1850-2005_fltbwPAD30_ym.nc"
Lcfile3="/Users/stergios/lpool/CESM-LEM/GHG/CESM-LEM_GHG_PREC_1850-2005_fltbwPAD30_ym.nc"

Lbfile="/Users/stergios/lpool/CESM-LEM/GHG/CESM-LEM_GHG_PREC_900-1200_tm.nc"

Lsfile4="/Users/stergios/lpool/CESM-LEM/CESM-LEM_tsi_850-2005_fltbwPAD"+sfilter2+"_ym.nc"
Lcfile4="/Users/stergios/lpool/CESM-LEM/CESM-LEM_CO2_850-2007_fltbwPAD30_ym.nc"

lmfile="/Users/stergios/lpool/obs/noaa-lmrv2.1/noaa-lmrv2.1_air_850-1800_ym.nc"
tsifile="/Users/stergios/lpool/CESM-LEM/CESM-LEM_tsi_850-1800_fltbwPAD150_ym.nc"

graph1="/Users/stergios/analyses/19_pdrmip/figures/FFig3"

graph3="/Users/stergios/analyses/19_pdrmip/figures/Efig3_pr_regtsi0"
graph4="/Users/stergios/analyses/19_pdrmip/figures/FFigSup_pr_MCALIA"

graph2="/Users/stergios/analyses/19_pdrmip/figures/FFigSup_precAll"
 rad   = 4.0*atan(1.0)/180.0
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

;read TSI
;=========
print("read "+Lsfile4)
ff1=addfile(Lsfile4,"r")
tsi=ff1->tsi(:,0,0)
tsi_orig=(/tsi-tsi(0)/)
tsi=(/tsi-avg(tsi)/)
delete(ff1)

;read TSI
;=========
print("read "+Lcfile4)
ff1=addfile(Lcfile4,"r")
co2=ff1->CO2
co2_orig=(/co2-co2(0)/)
co2=(/co2-avg(co2)/)
delete(ff1)

; READ PDRMIP
;++++++++++++++++++++++++++++++
sol_tas_fldm_all=new((/nmodels,100/),double)
co2_tas_fldm_all=new((/nmodels,100/),double)
sol_pr_all=new((/nmodels,100,60,90/),double)
co2_pr_all=new((/nmodels,100,60,90/),double)
sol_tas_all=new((/nmodels,100,60,90/),double)
co2_tas_all=new((/nmodels,100,60,90/),double)

do ii=0,nmodels-1
model=models(ii)
sfile1="/Users/stergios/lpool/pdrmip/"+model+"/tas_"+model+"_solar_minus_base_coupled_100yr_fltbw"+sfilter+"_timeseries_fldm.nc"
cfile1="/Users/stergios/lpool/pdrmip/"+model+"/tas_"+model+"_co2_minus_base_coupled_100yr_fltbw"+sfilter+"_timeseries_fldm.nc"

sfile2="/Users/stergios/lpool/pdrmip/"+model+"/pr_"+model+"_SOL_minus_base_coupled_100yr_fltbw"+sfilter+"_timeseries.nc"
cfile2="/Users/stergios/lpool/pdrmip/"+model+"/pr_"+model+"_CO2_minus_base_coupled_100yr_fltbw"+sfilter+"_timeseries.nc"

sfile3="/Users/stergios/lpool/pdrmip/"+model+"/tas_"+model+"_solar_minus_base_coupled_100yr_fltbw"+sfilter+"_timeseries.nc"
cfile3="/Users/stergios/lpool/pdrmip/"+model+"/tas_"+model+"_CO2_minus_base_coupled_100yr_fltbw"+sfilter+"_timeseries.nc"

bfile1="/Users/stergios/lpool/pdrmip/"+model+"/tas_"+model+"_base_coupled_100yr_timeseries.nc"
bfile2="/Users/stergios/lpool/pdrmip/"+model+"/pr_"+model+"_base_coupled_100yr_timeseries.nc"

 ;read PDRMIP solar
;==================
print("read "+sfile1)
ff1=addfile(sfile1,"r")
sol_tas_fldm=ff1->tas(:,0,0)
delete(ff1)

print("read "+sfile2)
ff1=addfile(sfile2,"r")
sol_pr=ff1->pr
delete(ff1)

print("read "+sfile3)
ff1=addfile(sfile3,"r")
sol_tas=ff1->tas
delete(ff1)

;read CO2
;===================
print("read "+cfile1)
ff1=addfile(cfile1,"r")
co2_tas_fldm=ff1->tas(:,0,0)
delete(ff1)

print("read "+cfile2)
ff1=addfile(cfile2,"r")
co2_pr=ff1->pr
delete(ff1)

print("read "+cfile3)
ff1=addfile(cfile3,"r")
co2_tas=ff1->tas
delete(ff1)


;read climatology
;===================
print("read "+bfile1)
ff1=addfile(bfile1,"r")
dd1=ff1->tas
clat   = cos(dd1&lat*rad)
base_tas=wgt_areaave_Wrap(dd1, clat, 1.0, 1)
base_tas_tm=dim_avg_n_Wrap(base_tas,0)
delete(clat)
delete(dd1)
delete(ff1)


print("read "+bfile2)
ff1=addfile(bfile2,"r")
base_pr=ff1->pr
dd2=ff1->pr(0,:,:)
base_pr_tm=dim_avg_n_Wrap(base_pr,0)
delete(ff1)


sol_tas_fldm_all(ii,:)=sol_tas_fldm
co2_tas_fldm_all(ii,:)=co2_tas_fldm

sol_pr_all(ii,:,:,:)=sol_pr
co2_pr_all(ii,:,:,:)=co2_pr

sol_tas_all(ii,:,:,:)=sol_tas
co2_tas_all(ii,:,:,:)=co2_tas

;end do
delete(sol_tas)
delete(co2_tas)
delete(sol_tas_fldm)
delete(co2_tas_fldm)
delete(sol_pr)
delete(co2_pr)
print("++++++++++++++++++++++++")
end do

;scale to mm/day
sol_pr_all=(/86400*sol_pr_all/)
co2_pr_all=(/86400*co2_pr_all/)


;PDRMIP ensemble mean
sol_tas_em=dim_avg_n_Wrap(sol_tas_all,0)
sol_tas_fldm_em=dim_avg_n_Wrap(sol_tas_fldm_all,0)
sol_pr_em=dim_avg_n_Wrap(sol_pr_all,0)

co2_tas_em=dim_avg_n_Wrap(co2_tas_all,0)
co2_tas_fldm_em=dim_avg_n_Wrap(co2_tas_fldm_all,0)
co2_pr_em=dim_avg_n_Wrap(co2_pr_all,0)

; regress precipiatation
sol_rg_em=regCoef_n(sol_tas_fldm_em,sol_pr_em,0,0)
co2_rg_em=regCoef_n(co2_tas_fldm_em,co2_pr_em,0,0)

; regress surface tempeature
sol_rg2_em=regCoef_n(sol_tas_fldm_em,sol_tas_em,0,0)
co2_rg2_em=regCoef_n(co2_tas_fldm_em,co2_tas_em,0,0)


copy_VarCoords(dd2,sol_rg_em)
copy_VarCoords(dd2,co2_rg_em)
copy_VarCoords(dd2,sol_rg2_em)
copy_VarCoords(dd2,co2_rg2_em)

;++++++++++++++++++++++++++++++++++++++++++++++++++
; CESM SECTION
;+++++++++++++++++++++++++++++++++++++++++++++++
;read CESM
;===========
print("read "+Lsfile1)
ff1=addfile(Lsfile1,"r")
lem_sol_tas_fldm=ff1->TS(:,0,0)
delete(ff1)

print("read "+Lsfile2)
ff1=addfile(Lsfile2,"r")
lem_sol_tas=ff1->TS
delete(ff1)

print("read "+Lsfile3)
ff1=addfile(Lsfile3,"r")
lem_sol_pr=ff1->PREC
delete(ff1)

;read unfilted
print("read "+Lsfile5)
ff1=addfile(Lsfile5,"r")
lem_sol_pr_unflt=ff1->PRECC
delete(ff1)

;read CO2
;==========
print("read "+Lcfile1)
ff1=addfile(Lcfile1,"r")
lem_co2_tas_fldm=ff1->TS(:,0,0)
delete(ff1)

print("read "+Lcfile2)
ff1=addfile(Lcfile2,"r")
lem_co2_tas=ff1->TS
delete(ff1)

print("read "+Lcfile3)
ff1=addfile(Lcfile3,"r")
lem_co2_pr=ff1->PREC
delete(ff1)

;read model climatology
;===================
print("read "+Lbfile)
ff1=addfile(Lbfile,"r")
lem_pr_clim=ff1->PRECC(0,:,:)
delete(ff1)


; convert to mm/day
lem_sol_pr=(/1000*lem_sol_pr*86400/)
lem_co2_pr=(/1000*lem_co2_pr*86400/)
lem_sol_pr_unflt=(/1000*lem_sol_pr_unflt*86400/)

lem_sol_tas_rc           = regCoef_n(lem_sol_tas_fldm ,lem_sol_tas,0,0)
lem_sol_pr_rc           = regCoef_n(lem_sol_tas_fldm ,lem_sol_pr,0,0)
copy_VarCoords(lem_sol_tas(0,:,:),lem_sol_tas_rc)
copy_VarCoords(lem_sol_tas(0,:,:),lem_sol_pr_rc)

lem_co2_tas_rc           = regCoef_n(lem_co2_tas_fldm ,lem_co2_tas,0,0)
lem_co2_pr_rc           = regCoef_n(lem_co2_tas_fldm ,lem_co2_pr,0,0)

copy_VarCoords(lem_co2_tas(0,:,:),lem_co2_tas_rc)
copy_VarCoords(lem_co2_pr(0,:,:),lem_co2_pr_rc)

;regress to solar cycle
lem_tsi_sol_tas_rc           = regCoef_n(tsi,lem_sol_tas,0,0)
lem_tsi_sol_pr_rc           = regCoef_n(tsi,lem_sol_pr,0,0)
copy_VarCoords(lem_co2_tas(0,:,:),lem_tsi_sol_tas_rc)
copy_VarCoords(lem_co2_pr(0,:,:),lem_tsi_sol_pr_rc)


;calculate lagged regression s coefficients
nlagr = 30 ; 0,1,2
rlag = new ((/6,96,144/),"double")
rlag_pr = rlag
Nend = dimsizes(tsi)-1
il=0
if calc_lag then
do n=0,nlagr-1,5
  print ("lag: "+n)
  rlag(il,:,:) = regCoef_n(tsi(0:Nend-n), lem_sol_tas(n:Nend,:,:),0,0)
  rlag_pr(il,:,:) = regCoef_n(tsi(0:Nend-n), lem_sol_pr(n:Nend,:,:),0,0)
  il=il+1
end do
end if
copy_VarCoords(lem_co2_tas(0,:,:),rlag(0,:,:))
copy_VarCoords(lem_co2_pr(0,:,:),rlag_pr(0,:,:))

; calculate differences between the MCA-LIA
lem_sol_tas&time=ispan(850,2005,1)
lem_sol_pr&time=ispan(850,2005,1)
lem_sol_pr_unflt&time=ispan(850,2005,1)

tas_MCALIA=lem_sol_tas(0,:,:)
pr_MCALIA=lem_sol_pr(0,:,:)
pr_MCALIA_unflt=lem_sol_pr_unflt(0,:,:)

tas_MCALIA=(/dim_avg_n_Wrap(lem_sol_tas({1100:1200},:,:),0)-dim_avg_n_Wrap(lem_sol_tas({1630:1730},:,:),0)/)
pr_MCALIA=(/dim_avg_n_Wrap(lem_sol_pr({1100:1200},:,:),0)-dim_avg_n_Wrap(lem_sol_pr({1630:1730},:,:),0)/)
pr_MCALIA_unflt=(/dim_avg_n_Wrap(lem_sol_pr_unflt({1100:1200},:,:),0)-dim_avg_n_Wrap(lem_sol_pr_unflt({1630:1730},:,:),0)/)

;caclulate significance
ddd1=lem_sol_pr_unflt({1100:1200},:,:)
ddd2=lem_sol_pr_unflt({1630:1730},:,:)
;print(dimsizes(ddd1))
sigr = 0.05
 xtmp = ddd1(lat|:,lon|:,time|:)
xEqv = equiv_sample_size (xtmp, sigr,0)
;print(xEqv)
;exit
 ytmp = ddd2(lat|:,lon|:,time|:)
yEqv = equiv_sample_size (ytmp, sigr,0)
pr_ttest_unf=ddd1(0,:,:)
pr_ttest_unf=(/ttest(dim_avg_n_Wrap(ddd1,0),dim_variance_n_Wrap(ddd1,0),xEqv,dim_avg_n_Wrap(ddd2,0),dim_variance_n_Wrap(ddd2,0),yEqv,False,False)/)
;pr_ttest_unf=(/ttest(dim_avg_n_Wrap(ddd1,0),dim_variance_n_Wrap(ddd1,0),101,dim_avg_n_Wrap(ddd2,0),dim_variance_n_Wrap(ddd2,0),101,False,False)/)


; read LEM reconstruction
print("read "+lmfile)
ff1=addfile(lmfile,"r")
lm_tas=ff1->air
delete(ff1)

print("read "+tsifile)
ff1=addfile(tsifile,"r")
lm_tsi=ff1->tsi(:,0,0)
delete(ff1)

lm_tsi_tas_rc           = regCoef_n(lm_tsi,lm_tas,0,0)
copy_VarCoords(lm_tas(0,:,:),lm_tsi_tas_rc (:,:))

print("=================")
print("plot "+graph1+".png")
plot = new(4,graphic)
plot_tt = new(4,graphic)
wks = gsn_open_wks("png",graph1+"dd")

res1                      = True
res1@gsnDraw              = False        ; don't draw yet
res1@gsnFrame             = False        ; don't advance frame yet

res1@gsnAddCyclic         = True        ; plotted dataa are not cyclic
res1@cnFillOn             = True         ; turn on color fill
res1@cnLinesOn            = False        ; True is default
res1@cnLineLabelsOn       = False        ; turn on line labels
res1@gsnContourZeroLineThicknessF = 3.5 
res1@mpCenterLonF         = 180.               ; center plot at 180
res1@mpMinLatF            = -30         ; zoom in on map
res1@mpMaxLatF            = 30
res1@mpMinLonF            = 40         ; zoom in on map
res1@mpMaxLonF            = 290

res1@tiMainFontHeightF  = 0.016
res1@tmXBLabelFontHeightF = 0.016
res1@tmYLLabelFontHeightF = 0.016
res1@lbLabelBarOn                 = False
res1@mpFillOn             = False        ; turn off map fill
res1@mpCenterLonF         = 180.               ; center plot at 180

res1@mpGeophysicalLineThicknessF =2.
 res1@cnLevelSelectionMode = "ManualLevels"
 res1@cnFillPalette   ="CBR_drywet"
res1@cnMinLevelValF          = -1.              ; min level
res1@cnMaxLevelValF          = 1.           ; max level
res1@cnLevelSpacingF         = .2                  ; contour level spacing
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

res2                      = True
res2@gsnDraw              = False        ; don't draw yet
res2@gsnFrame             = False        ; don't advance frame yet
res2@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res2@cnMinLevelValF       = -1.4              ; set min contour level
  res2@cnMaxLevelValF       =  1.4              ; set max contour level
  res2@cnLevelSpacingF      =  .2              ; set contour spacing
  res2@cnLineLabelsOn       = True
   res2@cnFillPalette   ="posneg_1"
   ;"hotcold_18lev"
res2@cnLineThicknessF     = 3                ; 3x default line thickness
  res2@cnMonoLineColor      = False            ; Use multiple colors
  res2@cnInfoLabelOn   = False 

res_ttest = True                            ; res2 probability plots
res_ttest@gsnDraw             = False       ; Do not draw plot
res_ttest@gsnFrame            = False       ; Do not advance frome

res_ttest@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
res_ttest@cnMinLevelValF      = 0.00        ; set min contour level
res_ttest@cnMaxLevelValF      = 1.        ; set max contour level
res_ttest@cnLevelSpacingF     = 0.01        ; set contour spacing
res_ttest@cnInfoLabelOn       = False
res_ttest@cnLinesOn           = False       ; do not draw contour lines
res_ttest@cnLineLabelsOn      = False       ; do not draw contour labels
res_ttest@cnFillScaleF        = 1.         ; add extra density
res_ttest@cnFillDotSizeF=0.002
res_ttest@gsnRightString               = ""
res_ttest@gsnLeftString                = ""
res_ttest@tiYAxisString   = ""
res_ttest@tiXAxisString   = ""




;PANEL 1
;+++++++++
res1@gsnLeftString                = "PDRMIP-SOL"
plot(0)=gsn_csm_contour_map(wks, sol_rg_em,res1)
plot1_sst  = gsn_csm_contour(wks, sol_rg2_em, res2 )          ; draw second plot
overlay(plot(0),plot1_sst)
dd1= sol_rg_em
dd1=(/onedtond(abs(sol_rg_em@tval),(/dimsizes(dd1&lat),dimsizes(dd1&lon)/))/)
prob =dd1
nlons=dimsizes(prob&lon)
nlats=dimsizes(prob&lat)
df=1 ;?? correct
do ilat=0,nlats-1
do  ilon= 0, nlons-1
  prob(ilat,ilon) = (/betainc( df/(df+dd1(ilat,ilon)^2), df/2.0, 0.5)/)
end do
end do

plot1 = gsn_csm_contour(wks,prob, res_ttest)
plot1   = ShadeLtContour(plot1, 0.05, 17)
overlay (plot(0), plot1)
delete(dd1)
delete(prob)
;plot1 = gsn_csm_contour(wks,dd1, res_ttest)
;plot1   = ShadeGtContour(plot1, 2.6, 17)
;overlay (plot(0), plot1)
;delete(dd1)

;PANEL 2
;+++++++++
res1@gsnLeftString                = "PDRMIP-CO2"
plot(1)=gsn_csm_contour_map(wks, co2_rg_em,res1)
plot2_sst  = gsn_csm_contour(wks, co2_rg2_em, res2 )          ; draw second plot
overlay(plot(1),plot2_sst)
dd1= co2_rg_em
dd1=(/onedtond(abs(co2_rg_em@tval),(/dimsizes(dd1&lat),dimsizes(dd1&lon)/))/)

prob =dd1 
nlons=dimsizes(prob&lon)
nlats=dimsizes(prob&lat)
df=1 ;?? correct
do ilat=0,nlats-1
do  ilon= 0, nlons-1
  prob(ilat,ilon) = (/betainc( df/(df+dd1(ilat,ilon)^2), df/2.0, 0.5)/)
end do
end do

plot2 = gsn_csm_contour(wks,prob, res_ttest)
plot2   = ShadeLtContour(plot2, 0.05, 17)
overlay (plot(1), plot2)
delete(prob)
delete(dd1)

;plot2 = gsn_csm_contour(wks,dd1, res_ttest)
;plot2   = ShadeGtContour(plot2, 2.6, 17)
;overlay (plot(1), plot2)
;delete(dd1)

;PANEL 3
;+++++++++
res1@gsnLeftString                = "CESM-SOL"
plot(2)=gsn_csm_contour_map(wks, lem_sol_pr_rc,res1)
plot3_sst=gsn_csm_contour(wks, lem_sol_tas_rc,res2)
overlay(plot(2),plot3_sst)
dd1= lem_sol_pr_rc
dd1=(/onedtond(abs(lem_sol_pr_rc@tval),(/dimsizes(dd1&lat),dimsizes(dd1&lon)/))/)

prob =dd1 
nlons=dimsizes(prob&lon)
nlats=dimsizes(prob&lat)
df=1 ;?? correct
do ilat=0,nlats-1
do  ilon= 0, nlons-1
  prob(ilat,ilon) = (/betainc( df/(df+dd1(ilat,ilon)^2), df/2.0, 0.5)/)
end do
end do

plot3 = gsn_csm_contour(wks,prob, res_ttest)
plot3   = ShadeLtContour(plot3, 0.05, 17)
overlay (plot(2), plot3)
delete(prob)
delete(dd1)

;plot3 = gsn_csm_contour(wks,dd1, res_ttest)
;plot3   = ShadeGtContour(plot3, 2.6, 17)
;overlay (plot(2), plot3)
;delete(dd1)

;PANEL 4
;+++++++++
res1@gsnLeftString                = "CESM-GHG"
plot(3)=gsn_csm_contour_map(wks, lem_co2_pr_rc,res1)
plot4_sst=gsn_csm_contour(wks, lem_co2_tas_rc,res2)
overlay(plot(3),plot4_sst)
dd1= lem_co2_pr_rc
dd1=(/onedtond(abs(lem_co2_pr_rc@tval),(/dimsizes(dd1&lat),dimsizes(dd1&lon)/))/)
prob =dd1 
nlons=dimsizes(prob&lon)
nlats=dimsizes(prob&lat)
df=1 ;?? correct
do ilat=0,nlats-1
do  ilon= 0, nlons-1
  prob(ilat,ilon) = (/betainc( df/(df+dd1(ilat,ilon)^2), df/2.0, 0.5)/)
end do
end do

plot4 = gsn_csm_contour(wks,prob, res_ttest)
plot4   = ShadeLtContour(plot4, 0.05, 17)
overlay (plot(3), plot4)
delete(prob)
delete(dd1)

;plot4 = gsn_csm_contour(wks,dd1, res_ttest)
;plot4   = ShadeGtContour(plot4, 2.6, 17)
;overlay (plot(3), plot4)
;delete(dd1)


resP                     = True         ; modify the panel plot
resP@gsnPanelFigureStrings= (/"a)","b)","c)","d)","e)","f)","g)","h)","i)","j)","k)","l)"/)
resP@gsnPanelFigureStringsFontHeightF=0.01
resP@amJust     = "TopLeft"
;resP@lbLabelFontHeightF =.015
resP@gsnPanelLabelBar = True
resP@gsnMaximize         = True         ; large format

gsn_panel(wks,plot,(/2,2/),resP)
print("convert to png file: "+graph1+".png")
system("convert -trim "+graph1+"dd.png"+" "+graph1+".png")
system("rm "+graph1+"dd.png")

delete(wks)
delete(plot)
delete(plot_tt)
print("--------------------------------------")
;***********************************************************************
print("plot "+graph3+".png")

plot = new(2,graphic)
wks = gsn_open_wks("png",graph3+"dd")
res1@cnMinLevelValF          = -0.2               ; min level
res1@cnMaxLevelValF          = 0.2           ; max level
res1@cnLevelSpacingF         = .05                  ; contour level spacing

plot(0)=gsn_csm_contour_map(wks, lem_tsi_sol_pr_rc,res1)
plot(1)=gsn_csm_contour_map(wks, lem_tsi_sol_tas_rc,res1)
;plot(0)=gsn_csm_contour_map(wks, rlag_pr(0,:,:),res1)
;plot(1)=gsn_csm_contour_map(wks, rlag(0,:,:),res1)

gsn_panel(wks,plot,(/1,2/),resP)
print("convert to png file: "+graph3+".png")
system("convert -trim "+graph3+"dd.png"+" "+graph3+".png")
system("rm "+graph3+"dd.png")

delete(wks)
delete(plot)

print("--------------------------------------")
print("plot "+graph4+".png")

res_ttest@cnMinLevelValF      = 0.00        ; set min contour level
res_ttest@cnMaxLevelValF      = 3.        ; set max contour level
res_ttest@cnLevelSpacingF     = 0.1        ; set contour spacing

plot = new(1,graphic)
wks = gsn_open_wks("png",graph4+"dd")
res1@cnMinLevelValF          = -0.2               ; min level
res1@cnMaxLevelValF          = 0.2           ; max level
res1@cnLevelSpacingF         = .05                  ; contour level spacing

res1@gsnLeftString                = "Unfiltered CESM-SOL  MCA-LIA"
res1@gsnRightString                = ""
plot(0)=gsn_csm_contour_map(wks, pr_MCALIA_unflt,res1)
;plot(1)=gsn_csm_contour_map(wks, tas_MCALIA,res1)
;printVarSummary(pr_ttest_unf)
;print(pr_ttest_unf({0},:))
plot444 = gsn_csm_contour(wks,pr_ttest_unf, res_ttest)
plot444   = ShadeLtContour(plot444, 0.1, 17)
overlay (plot(0), plot444)

gsn_panel(wks,plot,(/1,1/),resP)
print("convert to png file: "+graph4+".png")
system("convert -trim "+graph4+"dd.png"+" "+graph4+".png")
system("rm "+graph4+"dd.png")

delete(wks)
delete(plot)

print("--------------------------------------")

print("=================")
print("plot "+graph2+".png")
plot = new(2*nmodels,graphic)
plot_tt = new(2*nmodels,graphic)
wks = gsn_open_wks("png",graph2+"dd")
res1@cnMinLevelValF          = -1.              ; min level
res1@cnMaxLevelValF          = 1.           ; max level
res1@cnLevelSpacingF         = .2                  ; contour level spacing

do ii=0,nmodels-1

sol_rg_dd=regCoef_n(sol_tas_fldm_all(ii,:),sol_pr_all(ii,:,:,:),0,0)
co2_rg_dd=regCoef_n(co2_tas_fldm_all(ii,:),co2_pr_all(ii,:,:,:),0,0)

copy_VarCoords(dd2,sol_rg_dd(:,:))
copy_VarCoords(dd2,co2_rg_dd(:,:))

kk=2*ii
res1@gsnLeftString                = models(ii)+" SOL"
plot(kk)=gsn_csm_contour_map(wks, sol_rg_dd,res1)
dd1= sol_rg_dd
dd1=(/onedtond(abs(dd1@tval),(/dimsizes(dd1&lat),dimsizes(dd1&lon)/))/)
plot_tt(kk) = gsn_csm_contour(wks,dd1, res_ttest)
plot_tt(kk)   = ShadeGtContour(plot_tt(kk), 3., 17)
overlay (plot(kk), plot_tt(kk))
delete(dd1)

kk=2*ii+1
res1@gsnLeftString                = models(ii)+" CO2"
plot(kk)=gsn_csm_contour_map(wks, co2_rg_dd,res1)
dd1= co2_rg_dd
dd1=(/onedtond(abs(dd1@tval),(/dimsizes(dd1&lat),dimsizes(dd1&lon)/))/)
plot_tt(kk) = gsn_csm_contour(wks,dd1, res_ttest)
plot_tt(kk)   = ShadeGtContour(plot_tt(kk), 3., 17)
overlay (plot(kk), plot_tt(kk))
delete(dd1)
delete(sol_rg_dd)
delete(co2_rg_dd)
end do

gsn_panel(wks,plot,(/nmodels,2/),resP)
print("convert to png file: "+graph2+".png")
system("convert -trim "+graph2+"dd.png"+" "+graph2+".png")
system("rm "+graph2+"dd.png")


end