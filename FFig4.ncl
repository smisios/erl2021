load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin

lon1=30
lon2=300

;PDRMIP
sfilter=30
;models=(/"CanESM2","GISS-E2-R","HadGEM2","HadGEM3","IPSL-CM5A","MIROC-SPRINTARS","MPI-ESM","NCAR-CESM1-CAM5","NorESM1"/)

models=(/"CanESM2","GISS-E2-R","HadGEM2","HadGEM3","IPSL-CM5A","MIROC-SPRINTARS","MPI-ESM","NCAR-CESM1-CAM5","NorESM1"/)
nmodels=dimsizes(models)

Lsfile1="/Users/stergios/lpool/CESM-LEM/SOL/CESM-LEM_SOL_TS_850-2005_fltbwPAD150_fldm_ym.nc"
Lsfile2="/Users/stergios/lpool/CESM-LEM/SOL/CESM-LEM_SOL_TS_850-2005_fltbwPAD150_ym.nc"
Lsfile3="/Users/stergios/lpool/CESM-LEM/SOL/CESM-LEM_SOL_walker_850-2005_fltbwPAD150_merm_ym.nc"

Lcfile1="/Users/stergios/lpool/CESM-LEM/GHG/CESM-LEM_GHG_TS_1850-2005_fltbwPAD30_fldm_ym.nc"
Lcfile2="/Users/stergios/lpool/CESM-LEM/GHG/CESM-LEM_GHG_TS_1850-2005_fltbwPAD30_ym.nc"
Lcfile3="/Users/stergios/lpool/CESM-LEM/GHG/CESM-LEM_GHG_walker_1850-2005_fltbwPAD30_merm_ym.nc"

graph1="/Users/stergios/analyses/19_pdrmip/figures/FFig4"
graph2="/Users/stergios/analyses/19_pdrmip/figures/FFig4_all"


;basedir="/pool/lacie/pdrmip/"

 rad    = 4.0*atan(1.0)/180.0
; READ PDRMIP

;++++++++++++++++++++++++++++++
sol_tas_fldm_all=new((/nmodels,100/),double)
co2_tas_fldm_all=new((/nmodels,100/),double)

sol_ua_all=new((/nmodels,100,14,90/),float)
co2_ua_all=new((/nmodels,100,14,90/),float)

do ii=0,nmodels-1
model=models(ii)
 ;PDRMIP
sfile1="/Users/stergios/lpool/pdrmip/"+model+"/tas_"+model+"_solar_minus_base_coupled_100yr_fltbw"+sfilter+"_timeseries_fldm.nc"
cfile1="/Users/stergios/lpool/pdrmip/"+model+"/tas_"+model+"_co2_minus_base_coupled_100yr_fltbw"+sfilter+"_timeseries_fldm.nc"

if (model.eq."GISS-E2-R")

sfile2="/Users/stergios/lpool/pdrmip/"+model+"/ua_walker_"+model+"_pdrmip-solar_coupled_100yr_timeseries_ym.nc"
cfile2="/Users/stergios/lpool/pdrmip/"+model+"/ua_walker_"+model+"_pdrmip-co2x2_coupled_100yr_timeseries_ym.nc"
bfile2="/Users/stergios/lpool/pdrmip/"+model+"/ua_walker_"+model+"_pdrmip-base_coupled_100yr_timeseries_ym.nc"
else 
sfile2="/Users/stergios/lpool/pdrmip/"+model+"/ua_walker_"+model+"_SOL_minus_base_coupled_100yr_timeseries.nc"
cfile2="/Users/stergios/lpool/pdrmip/"+model+"/ua_walker_"+model+"_co2_minus_base_coupled_100yr_timeseries.nc"
bfile2="/Users/stergios/lpool/pdrmip/"+model+"/ua_walker_"+model+"_base_coupled_100yr_timeseries.nc"
end if
bfile1="/Users/stergios/lpool/pdrmip/"+model+"/tas_"+model+"_base_coupled_100yr_timeseries.nc"


 ;read solar
;==================
print("read "+sfile1)
ff1=addfile(sfile1,"r")
sol_tas_fldm=ff1->tas(:,0,0)
delete(ff1)

print("read "+sfile2)
ff1=addfile(sfile2,"r")
sol_ua=ff1->ua(:,:,0,:)
delete(ff1)

;read CO2
;===================

print("read "+cfile1)
ff1=addfile(cfile1,"r")
co2_tas_fldm=ff1->tas(:,0,0)
delete(ff1)

print("read "+cfile2)
ff1=addfile(cfile2,"r")
co2_ua=ff1->ua(:,:,0,:)
delete(ff1)


;read climatology
;===================
print("read "+bfile1)
ff1=addfile(bfile1,"r")
dd1=ff1->tas
clat   = cos(dd1&lat*rad)
base_tas_fldm=wgt_areaave_Wrap(dd1, clat, 1.0, 1)
base_tas_fldm_tm=dim_avg_n_Wrap(base_tas_fldm,0)
delete(clat)
delete(dd1)
delete(ff1)


print("read "+bfile2)
ff1=addfile(bfile2,"r")
base_ua=ff1->ua(:,:,0,:)
dd2=ff1->ua(0,:,0,:)
base_ua_tm=dim_avg_n_Wrap(base_ua,0)
delete(ff1)

sol_tas_fldm_all(ii,:)=sol_tas_fldm
co2_tas_fldm_all(ii,:)=co2_tas_fldm


do itt=0,99,1
sol_ua_all(ii,itt,:,:)=sol_ua(itt,:,:)
co2_ua_all(ii,itt,:,:)=co2_ua(itt,:,:)

end do
delete(sol_tas_fldm)
delete(co2_tas_fldm)
delete(sol_ua)
delete(co2_ua)
print("++++++++++++++++++++++++")
end do



;ensemble mean
sol_tas_fldm_em=dim_avg_n_Wrap(sol_tas_fldm_all,0)
sol_ua_em=dim_avg_n_Wrap(sol_ua_all,0)

co2_tas_fldm_em=dim_avg_n_Wrap(co2_tas_fldm_all,0)
co2_ua_em=dim_avg_n_Wrap(co2_ua_all,0)

;regress every model
sol_rg_em=regCoef_n(sol_tas_fldm_em,sol_ua_em,0,0)
co2_rg_em=regCoef_n(co2_tas_fldm_em,co2_ua_em,0,0)

copy_VarCoords(dd2,sol_rg_em)
copy_VarCoords(dd2,co2_rg_em)


; read LEM
;read model
;===========
print("read "+Lsfile1)
ff1=addfile(Lsfile1,"r")
lem_sol_tas_fldm=ff1->TS(:,0,0)
delete(ff1)

print("read "+Lsfile3)
ff1=addfile(Lsfile3,"r")
lem_sol_ua=ff1->ua
delete(ff1)
lem_sol_ua&lev@units="hPa"

print("read "+Lcfile1)
ff1=addfile(Lcfile1,"r")
lem_co2_tas_fldm=ff1->TS(:,0,0)
delete(ff1)

print("read "+Lcfile3)
ff1=addfile(Lcfile3,"r")
lem_co2_ua=ff1->ua
delete(ff1)
lem_co2_ua&lev@units="hPa"

lem_sol_ua_rc           = regCoef_n(lem_sol_tas_fldm ,lem_sol_ua,0,0)
copy_VarCoords(lem_sol_ua(0,:,:,:),lem_sol_ua_rc )
printVarSummary(lem_sol_ua)
;efdof=equiv_sample_size(lem_sol_ua(lev|:,lat|:,lon|:,time|:),0.05,0)
;printVarSummary(sol_ua_em)
;sol_ua_em!0="time"
;efdof2=equiv_sample_size(sol_ua_em(lev|:,lon|:,time|:),0.05,0)
;print(max(efdof))
;print(efdof2)
;exit


lem_co2_ua_rc           = regCoef_n(lem_co2_tas_fldm ,lem_co2_ua,0,0)
copy_VarCoords(lem_co2_ua(0,:,:,:),lem_co2_ua_rc )


print("plot "+graph1+".png")
wks = gsn_open_wks("png",graph1+"dd")
plot = new(4,graphic)
plot_tt = new(4,graphic) 

res                              = True
res@gsnFrame                    = False        ; don't draw yet
res@gsnDraw                    = False        ; don't draw yet


res@cnFillOn                = True             ; turn on color fill
res@cnLevelSelectionMode         = "ManualLevels"        ; manual contour levels

res@cnLinesOn                    = False
res@cnLineLabelsOn               = False                  ; turn on line labels
res@trYLog                       = False
res@cnInfoLabelOn                = False
res@lbLabelFontHeightF=0.013
res@tmYLLabelFontHeightF = 0.017	
res@trYMinF=90
;-------
; labels
;-------

res@gsnRightString               = " " 
res@tiYAxisString                  = "Pressure (hPa)"
res@tiMainString                 = ""
res@tiXAxisString                  = ""
res@cnLevelSpacingF              = 0.2           ; contour interval
res@cnMinLevelValF               = -1.                  ; min level
res@cnMaxLevelValF               = 1.                  ; max level

;----
;axis
;----
res@tmXTMinorOn                  = False
res@tmYLMode                    = "Explicit"
res@tmYLLabels                  = (/"1000","800","600","400","200","100","50"/)
res@tmYLValues                  = (/1000,800,600,400,200,100,50/)
res@tmYRMode                    = "Automatic"          ; turn off special labels on right axis


;print atmosphere
;-----------------

res@cnFillPalette       = "nrl_sirkes"
res@cnFillOn            = True             ; color fill plot
res@cnLinesOn           = False
res@cnLineLabelsOn      = False
res@cnInfoLabelOn       = False
res@cnSpanFillPalette=True
res@gsnLeftStringFontHeightF= 0.035
res@gsnRightStringFontHeightF= 0.025
res@tmXBLabelFontHeightF = 0.035
res@tmYLLabelFontHeightF = 0.035

res@tiMainFontHeightF=0.035


;==================================================================================
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



res@gsnLeftString               = "PDRMIP-SOL"
plot(0)=gsn_csm_pres_hgt(wks,sol_rg_em(:,{lon1:lon2}),res)
dd1= sol_rg_em
dd1=(/onedtond(abs(dd1@tval),(/dimsizes(dd1&lev),dimsizes(dd1&lon)/))/)

prob =dd1
nlons=dimsizes(prob&lon)
nlevs=dimsizes(prob&lev)

df=4 ;?? correct
do ilev=0,nlevs-1
do  ilon= 0, nlons-1
  prob(ilev,ilon) = (/betainc( df/(df+dd1(ilev,ilon)^2), df/2.0, 0.5)/)
end do
end do

plot_tt(0) = gsn_csm_pres_hgt(wks,prob, res_ttest)
plot_tt(0)   = ShadeLtContour(plot_tt(0), .05, 17)
overlay (plot(0), plot_tt(0))
delete(dd1)
delete(prob)


;plot_tt(0) = gsn_csm_pres_hgt(wks,dd1, res_ttest)
;plot_tt(0)   = ShadeGtContour(plot_tt(0), 2., 17)
;overlay (plot(0), plot_tt(0))
;delete(dd1)




res@gsnLeftString               = "PDRMIP-CO2"
plot(1)=gsn_csm_pres_hgt(wks,co2_rg_em(:,{lon1:lon2}),res)
dd1= co2_rg_em
dd1=(/onedtond(abs(dd1@tval),(/dimsizes(dd1&lev),dimsizes(dd1&lon)/))/)
prob =dd1
nlons=dimsizes(prob&lon)
nlevs=dimsizes(prob&lev)

do ilev=0,nlevs-1
do  ilon= 0, nlons-1
  prob(ilev,ilon) = (/betainc( df/(df+dd1(ilev,ilon)^2), df/2.0, 0.5)/)
end do
end do
plot_tt(1) = gsn_csm_pres_hgt(wks,prob, res_ttest)
plot_tt(1)   = ShadeLtContour(plot_tt(1), .05, 17)
overlay (plot(1), plot_tt(1))
delete(dd1)
delete(prob)

;plot_tt(1) = gsn_csm_pres_hgt(wks,dd1, res_ttest)
;plot_tt(1)   = ShadeGtContour(plot_tt(1), 2., 17)
;overlay (plot(1), plot_tt(1))
;delete(dd1)



res@gsnLeftString               = "CESM-SOL"
plot(2)=gsn_csm_pres_hgt(wks,lem_sol_ua_rc(:,0,{lon1:lon2}),res)
dd1= lem_sol_ua_rc(:,0,:)
dd1=(/onedtond(abs(dd1@tval),(/dimsizes(dd1&lev),dimsizes(dd1&lon)/))/)
prob =dd1
nlons=dimsizes(prob&lon)
nlevs=dimsizes(prob&lev)

do ilev=0,nlevs-1
do  ilon= 0, nlons-1
  prob(ilev,ilon) = (/betainc( df/(df+dd1(ilev,ilon)^2), df/2.0, 0.5)/)
end do
end do
plot_tt(2) = gsn_csm_pres_hgt(wks,prob, res_ttest)
plot_tt(2)   = ShadeLtContour(plot_tt(2), .05, 17)
overlay (plot(2), plot_tt(2))
delete(dd1)
delete(prob)

;plot_tt(2) = gsn_csm_pres_hgt(wks,dd1, res_ttest)
;plot_tt(2)   = ShadeGtContour(plot_tt(2), 2., 17)
;overlay (plot(2), plot_tt(2))
;delete(dd1)

res@gsnLeftString               = "CESM-GHG"
plot(3)=gsn_csm_pres_hgt(wks,lem_co2_ua_rc(:,0,{lon1:lon2}),res)
dd1= lem_co2_ua_rc(:,0,:)
dd1=(/onedtond(abs(dd1@tval),(/dimsizes(dd1&lev),dimsizes(dd1&lon)/))/)
prob =dd1
nlons=dimsizes(prob&lon)
nlevs=dimsizes(prob&lev)

do ilev=0,nlevs-1
do  ilon= 0, nlons-1
  prob(ilev,ilon) = (/betainc( df/(df+dd1(ilev,ilon)^2), df/2.0, 0.5)/)
end do
end do
plot_tt(3) = gsn_csm_pres_hgt(wks,prob, res_ttest)
plot_tt(3)   = ShadeLtContour(plot_tt(3), .05, 17)
overlay (plot(3), plot_tt(3))
delete(dd1)
delete(prob)


;plot_tt(3) = gsn_csm_pres_hgt(wks,dd1, res_ttest)
;plot_tt(3)   = ShadeGtContour(plot_tt(3), 2., 17)
;overlay (plot(3), plot_tt(3))
;delete(dd1)


resP                     = True         ; modify the panel plot
resP@gsnMaximize         = True         ; large format

gsn_panel(wks,plot,(/2,2/),resP)
system("convert -trim "+graph1+"dd.png"+" "+graph1+".png")
system("rm "+graph1+"dd.png")
;==================================================================================
delete(wks)
delete(plot)
delete(plot_tt)
print("=================")

res_ttest@cnMinLevelValF      = 0.00        ; set min contour level
res_ttest@cnMaxLevelValF      = 3.        ; set max contour level
res_ttest@cnLevelSpacingF     = 0.1        ; set contour spacing

print("plot "+graph2+".png")
plot = new(2*nmodels,graphic)
plot_tt = new(2*nmodels,graphic)
wks = gsn_open_wks("png",graph2+"dd")

iii=0
do imod=0,nmodels-1

sol_rg_dd=regCoef_n(sol_tas_fldm_all(imod,:),sol_ua_all(imod,:,:,:),0,0)
co2_rg_dd=regCoef_n(co2_tas_fldm_all(imod,:),co2_ua_all(imod,:,:,:),0,0)

copy_VarCoords(dd2,sol_rg_dd(:,:))
copy_VarCoords(dd2,co2_rg_dd(:,:))



res@gsnLeftString               = models(imod)+" SOL"
plot(iii)=gsn_csm_pres_hgt(wks,sol_rg_dd(:,{lon1:lon2}),res)
dd1= sol_rg_dd
dd1=(/onedtond(abs(dd1@tval),(/dimsizes(dd1&lev),dimsizes(dd1&lon)/))/)

plot_tt(iii) = gsn_csm_pres_hgt(wks,dd1, res_ttest)
plot_tt(iii)   = ShadeGtContour(plot_tt(iii), 2., 17)
overlay (plot(iii), plot_tt(iii))
delete(dd1)


iii=iii+1
res@gsnLeftString               = models(imod)+" CO2"
plot(iii)=gsn_csm_pres_hgt(wks,co2_rg_dd(:,{lon1:lon2}),res)
dd1= co2_rg_dd
dd1=(/onedtond(abs(dd1@tval),(/dimsizes(dd1&lev),dimsizes(dd1&lon)/))/)

plot_tt(iii) = gsn_csm_pres_hgt(wks,dd1, res_ttest)
plot_tt(iii)   = ShadeGtContour(plot_tt(iii), 2., 17)
overlay (plot(iii), plot_tt(iii))
delete(dd1)


iii=iii+1
delete(sol_rg_dd)
delete(co2_rg_dd)
end do

resP                     = True         ; modify the panel plot
resP@gsnMaximize         = True         ; large format

gsn_panel(wks,plot,(/nmodels,2/),resP)
system("convert -trim "+graph2+"dd.png"+" "+graph2+".png")
system("rm "+graph2+"dd.png")




end
