load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin

sfilter2=150
;models=(/"CAM4","CAM5","CanESM2","GISS","HadGEM2","HadGEM3","IPSL","MIROC","MPI","NorESM"/)
models=(/"CAM4","CAM5","CanESM2","HadGEM2","HadGEM3","IPSL","MPI"/)
nmodels=dimsizes(models)

graph1="/Users/stergios/analyses/19_pdrmip/figures/FF_supFig1"
graph2="/Users/stergios/analyses/19_pdrmip/figures/fig_TempPDRMIP_all"



sfile=True
cfile=True
nfile2=True
pconv=True

;solar
sfile@CAM4="ta_NCAR-CESM1-CAM4_SOL_minus_base_coupled_50yavg_zonal_mean.nc"
sfile@CAM5="ta_NCAR-CESM1-CAM5_SOL_minus_base_coupled_50yavg_zonal_mean.nc"

sfile@CanESM2="ta_CanESM2_SOL_minus_base_coupled_50yavg_zonal_mean.nc"
sfile@GISS="ta_GISS-E2-R_SOL_minus_base_coupled_50yavg_zonal_mean.nc"
sfile@HadGEM2="ta_HadGEM2_SOL_minus_base_coupled_50yavg_zonal_mean.nc"
sfile@HadGEM3="ta_HadGEM3_SOL_minus_base_coupled_50yavg_zonal_mean.nc"
sfile@IPSL="ta_IPSL-CM5A_SOL_minus_base_coupled_50yavg_zonal_mean.nc"
sfile@MIROC="ta_MIROC-SPRINTARS_SOL_minus_base_coupled_50yavg_zonal_mean.nc"
sfile@MPI="ta_MPI-ESM_SOL_minus_base_coupled_50yavg_zonal_mean.nc"
sfile@NorESM="ta_NorESM1_SOL_minus_base_coupled_50yavg_zonal_mean.nc"

;co2
cfile@CAM4="ta_NCAR-CESM1-CAM4_CO2_minus_base_coupled_50yavg_zonal_mean.nc"
cfile@CAM5="ta_NCAR-CESM1-CAM5_CO2_minus_base_coupled_50yavg_zonal_mean.nc"
cfile@CanESM2="ta_CanESM2_CO2_minus_base_coupled_50yavg_zonal_mean.nc"
cfile@GISS="ta_GISS-E2-R_CO2_minus_base_coupled_50yavg_zonal_mean.nc"
cfile@HadGEM2="ta_HadGEM2_CO2_minus_base_coupled_50yavg_zonal_mean.nc"
cfile@HadGEM3="ta_HadGEM3_CO2_minus_base_coupled_50yavg_zonal_mean.nc"
cfile@IPSL="ta_IPSL-CM5A_CO2_minus_base_coupled_50yavg_zonal_mean.nc"
cfile@MIROC="ta_MIROC-SPRINTARS_CO2_minus_base_coupled_50yavg_zonal_mean.nc"
cfile@MPI="ta_MPI-ESM_CO2_minus_base_coupled_50yavg_zonal_mean.nc"
cfile@NorESM="ta_NorESM1_CO2_minus_base_coupled_50yavg_zonal_mean.nc"

nfile2@CAM4="NCAR-CESM1-CAM4"
nfile2@CAM5="NCAR-CESM1-CAM5"
nfile2@CanESM2="CanESM2"
nfile2@GISS="GISS-E2-R"
nfile2@HadGEM2="HadGEM2"
nfile2@HadGEM3="HadGEM3"
nfile2@IPSL="IPSL-CM5A"
nfile2@MIROC="MIROC-SPRINTARS"
nfile2@MPI="MPI-ESM"
nfile2@NorESM="NorESM1"


Lsfile1="/Users/stergios/lpool/CESM-LEM/SOL/CESM-LEM_SOL_TS_850-2005_fltbwPAD"+sfilter2+"_fldm_ym.nc"
Lsfile2="/Users/stergios/lpool/CESM-LEM/SOL/CESM-LEM_SOL_T_850-2005_fltbwPAD"+sfilter2+"_zm_ym.nc"
Lsfile4="/Users/stergios/lpool/CESM-LEM/CESM-LEM_tsi_850-2005_fltbwPAD"+sfilter2+"_ym.nc"

Lcfile1="/Users/stergios/lpool/CESM-LEM/GHG/CESM-LEM_GHG_TS_1850-2005_fltbwPAD30_fldm_ym.nc"
Lcfile2="/Users/stergios/lpool/CESM-LEM/GHG/CESM-LEM_GHG_zmta_1850-2005_fltbwPAD30_zm_ym.nc"
;Lcfile2="/Users/stergios/lpool/CESM-LEM/GHG/CESM-LEM_GHG_zmta_1850-2005_zm_ym.nc"


pconv@CAM4=1.
pconv@CAM5=1.
pconv@CanESM2=1.
pconv@GISS=9999
pconv@HadGEM2=100.
pconv@HadGEM3=100.
pconv@IPSL=1.
pconv@MIROC=999.
pconv@MPI=1.
pconv@NorESM=999.


basedir="/Users/stergios/lpool/pdrmip/"

;svar=new((/nmodels,nlevs,nlats/),float)

newlats=ispan(-90,90,2)
newlevs=(/100000, 92500, 85000, 60000, 50000, 40000, 30000, 20000, 10000, 7000, 5000, 3000, 2000, 1000/)
nnewlats=dimsizes(newlats)
nnewlevs=dimsizes(newlevs)
svar=new((/nmodels,nnewlevs,nnewlats/),float)
cvar=svar
cvar2=svar
svar2=svar
;aggrement
svar_agr=svar
cvar_agr=svar

;++++++++++++++++++++++++++++++++++++++++
; READ PRDMIP
;++++++++++++++++++++++++++++++++++++++++

do imod=0,nmodels-1
model=models(imod)
;read solar
ifile=basedir+nfile2@$model$+"/"+sfile@$model$
print("read "+ifile)
ff1=addfile(ifile,"r")
dd1=ff1->ta

dd2=linint1_n_Wrap(dd1&$dd1!1$,dd1,False,newlats,0,1)
dd3=int2p_n_Wrap(pconv@$model$*dd2&$dd2!0$,dd2,newlevs,1,0)

 rad    = 4.0*atan(1.0)/180.0
clat   = cos(newlats*rad)
svar(imod,:,:)=(/dd3/avg(dd3({100000},:)*clat)/)
svar2(imod,:,:)=dd3
delete(dd1)
delete(dd2)
delete(dd3)
delete(ifile)
delete(ff1)

;read co2
ifile=basedir+nfile2@$model$+"/"+cfile@$model$
print("read "+ifile)
ff1=addfile(ifile,"r")
dd1=ff1->ta
dd2=linint1_n_Wrap(dd1&$dd1!1$,dd1,False,newlats,0,1)
dd3=int2p_n_Wrap(pconv@$model$*dd2&$dd2!0$,dd2,newlevs,1,0)
 rad    = 4.0*atan(1.0)/180.0
clat   = cos(newlats*rad)
cvar(imod,:,:)=(/dd3/avg(dd3({100000},:)*clat)/)
cvar2(imod,:,:)=dd3
delete(dd1)
delete(dd2)
delete(dd3)
delete(ifile)
delete(ff1)

end do
print("+++++++++++++++++++++++++++++++++++++++")

cvar2&LEV@units="hPa"
svar2&LEV@units="hPa"

svar!1="lev"

svar&lev=newlevs
svar&lev@units="Pa"
svar!2="lat"
svar&lat=newlats
svar&lat@units="degrees_north"

copy_VarCoords(svar(0,:,:),cvar(0,:,:))


;calculate ensemble mean 
svar_em=dim_avg_n_Wrap(svar,0)
cvar_em=dim_avg_n_Wrap(cvar,0)

; find agreement
do ii = 0,nmodels-1
  svar_agr(ii,:,:)=(/where((svar_em*svar(ii,:,:)).gt.0,1,0)/)
  cvar_agr(ii,:,:)=(/where((cvar_em*cvar(ii,:,:)).gt.0,1,0)/)
end do

sperm=100*(dim_sum_n(svar_agr,0)/nmodels)
cperm=100*(dim_sum_n(cvar_agr,0)/nmodels)
copy_VarCoords(svar_em,sperm)
copy_VarCoords(svar_em,cperm)

;++++++++++++++++++++++++++++++++++++++++++++
;read CESM LEM-SOL
;+++++++++++++++++++++++++++++++++++++++++++
print("read "+Lsfile1)
ff1=addfile(Lsfile1,"r")
Lstas=ff1->TS(:,0,0)
delete(ff1)

print("read "+Lsfile2)
ff1=addfile(Lsfile2,"r")
Lsta=ff1->T(:,:,:,0)
delete(ff1)


;read CESM LEM-GHG
print("read "+Lcfile1)
ff1=addfile(Lcfile1,"r")
Lctas=ff1->TS(:,0,0)
delete(ff1)

print("read "+Lcfile2)
ff1=addfile(Lcfile2,"r")
Lcta=ff1->zmta(:,:,:,0)
delete(ff1)

print("read "+Lsfile4)
ff1=addfile(Lsfile4,"r")
tsi=ff1->tsi(:,0,0)
delete(ff1)

;CESM calculate regressions
;SOLAR
flt_lem_sol_rg=regCoef_n(Lstas,Lsta,0,0)
flt2_lem_sol_rg=regCoef_n(tsi,Lsta,0,0)
copy_VarCoords(Lsta(0,:,:),flt_lem_sol_rg)
copy_VarCoords(Lsta(0,:,:),flt2_lem_sol_rg)
flt_lem_sol_rg&lev@units="hPa"
flt2_lem_sol_rg&lev@units="hPa"

;GHG
flt_lem_co2_rg=regCoef_n(Lctas,Lcta,0,0)
copy_VarCoords(Lcta(0,:,:),flt_lem_co2_rg)
flt_lem_co2_rg&lev@units="hPa"


;!!!!!!!!!!!!!!!!!!!
print("print "+graph1+".png")
wks = gsn_open_wks("png",graph1+"dd")
plot = new(4,graphic)
plot_tt = new(4,graphic)

res                      = True                  ; plot mods desired
res@gsnFrame          = False
res@gsnDraw          = False

res@cnFillOn             = True                  ; turn on color
res@cnLinesOn            = False                 ; no contour lines
res@cnLineLabelsOn       = False
res@lbLabelAutoStride    = True                  ; optimal labels
;res@lbLabelBarOn = False


res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF       = -2.
res@cnMaxLevelValF       =  2.
res@cnLevelSpacingF      = 0.5                  ; set contour spacing

res@tmYLMode                    = "Explicit"
res@tmYLLabels                  =(/"1000","","","","","500","","300","","100","","","","","50","","30","","10"/)
res@tmYLValues                  =(/1000,900,800,700,600,500,400,300,200,100,90,80,70,60,50,40,30,20,10/)


res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025
res@gsnLeftStringFontHeightF= 0.025
res@tiMainFontHeightF=0.02


res@tiXAxisFontHeightF=0.025
res@tiYAxisFontHeightF=0.025
res@cnInfoLabelOn  = False
res@tmYRMode                    = "Automatic"


res@gsnLeftString               = "SOLAR"
res@gsnRightString               = ""
res@tiYAxisString                  = "Pressure [hPa]"


;========================================================
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

res@cnFillPalette="BlueDarkRed18"

res@gsnLeftString               = "CESM-SOL"
plot(0) = gsn_csm_pres_hgt (wks,flt_lem_sol_rg,res)
df=1 ;?? correct
dd1= flt_lem_sol_rg
dd1=(/onedtond(abs(dd1@tval),(/dimsizes(dd1&lev),dimsizes(dd1&lat)/))/)
prob =dd1
nlats=dimsizes(prob&lat)
nlevs=dimsizes(prob&lev)


do ilev=0,nlevs-1
do  ilat= 0, nlats-1
  prob(ilev,ilat) = (/betainc( df/(df+dd1(ilev,ilat)^2), df/2.0, 0.5)/)
end do
end do
plot_tt(0) = gsn_csm_pres_hgt(wks,prob, res_ttest)
plot_tt(0)   = ShadeLtContour(plot_tt(0), 0.05, 17)
overlay (plot(0), plot_tt(0))
delete(dd1)
delete(prob)


res@gsnLeftString               = "CESM-GHG"

plot(1) = gsn_csm_pres_hgt (wks,flt_lem_co2_rg,res)
dd1= flt_lem_co2_rg
dd1=(/onedtond(abs(dd1@tval),(/dimsizes(dd1&lev),dimsizes(dd1&lat)/))/)
prob =dd1
nlats=dimsizes(prob&lat)
nlevs=dimsizes(prob&lev)

do ilev=0,nlevs-1
do  ilat= 0, nlats-1
  prob(ilev,ilat) = (/betainc( df/(df+dd1(ilev,ilat)^2), df/2.0, 0.5)/)
end do
end do

plot_tt(1) = gsn_csm_pres_hgt(wks,prob, res_ttest)
plot_tt(1)   = ShadeLtContour(plot_tt(1), 0.05, 17)
overlay (plot(1), plot_tt(1))
delete(dd1)
delete(prob)


res@gsnLeftString               = "CESM-SOL regressed to TSI"
res@cnMinLevelValF       = -.2
res@cnMaxLevelValF       =  .2
res@cnLevelSpacingF      = 0.05                  ; set contour spacing
plot(2) = gsn_csm_pres_hgt (wks,flt2_lem_sol_rg,res)

dd1= flt2_lem_sol_rg
dd1=(/onedtond(abs(dd1@tval),(/dimsizes(dd1&lev),dimsizes(dd1&lat)/))/)
prob =dd1
nlats=dimsizes(prob&lat)
nlevs=dimsizes(prob&lev)


do ilev=0,nlevs-1
do  ilat= 0, nlats-1
  prob(ilev,ilat) = (/betainc( df/(df+dd1(ilev,ilat)^2), df/2.0, 0.5)/)
end do
end do

plot_tt(2) = gsn_csm_pres_hgt(wks,prob, res_ttest)
plot_tt(2)   = ShadeLtContour(plot_tt(2), 0.05, 17)
overlay (plot(2), plot_tt(2))
delete(dd1)



res@cnMinLevelValF       = -5.
res@cnMaxLevelValF       =  5.
res@cnLevelSpacingF      = 1.                  ; set contour spacing

resP                     = True         ; modify the panel plot
resP@gsnMaximize         = True         ; large format
;resP@gsnPanelLabelBar = True
resP@lbLabelFontHeightF =.02
resP@gsnPanelFigureStrings= (/"a)","b)","c)","d)","e)","f)"/)
resP@gsnPanelFigureStringsFontHeightF=0.011
resP@amJust     = "TopLeft"
resP@lbLabelFontHeightF =.01

gsn_panel(wks,plot,(/2,2/),resP)
system("convert -trim "+graph1+"dd.png"+" "+graph1+".png")
system("rm "+graph1+"dd.png")

delete(plot)
delete(wks)
;++++++++++++++++++++++++++++++++++++++++++++
print("print "+graph2+".png")
wks = gsn_open_wks("png",graph2+"dd")
plot = new(2*nmodels,graphic)


do imod=0,nmodels-1
model=models(imod)
res@gsnLeftString               = model+" SOL"
plot(imod*2) = gsn_csm_pres_hgt (wks,svar2(imod,:,:),res)
res@gsnLeftString               = model+" CO2"
plot(imod*2+1) = gsn_csm_pres_hgt (wks,cvar2(imod,:,:),res)
end do

gsn_panel(wks,plot,(/nmodels,2/),resP)
system("convert -trim "+graph2+"dd.png"+" "+graph2+".png")
system("rm "+graph2+"dd.png")

delete(plot)
delete(wks)


end
