load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin

lemflt_sol=150
lemflt_co2=30
sfilter=30
bwd_sol=75
bwd_co2=15
selyear1=900
selyear2=2005

;PDRMIP
models=(/"CanESM2","GISS-E2-R","HadGEM2","HadGEM3","IPSL-CM5A","MIROC-SPRINTARS","MPI-ESM","NCAR-CESM1-CAM5","NorESM1"/)

nmodels=dimsizes(models)

bfile1="/Users/stergios/lpool/pdrmip/solar_paper/precipitation/tas_ensemble_mean_base_coupled_100yr_timeseries_fldm_tm.nc"
bfile2="/Users/stergios/lpool/pdrmip/solar_paper/precipitation/pr_ensemble_mean_base_coupled_100yr_timeseries_fldm.nc"

;CESM filtered with a butterworth filter
Lsfile1="/Users/stergios/lpool//CESM-LEM/SOL/CESM-LEM_SOL_TS_850-2005_fltbwPAD"+lemflt_sol+"_fldm_ym.nc"
Lsfile2="/Users/stergios/lpool/CESM-LEM/SOL/CESM-LEM_SOL_PREC_850-2005_fltbwPAD"+lemflt_sol+"_fldm_ym.nc"

Lcfile1="/Users/stergios/lpool/CESM-LEM/GHG/CESM-LEM_GHG_TS_1850-2005_fltbwPAD"+lemflt_co2+"_fldm_ym.nc"
Lcfile2="/Users/stergios/lpool/CESM-LEM/GHG/CESM-LEM_GHG_PREC_1850-2005_fltbwPAD"+lemflt_co2+"_fldm_ym.nc"

Lbfile1="/Users/stergios/lpool/CESM-LEM/GHG/CESM-LEM_GHG_RG4_TS_900-1200_fldm_tm.nc"
Lbfile2="/Users/stergios/lpool/CESM-LEM/GHG/CESM-LEM_GHG_PREC_900-1200_fldm_tm.nc"

graph1="/Users/stergios/analyses/19_pdrmip/figures/FFig2"

rad    = 4.0*atan(1.0)/180.0

; READ PDRMIP
;++++++++++++++++++++++++++++++
sol_tas_all=new((/nmodels,100/),double)
sol_pr_all=new((/nmodels,100/),double)
co2_tas_all=new((/nmodels,100/),double)
co2_pr_all=new((/nmodels,100/),double)

;READ PDRMIP
;------------
do ii=0,nmodels-1
 model=models(ii)
 ;surface temperature
 sfile1="/Users/stergios/lpool/pdrmip/"+model+"/tas_"+model+"_solar_minus_base_coupled_100yr_fltbw"+sfilter+"_timeseries_fldm.nc"
 cfile1="/Users/stergios/lpool/pdrmip/"+model+"/tas_"+model+"_co2_minus_base_coupled_100yr_fltbw"+sfilter+"_timeseries_fldm.nc"
 ;precipition
 sfile2="/Users/stergios/lpool/pdrmip/"+model+"/pr_"+model+"_SOL_minus_base_coupled_100yr_fltbw"+sfilter+"_timeseries_fldm.nc"
 cfile2="/Users/stergios/lpool/pdrmip/"+model+"/pr_"+model+"_CO2_minus_base_coupled_100yr_fltbw"+sfilter+"_timeseries_fldm.nc"
 ; climatology
 bfile1="/Users/stergios/lpool/pdrmip/"+model+"/tas_"+model+"_base_coupled_100yr_timeseries.nc"
 bfile2="/Users/stergios/lpool/pdrmip/"+model+"/pr_"+model+"_base_coupled_100yr_timeseries.nc"

 ;read solar
;==================
print("read "+sfile1)
ff1=addfile(sfile1,"r")
sol_tas=ff1->tas(:,0,0)
delete(ff1)
print("read "+sfile2)
ff1=addfile(sfile2,"r")
sol_pr=ff1->pr(:,0,0)
delete(ff1)

;read CO2
;===================
print("read "+cfile1)
ff1=addfile(cfile1,"r")
co2_tas=ff1->tas(:,0,0)
delete(ff1)
print("read "+cfile2)
ff1=addfile(cfile2,"r")
co2_pr=ff1->pr(:,0,0)
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
dd1=ff1->pr
clat   = cos(dd1&lat*rad)
base_pr=wgt_areaave_Wrap(dd1, clat, 1.0, 1)
base_pr_tm=dim_avg_n_Wrap(base_pr,0)
delete(clat)
delete(dd1)
delete(ff1)

sol_tas_all(ii,:)=sol_tas
co2_tas_all(ii,:)=co2_tas
; make percentages to climatology
sol_pr_all(ii,:)=100*(sol_pr)/base_pr_tm
co2_pr_all(ii,:)=100*(co2_pr)/base_pr_tm

delete(sol_tas)
delete(co2_tas)
delete(sol_pr)
delete(co2_pr)
print("++++++++++++++++++++++++")
end do

;ensemble mean
sol_tas_em=dim_avg_n_Wrap(sol_tas_all,0)
sol_pr_em=dim_avg_n_Wrap(sol_pr_all,0)

co2_tas_em=dim_avg_n_Wrap(co2_tas_all,0)
co2_pr_em=dim_avg_n_Wrap(co2_pr_all,0)

;make regressions
sol_rg=regline(sol_tas_em,sol_pr_em)
co2_rg=regline(co2_tas_em,co2_pr_em)

;regress every model
sol_rg_all=regCoef_n(sol_tas_all,sol_pr_all,1,1)
co2_rg_all=regCoef_n(co2_tas_all,co2_pr_all,1,1)

;print pdrmip
print ("++++++++++++++++++++++++++++++")
print ("PDRMIP")
print ("PDRMIP-SOL hydrological sens.: "+sol_rg)
print ("PDRMIP-CO2 hydrological sens.: "+co2_rg)
print ("PDRMIP-SOL intercept: "+sol_rg@yintercept)
print ("PDRMIP-CO2 intercept: "+co2_rg@yintercept)

print ("PDRMIP-SOL sens. MAX:"+max(sol_rg_all)+" min: "+min(sol_rg_all))
print ("PDRMIP-CO2 sens. MAX:"+max(co2_rg_all)+" min: "+min(co2_rg_all))
print ("++++++++++++++++++++++++++++++")

print (co2_rg_all)

;8888888888888888888888888888888888888888888888888888888888888
;read CESM
;8888888888888888888888888888888888888888888888888888888888888
;read SOL
;===================
print("read "+Lsfile1)
ff1=addfile(Lsfile1,"r")
lem_sol_tas=ff1->TS(:,0,0)
lntimes=dimsizes(lem_sol_tas)
delete(ff1)

print("read "+Lsfile2)
ff1=addfile(Lsfile2,"r")
lem_sol_pr=ff1->PREC(:,0,0)
delete(ff1)

;read CO2
;===================
print("read "+Lcfile1)
ff1=addfile(Lcfile1,"r")
lem_co2_tas=ff1->TS(:,0,0)
lntimes_co2=dimsizes(lem_co2_tas)
delete(ff1)
print("read "+Lcfile2)
ff1=addfile(Lcfile2,"r")
lem_co2_pr=ff1->PREC(:,0,0)
delete(ff1)

;read climatology
;===================
print("read "+Lbfile1)
ff1=addfile(Lbfile1,"r")
lem_tas_clim=ff1->TS(0,0,0)
delete(ff1)

print("read "+Lbfile2)
ff1=addfile(Lbfile2,"r")
lem_pr_clim=ff1->PRECC(0,0,0)
delete(ff1)

; scale
lem_co2_pr_scaled1=(/100*(lem_co2_pr-lem_pr_clim)/lem_pr_clim/)
lem_sol_pr_scaled1=(/100*(lem_sol_pr-lem_pr_clim)/lem_pr_clim/)

lem_sol_tas_scaled1=(/lem_sol_tas-lem_tas_clim/)
lem_co2_tas_scaled1=(/lem_co2_tas-lem_tas_clim/)

lem_sol_rg=regline(lem_sol_tas_scaled1,lem_sol_pr_scaled1)
lem_co2_rg=regline(lem_co2_tas_scaled1,lem_co2_pr_scaled1)

print(lem_sol_rg)
;print CESM
print ("++++++++++++++++++++++++++++++")
print ("CESM")
print ("CESM-SOL hydrological sens.: "+lem_sol_rg)
print ("CESM-GHG hydrological sens.: "+lem_co2_rg)
print ("CESM-SOL intercept: "+lem_sol_rg@yintercept)
print ("CESM-CO2 intercept: "+lem_co2_rg@yintercept)
print ("++++++++++++++++++++++++++++++")


;estrapolate CESM solid lines
xx_sol=fspan(-1,max(sol_tas_em)+1,10)
lem_sol_pr_extrap = lem_sol_rg*xx_sol + lem_sol_rg@yintercept              ; use solid line
xx_co2=fspan(-1,max(co2_tas_em)+1,10)
lem_co2_pr_extrap = lem_co2_rg*xx_sol + lem_co2_rg@yintercept  

;estrapolate individual PRMIP
sol_pr_extrap_all=new((/nmodels,dimsizes(xx_sol)/),double)
co2_pr_extrap_all=new((/nmodels,dimsizes(xx_co2)/),double)
sol_tas_extrap_all=new((/nmodels,dimsizes(xx_sol)/),double)
co2_tas_extrap_all=new((/nmodels,dimsizes(xx_co2)/),double)

do ii=0,nmodels-1
dd1=sol_rg_all@yintercept
sol_pr_extrap_all(ii,:)=sol_rg_all(ii)*xx_sol + dd1(ii)
sol_tas_extrap_all(ii,:)=xx_sol
delete(dd1)
dd1=co2_rg_all@yintercept
co2_pr_extrap_all(ii,:)=co2_rg_all(ii)*xx_co2 + dd1(ii)
co2_tas_extrap_all(ii,:)=xx_co2
delete(dd1)
end do

; add standard deviation
 sxp    = new( (/2*dimsizes(xx_sol)/), double )
 syp    = new( (/2*dimsizes(xx_sol)/), double )
  cxp    = new( (/2*dimsizes(xx_co2)/), double )
 cyp    = new( (/2*dimsizes(xx_co2)/), double )

 do k=0,dimsizes(xx_sol)-1
   dx             = 2*lem_sol_rg@rstd*xx_sol(k)
   syp(k)          = lem_sol_pr_extrap(k) + dx
   sxp(k)          = xx_sol(k)
   sxp(2*dimsizes(xx_sol)-1-k) = xx_sol(k)
   syp(2*dimsizes(xx_sol)-1-k) = lem_sol_pr_extrap(k) - dx
 end do

 
  do k=0,dimsizes(xx_co2)-1
   dx             = 2*lem_co2_rg@rstd*xx_co2(k)
   cyp(k)          = lem_co2_pr_extrap(k) + dx
   cxp(k)          = xx_co2(k)
   cxp(2*dimsizes(xx_co2)-1-k) = xx_co2(k)
   cyp(2*dimsizes(xx_co2)-1-k) = lem_co2_pr_extrap(k) - dx
 end do
 
print("=================")
print("plot "+graph1+".png")
wks = gsn_open_wks("png",graph1+"dd")

 res                     = True                   ; plot mods desired
 res@gsnDraw            =  False                   ; do not draw the plot
 res@gsnFrame           =  False                   ; do not advance the frame
 res@xyMarkLineModes     = (/"Markers","Markers"/)
 res@xyMarkers           = (/4,6/)
 res@xyMarkerColors       = (/"black","red"/)                  ; Marker color
 res@xyMarkerSizeF       = 0.02                  ; Marker size (default 0.01)
 res@xyMarkerThicknessF  = 5.
 res@trXMinF=-0.1
res@trXMaxF=3.3
 res@trYMinF=-3
res@trYMaxF=7
 
  res@tiYAxisString = "Precipitation Anomaly [%]"
 res@tiXAxisString = " Temperature Anomaly [K]"
 
 
 plot  = gsn_csm_xy (wks,sol_tas_em,sol_pr_em,res)
 res@xyMarkerColors       = (/"red","red"/)                  ; Marker color
 plot1  = gsn_csm_xy (wks,co2_tas_em,co2_pr_em,res)
 overlay (plot, plot1)
  res@xyMarkerColors       = (/"blue","red"/)                  ; Marker color
 res@xyMarkerSizeF       = 0.015                 ; Marker size (default 0.01)
  res@xyMarkers           = (/4,6/)
  plot2  = gsn_csm_xy (wks,lem_co2_tas_scaled1,lem_co2_pr_scaled1,res)
 overlay (plot, plot2)
 
   res@xyMarkers           = (/4,6/)
   res@xyMarkerColors       = (/"orange","red"/) 

  plot3  = gsn_csm_xy (wks,lem_sol_tas_scaled1,lem_sol_pr_scaled1,res)
 overlay (plot, plot3)
 
 res@xyMarkLineModes=(/"Lines","Lines"/)
   res@xyLineColor     = "red"
   res@xyLineThicknessF= 2.
   res@xyDashPattern=1
   plot6  = gsn_csm_xy (wks,co2_tas_extrap_all,co2_pr_extrap_all,res)
   overlay (plot, plot6)
   
   
   res@xyLineColor     = "black"
   plot7  = gsn_csm_xy (wks,sol_tas_extrap_all,sol_pr_extrap_all,res)
   overlay (plot, plot7)
   
   
  res@xyLineColors     = (/"orange","blue"/)
 ;res@xyLineThicknesses = (/7.0,7.0/)
  res@xyLineThicknessF= 9.
   res@xyDashPattern=0
   plot4  = gsn_csm_xy (wks,xx_sol,lem_sol_pr_extrap,res)
 overlay (plot, plot4)
 res@xyLineColors     = (/"blue","red"/)
    plot5  = gsn_csm_xy (wks,xx_co2,lem_co2_pr_extrap,res)
 overlay (plot, plot5)
 
 ; gsres                   = True                        ; poly res
 ;gsres@tfPolyDrawOrder   = "Predraw"                   ; draw this first
 ;gsres@gsFillColor       = "Grey"                 ; color chosen

 ;dummy = gsn_add_polygon (wks,plot,sxp,syp,gsres)
 ; dummy2 = gsn_add_polygon (wks,plot,cxp,cyp,gsres)
 
 
 ;add text
 txres               = True                            
 txres@txFontHeightF = 0.015             ; Set the font height
 txres@txFontColor="Black"
 text = gsn_add_text(wks,plot,"PDRMIP-SOL: "+sprintf("%6.2g",sol_rg)+" %/K",2.5,0.4 ,txres)
 txres@txFontColor="Red"
 text = gsn_add_text(wks,plot,"PDRMIP-CO2: "+sprintf("%6.2g",co2_rg)+" %/K",2.5,0.1 ,txres)
 txres@txFontColor="Orange"
 text = gsn_add_text(wks,plot,"CESM-SOL: "+sprintf("%6.2g",lem_sol_rg)+" %/K",2.55,-.2 ,txres)
  txres@txFontColor="Blue"
 text = gsn_add_text(wks,plot,"CESM-GHG: "+sprintf("%6.2g",lem_co2_rg)+" %/K",2.53,-0.5 ,txres)
 
resP                     = True         ; modify the panel plot
resP@gsnMaximize         = True         ; large format
;resP@gsnPanelFigureStrings= (/"a)","b)","c)","d)","e)","f)"/)
;resP@gsnPanelFigureStringsFontHeightF=0.011
;resP@amJust     = "TopLeft"
;resP@lbLabelFontHeightF =.01

gsn_panel(wks,plot,(/1,1/),resP)
system("convert -trim "+graph1+"dd.png"+" "+graph1+".png")
system("rm "+graph1+"dd.png")

delete(plot)
delete(wks)

end