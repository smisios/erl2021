load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin

;LEM

;PDRMIP
sfilter=30
sfilter2=150
models=(/"CanESM2","GISS-E2-R","HadGEM2","HadGEM3","IPSL-CM5A","MIROC-SPRINTARS","MPI-ESM","NCAR-CESM1-CAM5","NorESM1"/)
nmodels=dimsizes(models)

bfile2="/Users/stergios/lpool/pdrmip/solar_paper/precipitation/pr_ensemble_mean_base_coupled_100yr_timeseries_fldm_tm.nc"

Lsfile1="/Users/stergios/lpool/CESM-LEM/SOL/CESM-LEM_SOL_TS_850-2005_fltbwPAD"+sfilter2+"_fldm_ym.nc"
Lsfile11="/Users/stergios/lpool/CESM-LEM/SOL/CESM-LEM_SOL_TS_850-2005_fldm_ym.nc"
Lsfile2="/Users/stergios/lpool/CESM-LEM/SOL/CESM-LEM_SOL_PREC_850-2005_fltbwPAD"+sfilter2+"_fldm_ym.nc"
Lsfile3="/Users/stergios/lpool/CESM-LEM/CESM-LEM_tsi_850-2005_fltbwPAD"+sfilter2+"_ym.nc"

Lcfile1="/Users/stergios/lpool/CESM-LEM/GHG/CESM-LEM_GHG_TS_850-2005_fltbwPAD30_fldm_ym.nc"
Lcfile11="/Users/stergios/lpool/CESM-LEM/GHG/CESM-LEM_GHG_TS_850-2005_fldm_ym.nc"
Lcfile2="/Users/stergios/lpool/CESM-LEM/GHG/CESM-LEM_GHG_PREC_850-2005_fltbwPAD30_fldm_ym.nc"
Lcfile3="/Users/stergios/lpool/CESM-LEM/CESM-LEM_CO2_850-2007_fltbwPAD30_ym.nc"

Lbfile1="/Users/stergios/lpool/CESM-LEM/GHG/CESM-LEM_GHG_RG4_TS_850-1200_fldm_tm.nc"
Lbfile2="/Users/stergios/lpool/CESM-LEM/GHG/CESM-LEM_GHG_PREC_850-1200_fldm_tm.nc"

graph1="/Users/stergios/analyses/19_pdrmip/figures/FFig1"
;graph1="/Users/stergios/analyses/19_pdrmip/figures/FFig1"


;read TSI
;=========
print("read "+Lsfile3)
ff1=addfile(Lsfile3,"r")
tsi=ff1->tsi(:,0,0)
printVarSummary(tsi)
tsi_orig=(/tsi-tsi(0)/)
tsi=(/tsi-avg(tsi)/)
delete(ff1)

;read TSI
;=========
print("read "+Lcfile3)
ff1=addfile(Lcfile3,"r")
co2=ff1->CO2(0:1155)
co2_orig=(/co2-co2(0)/)
co2=(/co2-avg(co2)/)
delete(ff1)

;++++++++++++++++++++++++++++++
sol_tas_fldm_all=new((/nmodels,100/),double)
co2_tas_fldm_all=new((/nmodels,100/),double)
sol_pr_fldm_all=new((/nmodels,100/),double)
co2_pr_fldm_all=new((/nmodels,100/),double)

;PDRMIP
;===================
do ii=0,nmodels-1
model=models(ii)
 ;PDRMIP
sfile1="/Users/stergios/lpool/pdrmip/"+model+"/tas_"+model+"_solar_minus_base_coupled_100yr_fltbw"+sfilter+"_timeseries_fldm.nc"
cfile1="/Users/stergios/lpool/pdrmip/"+model+"/tas_"+model+"_co2_minus_base_coupled_100yr_fltbw"+sfilter+"_timeseries_fldm.nc"

sfile2="/Users/stergios/lpool/pdrmip/"+model+"/pr_"+model+"_SOL_minus_base_coupled_100yr_fltbw"+sfilter+"_timeseries_fldm.nc"
cfile2="/Users/stergios/lpool/pdrmip/"+model+"/pr_"+model+"_CO2_minus_base_coupled_100yr_fltbw"+sfilter+"_timeseries_fldm.nc"

 ;read solar
;==================
print("read "+sfile1)
ff1=addfile(sfile1,"r")
sol_tas_fldm_all(ii,:)=ff1->tas(:,0,0)
delete(ff1)
print("read "+sfile2)
ff1=addfile(sfile2,"r")
sol_pr_fldm_all(ii,:)=ff1->pr(:,0,0)
delete(ff1)

;read CO2
;===================
print("read "+cfile1)
ff1=addfile(cfile1,"r")
co2_tas_fldm_all(ii,:)=ff1->tas(:,0,0)
delete(ff1)
print("read "+cfile2)
ff1=addfile(cfile2,"r")
co2_pr_fldm_all(ii,:)=ff1->pr(:,0,0)
delete(ff1)
end do

sol_tas=dim_avg_n_Wrap(sol_tas_fldm_all,0)
co2_tas=dim_avg_n_Wrap(co2_tas_fldm_all,0)
sol_pr=dim_avg_n_Wrap(sol_pr_fldm_all,0)
co2_pr=dim_avg_n_Wrap(co2_pr_fldm_all,0)

time_pdrmip=ispan(1,dimsizes(sol_tas),1)

;read climatology
;===================
print("read "+bfile2)
ff1=addfile(bfile2,"r")
pr_clim=ff1->pr(0,0,0)
delete(ff1)
;**************************************

; READ LME
;read SOL
;===================
print("read "+Lsfile1)
ff1=addfile(Lsfile1,"r")
lem_sol_tas=ff1->TS(:,0,0)
sol_nsteps=dimsizes(lem_sol_tas)
time_sol=ispan(850,2005,1)
delete(ff1)

;read SOL unfiletred
;===================
print("read "+Lsfile11)
ff1=addfile(Lsfile11,"r")
lem_sol_tas_unf=ff1->TS(:,0,0)
sol_nsteps=dimsizes(lem_sol_tas_unf)
time_sol=ispan(850,2005,1)
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
time_co2=ispan(850,2005,1)
delete(ff1)

;read CO2
;===================
print("read "+Lcfile11)
ff1=addfile(Lcfile11,"r")
lem_co2_tas_unf=ff1->TS(:,0,0)
time_co2_unf=ispan(850,2005,1)
delete(ff1)

print("read "+Lcfile2)
ff1=addfile(Lcfile2,"r")
lem_co2_pr=ff1->PREC(:,0,0)
delete(ff1)

print("read "+Lbfile1)
ff1=addfile(Lbfile1,"r")
lem_tas_clim=ff1->TS(0,0,0)
delete(ff1)

print("read "+Lbfile2)
ff1=addfile(Lbfile2,"r")
lem_pr_clim=ff1->PRECC(0,0,0)
delete(ff1)

;make percentages PDRMIP
;===================
sol_pr_scaled=sol_pr
sol_pr_scaled=(/100*sol_pr/pr_clim/)

co2_pr_scaled=co2_pr
co2_pr_scaled=(/100*co2_pr/pr_clim/)

; play with CES-LEM
;===================
lem_co2_pr_scaled=(/100*(lem_co2_pr-lem_pr_clim)/lem_pr_clim/)
lem_sol_pr_scaled=(/100*(lem_sol_pr-lem_pr_clim)/lem_pr_clim/)

lem_sol_tas=(/lem_sol_tas-lem_tas_clim/)
lem_co2_tas=(/lem_co2_tas-lem_tas_clim/)

lem_sol_tas_unf=(/lem_sol_tas_unf-lem_tas_clim/)
lem_co2_tas_unf=(/lem_co2_tas_unf-lem_tas_clim/)

neqv=equiv_sample_size(tsi,0.01,0)
ttpr     = rtest(0.6, 8, 0) 
print(ttpr)
;print(lem_sol_tas)
print (neqv)


nn1=0
nn2=sol_nsteps-1
cor_tas=esccr(tsi(nn1:nn2),lem_sol_tas(nn1:nn2),50)
;cor_tas=esccr(lem_sol_tas(nn1:nn2),lem_sol_tas(nn1:nn2),50)
print(cor_tas)
neqv=equiv_sample_size(tsi,0.01,0)
ttpr     = rtest(0.6, 8, 0) 
print(ttpr)


cor_pr=esccr(tsi(nn1:nn2),lem_sol_pr_scaled(nn1:nn2),50)
cor_prtas=esccr(lem_sol_pr(nn1:nn2),lem_sol_tas(nn1:nn2),50)

nn2=sol_nsteps-1
;make regressions
lem_sol_pr_rg=regline(tsi(49:nn2),lem_sol_pr_scaled(49:nn2))
lem_sol_tas_rg=regline(tsi(49:nn2),lem_sol_tas(49:nn2))

;lem_co2_rg=regline(lem_co2_tas,lem_co2_pr_scaled)
nlagr = 30 ; 0,1,2
rlag = new (nlagr, "double")
rlag_pr = new (nlagr, "double")
Nend = nn2-1
do n=0,nlagr-1
  rlag(n) = regline(tsi(49:Nend-n), lem_sol_tas(49+n:Nend) )
  rlag_pr(n) = regline(tsi(49:Nend-n), lem_sol_pr_scaled(49+n:Nend) )
end do
print(rlag) 
print(rlag_pr) 

     
print(lem_sol_pr_rg)
;print(lem_sol_tas_rg)
;exit


print (cor_tas)
;print (cor_pr)
;print (cor_prtas)
;exit

print("plot "+graph1+".png")
wks = gsn_open_wks("png",graph1+"dd")
plot = new(2,graphic)

res                     = True                   ; plot mods desired
res@vpWidthF = 0.7
res@vpHeightF = 0.4
res@gsnDraw            =  False                   ; do not draw the plot
res@gsnFrame           =  False                   ; do not advance the frame
res@xyLineThicknesses      = 10    
res@xyLineColor      = "red"
res@tiYAxisString = "Temperature [K]"
res@tiYAxisFontColor = res@xyLineColor
 
resR                        = True                      
resR@xyDashPatterns         = 1                   ; dashed line for 2nd
resR@xyLineThicknesses      = 10                   ; thicker line
resR@xyLineColor      = "blue"
resR@tiYAxisString = "Precipitation [%]"
resR@tiYAxisFontColor = resR@xyLineColor

res@trYMinF=0.
res@trYMaxF=3.
resR@trYMinF=-1.
resR@trYMaxF=7.

; for PDRMIP
;plot(0) = gsn_csm_xy2(wks,time_pdrmip,sol_tas,sol_pr_scaled,res,resR)
;plot(1) = gsn_csm_xy2(wks,time_pdrmip,co2_tas,co2_pr_scaled,res,resR)

res@trYMinF=-0.25
res@trYMaxF=.25
resR@trYMinF=-.5
resR@trYMaxF=.5
res@trXMinF=850
res@trXMaxF=2005
dd1_data=new((/2,dimsizes(time_sol)/),float)
dd1_data(0,:)=dble2flt((/lem_sol_tas/))
dd1_data(1,:)=lem_sol_tas_unf
;plot(0) = gsn_csm_xy2(wks,time_sol,lem_sol_tas,lem_sol_pr_scaled,res,resR)
plot(0) = gsn_csm_xy2(wks,time_sol,dd1_data,lem_sol_pr_scaled,res,resR)
delete(dd1_data)
res@xyLineColors      = "grey"
res@xyLineThicknesses      = 10                   ; thicker line

tsi=tsi/10.-0.135
 plot1= gsn_csm_xy(wks,time_sol,tsi,res)
 overlay(plot(0),plot1)
 res@trYMinF=-0.2
res@trYMaxF=1.5
resR@trYMinF=-.5
resR@trYMaxF=2

  res@trXMinF=850
res@trXMaxF=2005
;      resR@pmLegendDisplayMode    = "Always"
;  resR@pmLegendSide           = "Top"               ; Change location of 
;  resR@pmLegendParallelPosF   = .70                 ; move units right
;  resR@pmLegendOrthogonalPosF = -0.9                ; move units down

;  resR@pmLegendWidthF         = 0.12                ; Change width and
;  resR@pmLegendHeightF        = 0.05                ; height of legend.
;  resR@lgLabelFontHeightF     = .018                 ; change font height
;  resR@lgPerimOn              = False               ; no box around
  res@xyLineColors      = "red"
  resR@xyLineColor      = "blue"
 ; resR@xyExplicitLegendLabels = (/" Precipitation","Temperature"/)  
 dd1_data=new((/2,dimsizes(time_co2)/),float)
dd1_data(0,:)=dble2flt((/lem_co2_tas/))
dd1_data(1,:)=lem_co2_tas_unf

 ;plot(1) = gsn_csm_xy2(wks,time_co2,lem_co2_tas,lem_co2_pr_scaled,res,resR)
  plot(1) = gsn_csm_xy2(wks,time_co2,dd1_data,lem_co2_pr_scaled,res,resR)

 res@xyLineColors      = "grey"
res@xyLineThicknesses      = 10                   ; thicker line
co2=co2/100.-0.07
 plot2= gsn_csm_xy(wks,time_co2,co2,res)
 overlay(plot(1),plot2)
  
 txres               = True                            
 txres@txFontHeightF = 0.02             ; Set the font height
text = gsn_add_text(wks,plot(0),"CESM-SOL",1800.,0.2 ,txres)
text2 = gsn_add_text(wks,plot(0),"150-yr filtered",1800.,0.17 ,txres)
text1 = gsn_add_text(wks,plot(1),"CESM-GHG ",1840.,1.3 ,txres)
text3 = gsn_add_text(wks,plot(1),"30-yr filtered",1840.,1.18 ,txres)
resP                     = True         ; modify the panel plot
resP@gsnMaximize         = True         ; large format
resP@gsnPanelFigureStrings= (/"a)","b)"/)
resP@gsnPanelFigureStringsFontHeightF=0.01
resP@amJust     = "TopLeft"


gsn_panel(wks,plot,(/1,2/),resP)
system("convert -trim "+graph1+"dd.png"+" "+graph1+".png")
system("rm "+graph1+"dd.png")

delete(plot)
delete(wks)

end